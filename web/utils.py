import json
import re
import subprocess
import sys
from subprocess import SubprocessError
from typing import Literal, Optional, Sequence, TypeVar

from flask import Flask, render_template, request
from pydantic import BaseModel, ValidationInfo, field_validator
from returns.converters import maybe_to_result
from returns.curry import partial
from returns.io import IOFailure, IOResult, IOResultE, IOSuccess
from returns.iterables import Fold
from returns.maybe import Maybe, Nothing, Some, maybe
from returns.pipeline import flow
from returns.pointfree import bind, cond, map_
from returns.result import Failure, Result, ResultE, Success, safe


@maybe
def get_request(key: str) -> Optional[str]:
    """
    Get a value from the query string of a Flask request.
    """
    return request.args.get(key)


## ----------------------------------------------------------------------------
## Genomic region input validation
## ----------------------------------------------------------------------------

# basic chromosome name variants for human reference genomes
accepted_chromosome_names = (
    {"chr" + str(i) for i in range(1, 23)}
    | {str(i) for i in range(1, 23)}
    | {"chrX", "chrY", "chrM", "X", "Y", "M"}
)


class GenomicRegion(BaseModel):
    """
    A class to represent a genomic region, with basic validation
    of the coordinates and chromosome name.
    """

    chrom: str
    start: int
    end: int

    @field_validator("start", "end")
    def validate_coordinates(cls, value: int):
        if value < 0:
            raise ValueError("Coordinates must be positive integers")
        return value

    @field_validator("end")
    def validate_end(cls, value: int, info: ValidationInfo):
        if value < info.data["start"]:
            raise ValueError(
                "End coordinate must be greater than or equal to start coordinate"
            )
        return value

    @field_validator("chrom")
    def validate_chrom(cls, value: str):
        if value not in accepted_chromosome_names:
            raise ValueError(f"Invalid chromosome name: {value}")
        return value

    @classmethod
    def from_str(cls, region: str) -> "GenomicRegion":
        """parse from genomic region string with format chrom:start-end"""
        chrom, start, end = re.split(r"[:-]", region)
        return cls(chrom=chrom, start=int(start), end=int(end))

    def __str__(self):
        return f"{self.chrom}:{self.start}-{self.end}"

    def __repr__(self):
        return f"GenomicRegion(chrom={self.chrom}, start={self.start}, end={self.end})"

    @safe
    @staticmethod
    def safe_parse(region: str) -> "GenomicRegion":
        return GenomicRegion.from_str(region)


## ----------------------------------------------------------------------------
## SVType input validation
## ----------------------------------------------------------------------------
class SVType(BaseModel):
    value: str

    @field_validator("value")
    def validate_sv_type(cls, value: str):
        if value not in ["DEL", "DUP", "INV", "BND", "INS"]:
            raise ValueError(f"Invalid SV type: {value}")
        return value

    def __str__(self):
        return self.value

    @safe
    @staticmethod
    def safe_parse(value: str) -> "SVType":
        return SVType(value=value)


## ----------------------------------------------------------------------------
## slop input validation: overkill...
## ----------------------------------------------------------------------------
class Slop(BaseModel):
    value: int

    @field_validator("value")
    def validate_slop(cls, value: int):
        if value < 0:
            raise ValueError(f"Slop must be a non-negative integer")
        return value

    def __str__(self):
        return str(self.value)

    @safe
    @staticmethod
    def safe_parse(value: int) -> "Slop":
        return Slop(value=value)


## ----------------------------------------------------------------------------
## Safe retrieval of validated inputs
## ----------------------------------------------------------------------------
def get_genomic_region(side: Literal["left", "right"]) -> ResultE[GenomicRegion]:
    """
    Retrieve genomic region from request and convert to validated GenomicRegion.
    """
    return flow(
        side,
        get_request,
        bind(GenomicRegion.safe_parse),
    )


def get_sv_type() -> ResultE[SVType]:
    """
    Retrieve SV type from request and convert to validated SVType.
    """
    return flow(
        "sv_type",
        get_request,
        bind(SVType.safe_parse),
    )


def get_slop() -> ResultE[int]:
    """
    Retrieve slop from request and convert to int.
    """
    return flow(
        "slop",
        get_request,
        bind(Slop.safe_parse),
    )


def get_raw() -> bool:
    """
    Retrieve raw from request. Treat *any* presence of raw as True.
    """
    match get_request("raw"):
        case Some(_):
            return True
        case _:
            return False


## ----------------------------------------------------------------------------
## Construct stix query commands
## ----------------------------------------------------------------------------
def construct_stix_query_command(
    stix_path: str,
    stix_index: str,
    stix_db: str,
    left_region: ResultE[GenomicRegion],
    right_region: ResultE[GenomicRegion],
    svtype: ResultE[SVType],
    slop: ResultE[int] = Success(100),
) -> ResultE[str]:
    """
    Construct a stix query command from validated inputs.
    Returns Failure if any of the inputs are invalid.
    """

    # get value of slop or Failure with default value
    sl = slop.value_or(100)

    return Result.do(
        f"{stix_path} -i {stix_index} -d {stix_db} -l {left} -r {right} -t {svt} -s {sl} -j"
        for left in left_region
        for right in right_region
        for svt in svtype
    )


def construct_sharded_commands(
    stix_path: str,
    index_shards: list[str],
    db_shards: list[str],
    left_region: ResultE[GenomicRegion],
    right_region: ResultE[GenomicRegion],
    svtype: ResultE[SVType],
    slop: ResultE[int] = Success(100),
) -> list[ResultE[str]]:
    """
    Construct a list of stix query commands for each shard from validated inputs.
    Returns Failure if any of the inputs are invalid
    """

    return [
        construct_stix_query_command(
            stix_path=stix_path,
            left_region=left_region,
            right_region=right_region,
            svtype=svtype,
            slop=slop,
            stix_index=index_shard,
            stix_db=db_shard,
        )
        for index_shard, db_shard in zip(index_shards, db_shards)
    ]


@safe
def stix_query(cmd: str) -> bytes:
    """
    Run a stix query command and return the output.
    """
    return subprocess.Popen(
        cmd.split(), stdout=subprocess.PIPE, stderr=sys.stderr
    ).communicate()[0]


def sharded_stix_queries(
    cmds: list[str],
) -> list[ResultE[bytes]]:
    """
    Run a list of stix query commands and return the output.
    """
    return [stix_query(cmd) for cmd in cmds]


def parse_result(result: bytes | str) -> list[dict]:
    """
    Parse the JSON output of a stix query.
    Each stix query output is a JSON string with fields:
        {
            "results": {
                "summary": [
                    {
                        ... # I'm going to ignore the summary
                    }
                ],
                "samples": [
                    {
                        "Giggle_File_Id": 0,
                        "Sample": "HG00118",
                        "Alt_File": "HG00118.bed.gz",
                        "Pairend": "0",
                        "Split": "0"
                    }
                ]
            }
        }

    We are only concerned with the list of samples, so for each sharded query output,
    just extract the list of samples and append to a single list.
    """
    return json.loads(result)["results"]["samples"]


def merge_sharded_results(results: Sequence[ResultE[bytes | str]]) -> ResultE[str]:
    """
    Merges results of sharded stix queries into single output JSON string.
    Returns JSON string where "samples" field is the merged list (summary discarded).
    """
    return (
        Fold.collect(results, Success(()))  # ResultE[tuple[str]]
        .map(lambda results: [parse_result(r) for r in results])  # type: ignore
        .map(lambda results: [item for sublist in results for item in sublist])
        .map(lambda samples: json.dumps({"results": {"samples": samples}}))
    )


## ----------------------------------------------------------------------------
## Shown when raw data is not requested
## ----------------------------------------------------------------------------
def render_sv_view(
    left: ResultE[GenomicRegion],
    right: ResultE[GenomicRegion],
    svtype: ResultE[SVType],
):
    """
    Render a view of the SV query results. If the input data is invalid,
    replace with a default value.
    """
    data = Result.value_or(
        Result.do(
            {"svtype": str(svt), "left": str(l), "right": str(r)}
            for svt in svtype
            for l in left
            for r in right
        ),
        {
            "svtype": "DEL",
            "left": "10:105053143-105053143",
            "right": "10-105054173-105054173",
        },
    )
    return render_template("sv.html", data=data)
