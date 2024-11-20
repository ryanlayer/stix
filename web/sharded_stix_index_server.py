import argparse
import subprocess
import sys
from typing import Optional

from flask import Flask, render_template, request
from returns.maybe import Maybe, Nothing, Some, maybe
from returns.pipeline import flow
from returns.result import Failure, Result, Success, safe

from utils import (
    construct_sharded_commands,
    get_genomic_region,
    get_raw,
    get_slop,
    get_sv_type,
    merge_sharded_results,
    render_sv_view,
    stix_query,
)

app = Flask(__name__)


## Args for starting the web app ----------------------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="STIX web server.")
    # TODO
    parser.add_argument(
        "--stix", dest="stix_path", help="Path to stix", type=str, required=True
    )
    parser.add_argument(
        "--port", dest="port", help="Port to run on", type=int, default=5000
    )
    parser.add_argument(
        "--db_shards",
        dest="db_shards",
        help="sharded STIX ped db's",
        narg="+",
        type=list[str],
        required=True,
    )
    parser.add_argument(
        "--index_shards",
        dest="index_shards",
        help="sharded STIX indices",
        narg="+",
        type=list[str],
        required=True,
    )
    return parser.parse_args()


## App entry point ------------------------------------------------------------
@app.route("/test")
def main():
    ## stix index info
    stix_path = app.config["stix_path"]
    index_shards = app.config["index_shards"]
    db_shards = app.config["db_shards"]

    ## stix query params
    left_bp = get_genomic_region("left")
    right_bp = get_genomic_region("right")
    svtype = get_sv_type()
    slop = get_slop()
    raw = get_raw()

    ## query the sharded index and return as JSON
    if raw:
        stix_shard_cmds = construct_sharded_commands(
            stix_path=stix_path,
            index_shards=index_shards,
            db_shards=db_shards,
            left_region=left_bp,
            right_region=right_bp,
            svtype=svtype,
            slop=slop,
        )
        sharded_query_results = [cmd.bind(stix_query) for cmd in stix_shard_cmds]
        return merge_sharded_results(sharded_query_results).value_or("{}")

    ## OR pass the args to the render template
    else:
        return render_sv_view(left=left_bp, right=right_bp, svtype=svtype)


if __name__ == "__main__":

    args = parse_args()
    app.config["stix_path"] = args.stix_path
    app.config["index_shards"] = args.indices
    app.config["db_shards"] = args.db_shards
    print(args.db_path)
    app.run(host="0.0.0.0", port=args.port)
