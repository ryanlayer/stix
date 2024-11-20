import sys
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pprint import pprint
from unittest.mock import patch

from returns.maybe import Some
from returns.result import Failure, Success, safe

from utils import *


class TestUtils(unittest.TestCase):

    def setup(self):
        self.held, sys.stdout = sys.stdout, StringIO()

    def teardown(self):
        self.output = sys.stdout.getvalue()
        sys.stdout = self.held
        print(self.output)

    ## ------------------------------------------------------------------------
    ## get function tests
    ## ------------------------------------------------------------------------
    @patch("utils.get_request")
    def test_get_genomic_region1(self, mock_request):
        ## success case
        mock_request.return_value = Some("chr1:100-200")
        result = get_genomic_region("left")

        match result:
            case Success(region):
                self.assertEqual(region.chrom, "chr1")
                self.assertEqual(region.start, 100)
                self.assertEqual(region.end, 200)
            case Failure(error):
                self.fail(
                    f"Expected Success with input of 'chr1:100-200' but got\n: {error}"
                )
            case _:
                self.fail(f"unexpected case: {result}")

    @patch("utils.get_request")
    def test_get_genomic_region2(self, mock_request):
        ## fail case: coordinates not formatted correctly
        mock_request.return_value = Some("chr1:100-200-300")
        result = get_genomic_region("left")
        self.assertIsInstance(result, Failure)

    @patch("utils.get_request")
    def test_get_genomic_region3(self, mock_request):
        ## fail case: non-accepted chromosome
        mock_request.return_value = Some("hellooo:100-200")
        result = get_genomic_region("right")
        self.assertIsInstance(result, Failure)

    @patch("utils.get_request")
    def test_get_genomic_region4(self, mock_request):
        ## fail case: end < start
        mock_request.return_value = Some("1:100-10")
        result = get_genomic_region("left")
        self.assertIsInstance(result, Failure)

    ## ------------------------------------------------------------------------
    ## construct_sharded_commands tests
    ## ------------------------------------------------------------------------
    # success case
    def test_construct_sharded_commands1(self):
        stix_path = "stix"
        index_shards = ["index1", "index2"]
        db_shards = ["db1", "db2"]
        svtype = SVType.safe_parse("DEL")
        left_bp = GenomicRegion.safe_parse("chr1:100-200")
        right_bp = GenomicRegion.safe_parse("chr1:300-400")

        result = construct_sharded_commands(
            stix_path=stix_path,
            db_shards=db_shards,
            index_shards=index_shards,
            left_region=left_bp,
            right_region=right_bp,
            svtype=svtype,
        )
        actual = [
            Success(
                "stix -i index1 -d db1 -l chr1:100-200 -r chr1:300-400 -t DEL -s 100 -j"
            ),
            Success(
                "stix -i index2 -d db2 -l chr1:100-200 -r chr1:300-400 -t DEL -s 100 -j"
            ),
        ]
        self.assertEqual(result, actual, f"\n\n{result=}\n\n{actual=}")

    # fail case: invalid svtype
    def test_construct_sharded_commands2(self):
        stix_path = "stix"
        index_shards = ["index1", "index2"]
        db_shards = ["db1", "db2"]
        svtype = SVType.safe_parse("HELLO")
        left_bp = GenomicRegion.safe_parse("chr1:100-200")
        right_bp = GenomicRegion.safe_parse("chr1:300-400")

        result = construct_sharded_commands(
            stix_path=stix_path,
            db_shards=db_shards,
            index_shards=index_shards,
            left_region=left_bp,
            right_region=right_bp,
            svtype=svtype,
        )
        for r in result:
            self.assertIsInstance(r, Failure)

    # fail case: invalid genomic region
    def test_construct_sharded_commands3(self):
        stix_path = "stix"
        index_shards = ["index1", "index2"]
        db_shards = ["db1", "db2"]
        svtype = SVType.safe_parse("DEL")
        left_bp = GenomicRegion.safe_parse("chr1:100-200")
        right_bp = GenomicRegion.safe_parse("chr1\t300\t400")

        result = construct_sharded_commands(
            stix_path=stix_path,
            db_shards=db_shards,
            index_shards=index_shards,
            left_region=left_bp,
            right_region=right_bp,
            svtype=svtype,
        )
        for r in result:
            self.assertIsInstance(r, Failure)

    ## ------------------------------------------------------------------------
    ## parse_result tests
    ## ------------------------------------------------------------------------
    def test_parse_result1(self):
        input = """{
            "results": {
                "summary": [ { } ],
                "samples": [
                    {
                        "Giggle_File_Id": "0",
                        "Sample": "HG00118",
                        "Alt_File": "HG00118.bed.gz",
                        "Pairend": "0",
                        "Split": "0"
                    }
                ]
            }
        }"""

        try:
            result = Success(parse_result(input))
        except Exception as e:
            result = Failure(e)

        expected = Success(
            [
                {
                    "Giggle_File_Id": "0",
                    "Sample": "HG00118",
                    "Alt_File": "HG00118.bed.gz",
                    "Pairend": "0",
                    "Split": "0",
                }
            ]
        )
        self.assertEqual(result, expected)

    def test_parse_result2(self):
        input = """{
            "results": {
                "summary": [ { } ],
                "samples": [
                    {
                        "Giggle_File_Id": "0",
                        "Sample": "HG00118",
                        "Alt_File": "HG00118.bed.gz",
                        "Pairend": "0",
                        "Split": "0"
                    }
                ], # trailing comma should make parsing fail
            }
        }"""
        try:
            result = Success(parse_result(input))
        except Exception as e:
            result = Failure(e)
        self.assertIsInstance(result, Failure)

    ## ------------------------------------------------------------------------
    ## merged_sharded_results tests
    ## ------------------------------------------------------------------------
    def test_merge_sharded_results(self):
        input = [
            Success(
                """{
            "results": {
                "summary": [ { } ],
                "samples": [
                    {
                        "Giggle_File_Id": "0",
                        "Sample": "HG00118",
                        "Alt_File": "HG00118.bed.gz",
                        "Pairend": "0",
                        "Split": "0"
                    }
                ]
            }
        }"""
            ),
            Success(
                """{
            "results": {
                "summary": [ { } ],
                "samples": [
                    {
                        "Giggle_File_Id": "1",
                        "Sample": "HG00119",
                        "Alt_File": "HG00119.bed.gz",
                        "Pairend": "0",
                        "Split": "0"
                    }
                ]
            }
        }"""
            ),
        ]
        result = merge_sharded_results(input)  # type: ignore
        expected = Success(
            json.dumps(  # bleh
                {
                    "results": {
                        "samples": [
                            {
                                "Giggle_File_Id": "0",
                                "Sample": "HG00118",
                                "Alt_File": "HG00118.bed.gz",
                                "Pairend": "0",
                                "Split": "0",
                            },
                            {
                                "Giggle_File_Id": "1",
                                "Sample": "HG00119",
                                "Alt_File": "HG00119.bed.gz",
                                "Pairend": "0",
                                "Split": "0",
                            },
                        ]
                    }
                }
            )
        )
        self.assertEqual(result, expected)


if __name__ == "__main__":
    unittest.main()
