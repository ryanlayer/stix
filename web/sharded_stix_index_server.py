import argparse
from pprint import pprint

from flask import Flask

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
        "--host", dest="host", help="IP to run on", type=str, default="127.0.0.1"
    )
    parser.add_argument(
        "--db_shards",
        dest="db_shards",
        help="sharded STIX ped db's",
        nargs="+",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--index_shards",
        dest="index_shards",
        help="sharded STIX indices",
        nargs="+",
        type=str,
        required=True,
    )
    return parser.parse_args()


## App entry point ------------------------------------------------------------
@app.route("/")
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
    else:
        return "only raw queries supported"

    ## OR pass the args to the render template
    # else:
    #     return render_sv_view(left=left_bp, right=right_bp, svtype=svtype)


if __name__ == "__main__":
    args = parse_args()
    app.config["stix_path"] = args.stix_path
    app.config["index_shards"] = args.index_shards
    app.config["db_shards"] = args.db_shards

    pprint(args.db_shards)
    pprint(args.index_shards)

    app.run(host=args.host, port=args.port)
