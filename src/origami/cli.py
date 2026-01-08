import argparse
import os
from origami.workflows import run_full, run_trim_only, run_oric_only
from origami.logging_setup import setup_logging
from importlib.resources import files

# Default DB paths (relative to project root)
DNAA_PREFIX = "dnaA_clustered100_DB"
ORIC_PREFIX = "OriC_clustered100_DB"

def get_db_prefixes():
    db_root = files("origami.data.clustered_DB")
    dnaA_prefix = db_root / "dnaA_clustered100_DB"
    oric_prefix = db_root / "OriC_clustered100_DB"
    return str(dnaA_prefix), str(oric_prefix)


def main():
    parser = argparse.ArgumentParser(
        prog="origenomi",
        description="Trim circular overlaps, find dnaA/oriC, rotate, and report."
    )
    subparsers = parser.add_subparsers(dest="command")

    def add_common(p):
        p.add_argument("-i", "--input", required=True, help="Input FASTA (single or multi-FASTA)")
        p.add_argument("-o","--out-dir", default=".", help="Output directory (default: .)")
        p.add_argument("--keep-temp", action="store_true", help="Keep temp files")
        p.add_argument("--verbose", action="store_true", help="Verbose logging")
        p.add_argument("--keep-plot", action="store_true", help="Keep the plots")

    # Full pipeline
    p_run = subparsers.add_parser("run", help="Trim + dnaA/oriC + rotate + report")
    add_common(p_run)

    # Trim-only
    p_trim = subparsers.add_parser("trim", help="Trim duplicated terminal overlaps only")
    add_common(p_trim)

    # OriC-only
    p_oric = subparsers.add_parser("oric", help="dnaA/oriC + rotate on original input")
    add_common(p_oric)

    args = parser.parse_args()

    # Inject default DB paths automatically
    args.dnaA_db, args.oric_db = get_db_prefixes()

    args.out_dir = os.path.abspath(args.out_dir)
    os.makedirs(args.out_dir, exist_ok=True)

    logger = setup_logging(args.verbose)

    if args.command == "run":
        run_full(args, logger)
    elif args.command == "trim":
        run_trim_only(args, logger)
    elif args.command == "oric":
        run_oric_only(args, logger)
    else:
        parser.print_help()
