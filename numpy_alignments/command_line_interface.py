import logging
logging.basicConfig(level=logging.INFO)
import argparse
from .numpy_alignments import NumpyAlignments
from .comparer import Comparer
import sys


def main():
    run_argument_parser(sys.argv[1:])


def store_alignments(args):
    if args.type == "pos":
        a = NumpyAlignments.from_pos(args.n_alignments)
    elif args.type == "sam":
        a = NumpyAlignments.from_sam(args.n_alignments)
    elif args.type == "truth":
        a = NumpyAlignments.from_truth(args.n_alignments)
    else:
        logging.error("Invalig type %s" % args.type)
        sys.exit()

    a.to_file(args.file_name)


def get_correct_rates(args):
    logging.info("Reading alignments from file")
    truth_alignments = NumpyAlignments.from_file(args.truth_alignments + ".npz")
    compare_alignments = {c: NumpyAlignments.from_file(c + ".npz") for c in args.compare_alignments.split(",")}

    logging.info("Comparing")
    comparer = Comparer(truth_alignments, compare_alignments)
    rates = comparer.get_correct_rates()
    for name, rate in rates.items():
        print(name, rate)


def compare_alignments(args):
    logging.info("Reading alignments from file")
    truth_alignments = NumpyAlignments.from_file(args.truth_alignments + ".npz")
    compare_alignments = {c: NumpyAlignments.from_file(c + ".npz") for c in args.compare_alignments.split(",")}

    logging.info("Comparing")
    for type in ["all", "variants", "nonvariants"]:
        comparer = Comparer(truth_alignments, compare_alignments, type=type)
        comparer.create_roc_plots(save_to_file=args.save_to_file + "_" + type + ".png")

    #comparer.get_wrong_alignments_correct_by_other("two_step_approach", "vg_chr20")
    #comparer.get_wrong_alignments_correct_by_other("bwa_10m_tuned", "vg_10m")



def run_argument_parser(args):
    parser = argparse.ArgumentParser(
        description='Numpy Alignments',
        prog='numpy_alignments',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=50, width=100))

    subparsers = parser.add_subparsers()

    # Store alignments
    store = subparsers.add_parser("store")
    store.add_argument("type", help="Type of alignments. Either sam, pos or truth.")
    store.add_argument("file_name", help="File name to store alignments to")
    store.add_argument("n_alignments", help="Must be >= number of alignments that is expected", type=int)
    store.set_defaults(func=store_alignments)

    # Compare alignments
    compare = subparsers.add_parser("compare")
    compare.add_argument("truth_alignments")
    compare.add_argument("compare_alignments", help="Comma-separated list of files to compare")
    compare.add_argument("-f", "--save-to-file", help="File name to save figure to (jpg)")
    compare.set_defaults(func=compare_alignments)

    # Compare (get correct rates)
    compare = subparsers.add_parser("get_correct_rates")
    compare.add_argument("truth_alignments")
    compare.add_argument("compare_alignments", help="Comma-separated list of files to compare")
    compare.set_defaults(func=get_correct_rates)

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)
    args.func(args)
