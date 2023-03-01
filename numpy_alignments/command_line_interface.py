import logging
logging.basicConfig(level=logging.INFO)
import argparse
from .numpy_alignments import NumpyAlignments, NumpyAlignments2
from .comparer import Comparer
import sys
from .htmlreport import make_report


def main():
    run_argument_parser(sys.argv[1:])


def set_correctness(args):
    logging.info("Reading alignments")
    truth = NumpyAlignments.from_file(args.truth_alignments)
    alignments = NumpyAlignments.from_file(args.alignments)
    logging.info("Setting correctness")
    alignments.is_correct = None
    alignments.set_correctness(truth, force=True)
    alignments.to_file(args.alignments)
    logging.info("Correctness was set and file written to same file again")


def make_html_report_wrapper(args):
    from time import time
    import os


    if args.report_id is not None:
        report_id = args.report_id
    else:
        report_id = str(time()).split(".")[0]

    if not os.path.isdir(report_id):
        os.mkdir(report_id)

    logging.info("Report will be in directory %s" % report_id)

    # Make
    truth_alignments = NumpyAlignments.from_file(args.truth_alignments)
    ids = args.compare_alignments.split(",")
    compare_alignments = {c: NumpyAlignments.from_file(c) for c in args.compare_alignments.split(",")}

    if args.names is not None:
        names = args.names.split(",")
        names = {name: names[i] for i, name in enumerate(ids)}
    else:
        names = {name: name for name in ids}

    colors = args.colors.split(",")
    colors = {name: colors[i] for i, name in enumerate(ids)}

    if len(colors) != len(names) or len(names) != len(compare_alignments):
        logging.error("Not enough names/colors/alignments. Numbers do not match")
        sys.exit()

    logging.info("Making plots")
    for type in ["all", "variants", "nonvariants"]:
        comparer = Comparer(truth_alignments, compare_alignments, colors, type=type)
        comparer.create_roc_plots(save_to_file=report_id + "/" + type + ".html")

    html = make_report(report_id, ids, names, colors)
    with open(report_id + "/report.html", "w") as f:
        f.write(html)
    logging.info("Final report written to %s/report.html" % report_id)


def store_alignments(args):
    if args.type == "pos":
        a = NumpyAlignments.from_pos(args.n_alignments)
    elif args.type == "sam":
        a = NumpyAlignments.from_sam(args.n_alignments)
    elif args.type == "truth":
        a = NumpyAlignments.from_truth(args.n_alignments)
    elif args.type == "bed":
        a = NumpyAlignments.from_bed(args.n_alignments)
    elif args.type == "vgpos":
        a = NumpyAlignments.from_vgpos(args.n_alignments)
    elif args.type == "bam":
        if args.n_variants is not None:
            a = NumpyAlignments2.from_bam_and_nvariants_txt(args.input, args.n_variants)
        else:
            a = NumpyAlignments2.from_bam(args.input)
    else:
        logging.error("Invalid type %s" % args.type)
        sys.exit()

    a.to_file(args.file_name)


def get_correct_rates(args):
    logging.info("Reading alignments from file")
    try:
        truth_alignments = NumpyAlignments.from_file(args.truth_alignments)
    except KeyError:
        truth_alignments = NumpyAlignments2.from_file(args.truth_alignments)

    try:
        compare_alignments = {c: NumpyAlignments.from_file(c) for c in args.compare_alignments.split(",")}
    except KeyError:
        compare_alignments = {c: NumpyAlignments2.from_file(c) for c in args.compare_alignments.split(",")}

    type = args.type #edit
    
    logging.info("Comparing..")
    comparer = Comparer(truth_alignments, compare_alignments, type=type, allowed_mismatch=args.allowed_bp_mismatch) #edit
    rates = comparer.get_correct_rates(args.min_mapq)
    for name, rate in rates.items():
        recall = rate[0]
        one_minus_precision = rate[1]
        precision = 1 - one_minus_precision

        if args.report_type == "all":
            print(name, rate[0], rate[1])
        elif args.report_type == "recall":
            print(recall)
        elif args.report_type == "one_minus_precision":
            print(one_minus_precision)
        elif args.report_type == "f1_score":
            f1 = 2 * precision * recall / (precision + recall)
            print(f1)
        else:
            raise Exception("Invalid report type")

def get_correct_rates_multi(args):
    truth_alignments = NumpyAlignments.from_file(args.truth_alignments)

    n_correct = 0
    with open(args.compare_alignments) as f:
        reads_checked = set()
        for line in f:
            l = line.split()
            name = l[3]
            pos = int(l[1])

            correct_pos = truth_alignments.positions[int(name)]
            if abs(correct_pos-pos) < args.allowed_bp_mismatch and name not in reads_checked:
                n_correct += 1
                reads_checked.add(name)

    logging.info("N correct: %d" % n_correct)
    logging.info("Rate: %.3f" % (n_correct / len(truth_alignments.positions)))

def compare_alignments(args):
    logging.info("Reading alignments from file")
    truth_alignments = NumpyAlignments.from_file(args.truth_alignments)
    compare_alignments = {c: NumpyAlignments.from_file(c) for c in args.compare_alignments.split(",")}

    logging.info("Comparing")
    for type in ["all", "variants", "nonvariants"]:
        comparer = Comparer(truth_alignments, compare_alignments, type=type, allowed_mismatch=args.allowed_mismatch)
        save_to_file = None
        if args.save_to_file is not None:
            save_to_file = args.save_to_file + "_" + type + ".html"
        comparer.create_roc_plots(save_to_file=save_to_file, limit_comparison=args.limit_to_n_reads)

    #comparer.get_wrong_alignments_correct_by_other("two_step_approach", "vg_chr20")
    #comparer.get_wrong_alignments_correct_by_other("bwa_10m_tuned", "vg_10m")


def rename(args):
    assert args.posfile is not None or args.fq is not None, "Either --fq or --posfile must be specified"

    if args.posfile is not None:
        if args.posfile == "-":
            f = sys.stdin
        else:
            f = open(args.posfile)

        for i, line in enumerate(f):
            l = line.split()
            new_id = "%09d" % i
            l[0] = new_id
            print('\t'.join(l))
    elif args.fq is not None:
        if args.fq == "-":
            f = sys.stdin
        else:
            f = open(args.fq)

        for i, line in enumerate(f):
            if i % 4 != 0:
                print(line.strip())
            else:
                print("@%09d" % (i//4))

    logging.info("Done")



def run_argument_parser(args):
    parser = argparse.ArgumentParser(
        description='Numpy Alignments',
        prog='numpy_alignments',
        formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=50, width=100))

    subparsers = parser.add_subparsers()

    # Store alignments
    store = subparsers.add_parser("store")
    store.add_argument("-c", "--coordinate-map", required=False, help="If set, can use coordinate map to count variants (only supported for SAM-files)")
    store.add_argument("-i", "--input", required=False)
    store.add_argument("-n", "--n_variants", required=False)
    store.add_argument("type", help="Type of alignments. Either sam, pos or truth.")
    store.add_argument("file_name", help="File name to store alignments to")
    store.add_argument("n_alignments", help="Must be >= number of alignments that is expected", type=int)
    store.set_defaults(func=store_alignments)

    # Compare alignments
    compare = subparsers.add_parser("compare")
    compare.add_argument("truth_alignments")
    compare.add_argument("compare_alignments", help="Comma-separated list of files to compare")
    compare.add_argument("-f", "--save-to-file", help="File name to save figure to (jpg)")
    compare.add_argument("-m", "--allowed-mismatch", help="Maximum number of bp between read position and correct position in order for read to considered as correctly mapped.", type=int, default=150)
    compare.add_argument("-l", "--limit-to-n-reads", help="Limit comparison to max this number of reads in order to make things faster", required=False, type=int, default=None)
    compare.set_defaults(func=compare_alignments)

    # Compare (get correct rates)
    compare = subparsers.add_parser("get_correct_rates")
    compare.add_argument("truth_alignments")
    compare.add_argument("compare_alignments", help="Comma-separated list of files to compare")
    compare.add_argument("type")
    compare.add_argument("-m", "--min-mapq", type=int, default=0)
    compare.add_argument("-t", "--allowed-bp-mismatch", type=int, default=150)
    compare.add_argument("-r", "--report-type", default="all", help="all, recall, one_minus_precision, f1_score")
    compare.set_defaults(func=get_correct_rates)

    #
    compare = subparsers.add_parser("get_correct_rates_multi")
    compare.add_argument("truth_alignments")
    compare.add_argument("compare_alignments")
    compare.add_argument("type")
    compare.add_argument("-m", "--min-mapq", type=int, default=0)
    compare.add_argument("-t", "--allowed-bp-mismatch", type=int, default=150)
    compare.set_defaults(func=get_correct_rates_multi)

    # Make ROC html report
    cmd = subparsers.add_parser("make_report")
    cmd.add_argument("truth_alignments")
    cmd.add_argument("compare_alignments", help="Comma-separeted list of files to compare")
    cmd.add_argument("-n", "--names", help="Comma-separated pretty readable names (corresponding to compare_alignments)")
    cmd.add_argument("colors", help="Comma-separated list of colors (must work with html)")
    cmd.add_argument("-f", "--report-id", required=False, default=None, help="Will be generated if not specified")
    cmd.set_defaults(func=make_html_report_wrapper)

    # Set correctness
    cmd = subparsers.add_parser("set_correctness")
    cmd.add_argument("truth_alignments")
    cmd.add_argument("alignments")
    cmd.set_defaults(func=set_correctness)

    # rename fq file to numeric increasing ids
    cmd = subparsers.add_parser("rename")
    cmd.add_argument("-q", "--fq", required=False)
    cmd.add_argument("-p", "--posfile", required=False)
    cmd.set_defaults(func=rename)


    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args(args)
    args.func(args)
