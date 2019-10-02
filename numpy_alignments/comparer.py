import logging
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class Comparer:
    def __init__(self, truth_alignments, compare_alignments, type='all'):
        self.truth_alignments = truth_alignments
        self.compare_alignments = compare_alignments
        self.type = type


    def create_roc_plots(self, save_to_file=None):
        mapq_intervals = [60, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 46, 44, 42, 40, 37, 34, 31, 27, 25, 23, 20, 17, 14, 11, 7, 5, 3, 2, 1, 0]

        for name, alignments in self.compare_alignments.items():
            logging.info("Setting corectness for %s" % name)
            alignments.set_correctness(self.truth_alignments)


        recalls = {name: [] for name in self.compare_alignments.keys()}
        precision = {name: [] for name in self.compare_alignments.keys()}

        for name, alignments in self.compare_alignments.items():
            logging.info("Processing %s" % name)
            if self.type == "all":
                total = len(self.truth_alignments.n_variants)
            elif self.type == "variants":
                total = len(np.where(self.truth_alignments.n_variants > 0)[0])
            elif self.type == "nonvariants":
                total = len(np.where(self.truth_alignments.n_variants == 0)[0])
            else:
                raise Exception("Invalig type %s" % type)

            #total = len(np.where(self.truth_alignments.n_variants > -1)[0])
            recalled_total = 0
            n_wrong_total = 0

            prev_interval = 100
            for mapq in mapq_intervals:
                logging.info("mapq: %d" % mapq)
                upper_limit = prev_interval
                lower_limit = mapq

                alignments = self.compare_alignments[name]
                if self.type == "all":
                    selection = np.where((alignments.mapqs >= lower_limit) & (alignments.mapqs < upper_limit))[0]
                elif self.type == "variants":
                    selection = np.where((alignments.mapqs >= lower_limit) & (alignments.mapqs < upper_limit) & (self.truth_alignments.n_variants > 0))[0]
                elif self.type == "nonvariants":
                    selection = np.where((alignments.mapqs >= lower_limit) & (alignments.mapqs < upper_limit) & (self.truth_alignments.n_variants == 0))[0]

                print("N with mapq >= %d and < %d: %d" % (lower_limit, upper_limit, len(selection)))

                # Find number of recalled and number of wrong here
                n_correct = np.sum(alignments.is_correct[selection])
                n_wrong = len(selection) - n_correct

                recalled_total += n_correct
                n_wrong_total += n_wrong

                recalls[name].append(recalled_total / total)
                precision[name].append(n_wrong_total / total)

                prev_interval = mapq

        plt.figure()
        for name in self.compare_alignments.keys():
            plt.plot(precision[name], recalls[name], label=name)
            plt.xscale("log")
            plt.legend()

        if save_to_file is not None:
            plt.savefig(save_to_file)
            logging.info("Saved plot to %s" % save_to_file)
        else:
            plt.show()

        print(recalls)
        print(precision)

    def get_wrong_alignments_correct_by_other(self, wrong_by, correct_by):
        selection_correct = set(list(np.nonzero(self.compare_alignments[correct_by].is_correct)[0]))
        selection_wrong = set(list(np.where(self.compare_alignments[wrong_by].is_correct == 0)[0]))
        selection_mapq = set(list(np.where(self.compare_alignments[wrong_by].mapqs >= 30)[0]))

        selection = selection_correct.intersection(selection_wrong).intersection(selection_mapq)

        print(selection)
        print(len(selection))
