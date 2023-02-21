import logging
import numpy as np
import plotly.graph_objects as go
import plotly



class Comparer:
    def __init__(self, truth_alignments, compare_alignments, colors=None, type='all', allowed_mismatch=150):
        self.truth_alignments = truth_alignments
        self.compare_alignments = compare_alignments
        self.type = type
        self.mapq_intervals = [60, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 46, 44, 42, 40, 37, 34, 30, 27, 25, 23, 20, 17, 14, 10, 6, 3, 2, 1, 0]
        logging.info("Allowed mismatch in bp is %d" % allowed_mismatch)
        self.allowed_mismatch = allowed_mismatch

        default_colors = ["blue", "red", "green", "purple", "orange", "black"]
        if colors is None:
            colors = {c: default_colors[i] for i, c in enumerate(self.compare_alignments.keys())}

        self.colors = colors

    def get_correct_rates(self, min_mapq=0):

        for name, alignments in self.compare_alignments.items():
            logging.info("Setting corectness for %s, allowed mismatch: %d" % (name, self.allowed_mismatch))
            alignments.set_correctness(self.truth_alignments, allowed_mismatch=self.allowed_mismatch)

        if self.type == "all":
            n_alignments = len(self.truth_alignments.positions)
        elif self.type == "variants":
            n_alignments = len(np.where((self.truth_alignments.n_variants > 0))[0])
        elif self.type == "nonvariants":
            n_alignments = len(np.where(self.truth_alignments.n_variants == 0)[0])
        else:
            raise Exception("Invalid type (must be all, variants or nonvariants)")

        rates = {}
        for name, alignments in self.compare_alignments.items():
            logging.info("Processing %s" % name)
            compare = self.compare_alignments[name]
            if self.type == "all":
                selection = compare.positions[np.where((compare.is_correct == 1) & (compare.mapqs >= min_mapq))[0]]
                selection_wrong = compare.positions[np.where((compare.is_correct == 0) & (compare.mapqs >= min_mapq))[0]]
            elif self.type == "variants":
                selection = compare.positions[np.where((compare.is_correct == 1) & (compare.n_variants > 0) & (compare.mapqs >= min_mapq))[0]]
                selection_wrong = compare.positions[np.where((compare.is_correct == 0) & (compare.n_variants > 0) & (compare.mapqs >= min_mapq))[0]]
            elif self.type == "nonvariants":
                selection = compare.positions[np.where((compare.is_correct == 1) & (compare.n_variants == 0) & (compare.mapqs >= min_mapq))[0]]
                selection_wrong = compare.positions[np.where((compare.is_correct == 0) & (compare.n_variants == 0) & (compare.mapqs >= min_mapq))[0]]

            n_correct = len(selection)
            n_wrong = len(selection_wrong)

            try:
                rates[name] = (n_correct / n_alignments, (n_wrong / (n_wrong + n_correct)))  # np.sum(self.compare_alignments[name].is_correct) / len(self.truth_alignments.positions)
            except ZeroDivisionError:
                logging.error("Name: %s, type: %s" % (name, self.type))
                logging.error("Got zerodivision error. N correct: %d, n_alignments: %d. N wrong: %d" % (n_correct, n_alignments, n_wrong))
                rates[name] = (0, 0)
                #raise

        return rates

    def create_roc_plots(self, save_to_file=None, limit_comparison=None):
        mapq_intervals = self.mapq_intervals

        for name, alignments in self.compare_alignments.items():
            logging.info("Setting corectness for %s" % name)
            alignments.set_correctness(self.truth_alignments)


        recalls = {name: [] for name in self.compare_alignments.keys()}
        precision = {name: [] for name in self.compare_alignments.keys()}
        n_reads = {name: [] for name in self.compare_alignments.keys()}

        for name, alignments in self.compare_alignments.items():
            logging.info("Processing %s" % name)
            if self.type == "all":
                total = len(alignments.n_variants)
            elif self.type == "variants":
                total = len(np.where(alignments.n_variants > 0)[0])
            elif self.type == "nonvariants":
                total = len(np.where(alignments.n_variants == 0)[0])
            else:
                raise Exception("Invalig type %s" % type)

            #total = len(np.where(self.truth_alignments.n_variants > -1)[0])
            recalled_total = 0
            n_wrong_total = 0

            prev_interval = 100
            for mapq in mapq_intervals:
                #logging.info("mapq: %d" % mapq)
                upper_limit = prev_interval
                lower_limit = mapq

                alignments = self.compare_alignments[name]
                if self.type == "all":
                    selection = np.where((alignments.mapqs >= lower_limit) & (alignments.mapqs < upper_limit))[0]
                elif self.type == "variants":
                    selection = np.where((alignments.mapqs >= lower_limit) & (alignments.mapqs < upper_limit) & (alignments.n_variants > 0))[0]
                elif self.type == "nonvariants":
                    selection = np.where((alignments.mapqs >= lower_limit) & (alignments.mapqs < upper_limit) & (alignments.n_variants == 0))[0]

                if limit_comparison is not None:
                    selection = selection[0:limit_comparison]
                    logging.warning("Limiting comparison to max %d reads" % limit_comparison)

                #print("N with mapq >= %d and < %d: %d" % (lower_limit, upper_limit, len(selection)))

                # Find number of recalled and number of wrong here
                n_correct = np.sum(alignments.is_correct[selection])
                n_wrong = len(selection) - n_correct

                recalled_total += n_correct
                n_wrong_total += n_wrong

                recalls[name].append(recalled_total / total)
                precision[name].append((n_wrong_total + 1) / (n_wrong_total + recalled_total))
                n_reads[name].append(np.log(1 + len(selection)))


                prev_interval = mapq

        fig = go.Figure(
            layout=go.Layout(
                xaxis=dict(showgrid=True, zeroline=True),
            )
        )
        #fig.update_layout(title_text="Reads (%s)" % self.type)
        fig.update_layout(xaxis_type="log")
        fig.update_layout(
            xaxis = dict(
                showexponent = 'all',
                exponentformat = 'e',
                nticks=5,
                tickfont=dict(
                    size=16
                )
        ))
        fig.update_layout(
            yaxis = dict(
                tickfont=dict(
                    size=16
                )
            ),
            xaxis_title="#wrong mapped / #mapped",
            yaxis_title="#correctly mapped / total"
        )
        fig.update_layout(showlegend=False)



        ticker_positions = ["top left", "bottom left", "top left", "bottom left", "top left", "bottom left", "top left", "bottom left"]
        for i, name in enumerate(self.compare_alignments.keys()):
            ticker_text = [m if m % 10 == 0 else None for m in self.mapq_intervals]


            sizes = np.array(n_reads[name])
            fig.add_trace(go.Scatter(x=precision[name], y=recalls[name],
                                     text=ticker_text,
                                     mode='markers+lines+text',
                                     textposition=ticker_positions[i],
                                     name=name,
                                     textfont=dict(
                                        size=16
                                     ),
                                     marker=dict(
                                         sizemode='area',
                                         sizeref=4. * max(sizes) / (40. ** 2),
                                         sizemin=7,
                                         color=self.colors[name]
                                         )
                                     ),
                          )

        if save_to_file is not None:
            #plt.savefig(save_to_file)
            plotly.offline.plot(fig, filename=save_to_file, auto_open=False)
            logging.info("Saved plot to %s" % save_to_file)
        else:
            fig.show()


    def get_wrong_alignments_correct_by_other(self, wrong_by, correct_by):
        selection_correct = set(list(np.nonzero(self.compare_alignments[correct_by].is_correct)[0]))
        selection_wrong = set(list(np.where(self.compare_alignments[wrong_by].is_correct == 0)[0]))
        selection_mapq = set(list(np.where(self.compare_alignments[wrong_by].mapqs >= 30)[0]))

        selection = selection_correct.intersection(selection_wrong).intersection(selection_mapq)

