import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import sys
from tqdm import tqdm


class NumpyAlignments:
    def __init__(self, chromosomes, positions, n_variants, scores, mapqs, is_correct=None):
        self.chromosomes = chromosomes
        self.positions = positions
        self.scores = scores
        self.mapqs = mapqs
        self.n_variants = n_variants
        self.is_correct = None  # indexes of correct alignments


    def __getitem__(self, item):
        data = {
            "chromosome": self.chromosomes[item],
            "position": self.positions[item],
            "score": self.scores[item],
            "mapq": self.mapqs[item],
            "n_variants": self.n_variants[item],
            "is_correct": self.is_correct[item] if self.is_correct is not None else None,
        }
        return data

    def set_correctness(self, truth_alignments, allowed_mismatch=150):
        if self.is_correct is not None:
            logging.info("Not setting correctness. Is set before")
            return

        logging.info("Allowing %d base pairs mismatch" % allowed_mismatch)
        # Sets which alignments are correctly align by checking against another alignment set
        self.is_correct = np.zeros(len(self.chromosomes), dtype=np.uint8)
        #chromosome_match = set(np.where(self.chromosomes == truth_alignments.chromosomes)[0])
        #position_match = set(np.where(np.abs(self.positions - truth_alignments.positions) <= allowed_mismatch)[0])
        #match = np.array(list(chromosome_match.intersection(position_match)))

        match = np.where((self.chromosomes == truth_alignments.chromosomes) & (np.abs(self.positions - truth_alignments.positions) <= allowed_mismatch))[0]

        logging.info("Number of matches: %d" % len(match))
        self.is_correct[match] = 1
        logging.info("N correct: %d" % len(match))
        #self.is_correct[np.where((self.chromosomes == truth_alignments.chromosomes) &
        # (np.abs(self.positions - truth_alignments.positions) <= allowed_mismatch))[0]] = 1

    @classmethod
    def from_sam(cls, n_alignments):
        chromosomes = np.zeros(n_alignments, dtype=np.uint8)
        positions = np.zeros(n_alignments, dtype=np.int32)  # int and not uint so we can subtract positions later
        n_variants = np.zeros(n_alignments, dtype=np.uint8)
        scores = np.zeros(n_alignments, dtype=np.uint16)
        mapqs = np.zeros(n_alignments, dtype=np.uint8)

        for line in tqdm(sys.stdin, total=n_alignments):
            if line.startswith("@"):
                continue

            l = line.split()
            if len(l) < 2:
                logging.error("Cannot parse line")
                logging.error(line)
                continue

            if int(l[1]) >= 256:
                continue  # not primary mapping

            identifier = int(l[0])
            chromosome = l[2]
            if chromosome == "X":
                chromosome = 23
            elif chromosome == "Y":
                chromosome = 24
            elif chromosome == "*":
                # Read couldn't map, ignore
                continue
            else:
                chromosome = int(chromosome)

            try:
                score = int(l[13].replace("AS:i:", ""))
            except ValueError:
                score = 0
                #logging.error("Could not parsed score from line. Skipping")
                #logging.error(line)
                #continue
            except IndexError:
                logging.error("Could not parsed line. Skipping")
                logging.error(line)
                continue

            position = int(l[3])

            try:
                chromosomes[identifier] = chromosome
                positions[identifier] = position
                scores[identifier] = score
                mapqs[identifier] = int(l[4])
            except IndexError:
                logging.error("Got indexerror when parsing line. Skipping")
                logging.error(line)

        logging.info("Done getting alignments")
        return cls(chromosomes, positions, n_variants, scores, mapqs)

    @classmethod
    def from_truth(cls, n_alignments):
        chromosomes = np.zeros(n_alignments, dtype=np.uint8)
        positions = np.zeros(n_alignments, dtype=np.int32)
        n_variants = np.zeros(n_alignments, dtype=np.uint8)
        scores = np.zeros(n_alignments, dtype=np.uint16)
        mapqs = np.zeros(n_alignments, dtype=np.uint8)

        for line in tqdm(sys.stdin, total=n_alignments):

            l = line.split()
            try:
                identifier = int(l[0])
            except IndexError:
                logging.error("Cannot parse line %s" % line)
                raise
                sys.exit()

            chromosome = l[1]
            if chromosome == "X":
                chromosome = 23
            elif chromosome == "Y":
                chromosome = 24
            else:
                chromosome = int(chromosome)

            position = int(l[2])
            variants = int(l[7])

            chromosomes[identifier] = chromosome
            positions[identifier] = position
            n_variants[identifier] = variants

        return cls(chromosomes, positions, n_variants, scores, mapqs)

    @classmethod
    def from_pos(cls, n_alignments):
        chromosomes = np.zeros(n_alignments, dtype=np.uint8)
        positions = np.zeros(n_alignments, dtype=np.int32)
        n_variants = np.zeros(n_alignments, dtype=np.uint8)
        scores = np.zeros(n_alignments, dtype=np.uint16)
        mapqs = np.zeros(n_alignments, dtype=np.uint8)

        for line in tqdm(sys.stdin, total=n_alignments):

            l = line.split()
            identifier = int(l[0])
            chromosome = l[1]
            if chromosome == "X":
                chromosome = 23
            elif chromosome == "Y":
                chromosome = 24
            elif chromosome == "null":
                chromosome = 0
            else:
                chromosome = int(chromosome)

            if l[2] == "null":
                position = 0
            else:
                position = int(l[2])

            chromosomes[identifier] = chromosome
            positions[identifier] = position
            try:
                mapqs[identifier] = int(l[3])
                scores[identifier] = int(l[4])
            except IndexError:
                logging.warning("Could not get mapq or score on line %s" % line)

        return cls(chromosomes, positions, n_variants, scores, mapqs)

    def to_file(self, file_name):
        logging.info("Saving to file %s" % file_name)
        np.savez(file_name,
                chromosomes=self.chromosomes,
                positions=self.positions,
                scores=self.scores,
                mapqs=self.mapqs,
                n_variants=self.n_variants,
                is_correct=self.is_correct)
        logging.info("Saved to %s" % file_name)

    @classmethod
    def from_file(cls, file_name):
        data = np.load(file_name)
        is_correct = None
        if "is_correct" in data:
            is_correct = data["is_correct"]

        return cls(data["chromosomes"], data["positions"].astype(np.int32), data["n_variants"], data["scores"], data["mapqs"], is_correct)

    def compare(self, other):
        pass


if __name__ == "__main__":
    alignments = NumpyAlignments.from_pos_file(sys.argv[1], int(sys.argv[2])) 
    alignments.to_file(sys.argv[3])



