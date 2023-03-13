import logging
import numpy as np
import sys
from tqdm import tqdm
import bionumpy as bnp
from graph_read_simulator.simulation import MultiChromosomeCoordinateMap
from shared_memory_wrapper import from_file, to_file
from bionumpy.datatypes import BamEntry
from bionumpy.bnpdataclass import bnpdataclass


def encode_chromosome(chromosome):
    if chromosome.startswith("chr"):
        return encode_chromosome(chromosome.replace("chr", ""))

    if chromosome == "X":
        chromosome = 23
    elif chromosome == "Y":
        chromosome = 24
    elif chromosome == "*":
        # Read couldn't map, ignore
        return -1
    else:
        chromosome = int(chromosome)

    return chromosome


def name_to_id(name):
    if "/" in name:
        name = name.split("/")
        read_id = int(name[0])
        pair_id = int(name[1])
        return read_id * 2 + pair_id - 1
    else:
        return int(name)


class NumpyAlignments:
    def __init__(self, chromosomes, positions, n_variants, scores, mapqs, is_correct=None):
        self.chromosomes = chromosomes
        self.positions = positions
        self.scores = scores
        self.mapqs = mapqs
        self.n_variants = n_variants
        self.is_correct = is_correct # indexes of correct alignments


    def __getitem__(self, item):
        data = {
            "chromosome": self.chromosomes[item],
            "position": self.positions[item],
            "score": self.scores[item],
            "mapq": self.mapqs[item],
            "n_variants": self.n_variants[item],
            "is_correct": self.is_correct[item] if self.is_correct is not None and len(self.is_correct) > 0 else None,
        }
        return data

    def set_correctness(self, truth_alignments, force=False, allowed_mismatch=150):
        if not force and self.is_correct is not None and len(self.is_correct) == len(self.positions):
            logging.info("Not setting correctness. Is set before")
            return

        logging.info("Allowing %d base pairs mismatch" % allowed_mismatch)
        # Sets which alignments are correctly align by checking against another alignment set
        self.is_correct = np.zeros(len(self.chromosomes), dtype=np.uint8)
        self.n_variants = truth_alignments.n_variants
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
        show_error = True
        chromosomes = np.zeros(n_alignments, dtype=np.uint8)
        positions = np.zeros(n_alignments, dtype=np.int32)  # int and not uint so we can subtract positions later
        n_variants = np.zeros(n_alignments, dtype=np.uint8)
        scores = np.zeros(n_alignments, dtype=np.uint16)
        mapqs = np.zeros(n_alignments, dtype=np.uint8)

        is_paired_end = False

        i = 0
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

            if l[6] != "*":
                if not is_paired_end:
                    logging.info("Assuming sam is paired end. Will assign IDs automatically based on line number")
                is_paired_end = True

            if is_paired_end:
                # hacky, should be fixed
                if "/" in l[0]:
                    identifier = name_to_id(l[0])
                else:
                    # this is the id after mapping without the /
                    identifier = name_to_id(l[0])*2
                    flag = int(l[1])
                    if flag >= 128:  # second in pair
                        identifier += 1
            else:
                identifier = name_to_id(l[0])

            chromosome = encode_chromosome(l[2])


            try:
                score = int(l[13].replace("AS:i:", ""))
            except ValueError:
                score = 0
                #logging.error("Could not parsed score from line. Skipping")
                #logging.error(line)
                #continue
            except IndexError:
                score = 0
                if show_error:
                    logging.error("Could not get score from line. Setting score to 0")
                    logging.error(line)
                show_error = False

            position = int(l[3])

            try:
                chromosomes[identifier] = chromosome
                positions[identifier] = position
                scores[identifier] = score
                mapqs[identifier] = int(l[4])
            except IndexError:
                logging.error("Got indexerror when parsing line. Skipping")
                logging.error(line)

            if "NVARIANTS:" in line:
                has_variant = 1
                if "NVARIANTS:i:0" in line:
                    has_variant = 0
                n_variants[identifier] = has_variant

            i += 1

        logging.info("Done getting alignments")
        return cls(chromosomes, positions, n_variants, scores, mapqs)


    @classmethod
    def from_bed(cls, n_alignments):
        chromosomes = np.zeros(n_alignments, dtype=np.uint8)
        positions = np.zeros(n_alignments, dtype=np.int32)
        n_variants = np.zeros(n_alignments, dtype=np.uint8)
        scores = np.zeros(n_alignments, dtype=np.uint16)
        mapqs = np.zeros(n_alignments, dtype=np.uint8)

        for line in tqdm(sys.stdin, total=n_alignments):

            l = line.split()
            try:
                identifier = name_to_id(l[3])
            except IndexError:
                logging.error("Cannot parse line %s" % line)
                raise
                sys.exit()

            chromosome = encode_chromosome(l[0])

            position = int(l[1])
            chromosomes[identifier] = chromosome
            positions[identifier] = position

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
                identifier = name_to_id(l[0])
            except IndexError:
                logging.error("Cannot parse line %s" % line)
                raise
                sys.exit()

            chromosome = encode_chromosome(l[1])

            position = int(l[2])
            if len(l) > 7:
                variants = int(l[7])
            else:
                variants = 0

            chromosomes[identifier] = chromosome
            positions[identifier] = position
            n_variants[identifier] = variants

        return cls(chromosomes, positions, n_variants, scores, mapqs)

    @classmethod
    def from_vgpos(cls, n_alignments):
        chromosomes = np.zeros(n_alignments, dtype=np.uint8)
        positions = np.zeros(n_alignments, dtype=np.int32)
        n_variants = np.zeros(n_alignments, dtype=np.uint8)
        scores = np.zeros(n_alignments, dtype=np.uint16)
        mapqs = np.zeros(n_alignments, dtype=np.uint8)

        for i, line in enumerate(tqdm(sys.stdin, total=n_alignments)):

            l = line.split()
            identifier = i
            chromosome = l[2]
            if chromosome == "X":
                chromosome = 23
            elif chromosome == "Y":
                chromosome = 24
            elif chromosome == "null":
                chromosome = 0
            else:
                try:
                    chromosome = int(chromosome)
                except ValueError:
                    logging.error("Could not parse chromosome %s. Setting to 0" % (chromosome))
                    chromosome = 0

            if l[3] == "null":
                position = 0
            else:
                position = int(l[3])

            chromosomes[identifier] = chromosome
            positions[identifier] = position

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
            identifier = name_to_id(l[0])
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
                is_correct=self.is_correct if self.is_correct is not None else np.array([]))
        logging.info("Saved to %s" % file_name)

    @classmethod
    def from_file(cls, file_name):
        try:
            data = np.load(file_name)
        except FileNotFoundError:
            data = np.load(file_name + ".npz")

        is_correct = None
        if "is_correct" in data:
            is_correct = data["is_correct"]

        return cls(data["chromosomes"], data["positions"].astype(np.int32), data["n_variants"], data["scores"], data["mapqs"], is_correct)

    def compare(self, other):
        pass


# Stores alignments sorted by "base name" and primary/seconday
# enables comparisons of alignments from different sources since these will be
# sorted in the same way (e.g. some sources will omitt paired-end information)

@bnpdataclass
class CustomBamEntry(BamEntry):
    base_name: bnp.encodings.BaseEncoding
    pair_id: int


# Tmp class that will replace NumpyAlignments
class NumpyAlignments2(NumpyAlignments):
    def __init__(self, data, n_variants=None, is_preprocessed=False):
        self.data = data
        self.n_variants = n_variants
        self.is_correct = None

        self._is_preprocessed = is_preprocessed
        if self._is_preprocessed:
            return
        self.preprocess()

    @classmethod
    def from_bam(cls, bam_file_name):
        data = bnp.open(bam_file_name).read()
        logging.info("%d alignments in bam" % len(data))
        data = data[data.flag < 256]  # remove seconday alignments
        logging.info("%d alignments after removing secondary alignments" % len(data))
        return cls(data)

    @classmethod
    def from_bam_and_nvariants_txt(cls, bam_file_name, nvariants_file_name):
        data = bnp.open(bam_file_name).read()
        logging.info("%d alignments in bam" % len(data))
        data = data[data.flag < 256]  # remove seconday alignments
        logging.info("%d alignments after removing secondary alignments" % len(data))
        n_variants = np.array([int(line.strip()) for line in open(nvariants_file_name)])
        return cls(data, n_variants)

    def to_file(self, file_name):
        return to_file([self.data.shallow_tuple(), self.n_variants], file_name)

    @classmethod
    def from_file(cls, file_name):
        data, n_variants = from_file(file_name)
        return cls(CustomBamEntry(*data), n_variants, is_preprocessed=True)

    def preprocess(self):

        logging.info("Preprocessing %d alignments" % len(self.data))
        # add base_name field with base_name (not including paired-end information if available)
        base_names = [str(name).split("/")[0] for name in self.data.name]

        # add pair-id (0 or 1)
        binary_flags = [bin(flag)[2:] for flag in self.data.flag]
        pair_ids = [0 if len(flag) < 8 or flag[-8] == 0 else 1 for flag in binary_flags]
        #new_data = self.data.add_fields({"base_name": base_names, "pair_id": pair_ids},
        #                                {"base_name": bnp.encodings.BaseEncoding, "pair_id": int})
        fields = self.data.shallow_tuple()
        new_data = CustomBamEntry(*self.data.shallow_tuple(), base_names, pair_ids)

        # sort alignments on base name and pair-id
        # after this sorting, these alignments can be compared to any bam with the same alignments
        sorting = sorted(range(len(new_data)), key=lambda i: (str(new_data[i].base_name), new_data[i].pair_id))
        sorted_data = new_data[sorting]
        self.data = sorted_data
        if self.n_variants is not None:
            self.n_variants = self.n_variants[sorting]

        self._is_preprocessed = True

    # methods to implement NumpyAlignments infterface
    @property
    def chromosomes(self):
        return self.data.chromosome

    @property
    def positions(self):
        return self.data.position

    @property
    def mapqs(self):
        return self.data.mapq

    @property
    def scores(self):
        return NotImplemented

    def set_correctness(self, truth_alignments, force=False, allowed_mismatch=150):
        assert len(truth_alignments.data) == len(self.data), "Truth alignments does not have same number of alignments (%d != %d)" % (len(truth_alignments.data), len(self.data))
        if not force and self.is_correct is not None and len(self.is_correct) == len(self.positions):
            logging.info("Not setting correctness. Is set before")
            return

        logging.info("Allowing %d base pairs mismatch" % allowed_mismatch)
        self.is_correct = np.zeros(len(self.chromosomes), dtype=np.uint8)
        self.n_variants = truth_alignments.n_variants

        chromosome_match = bnp.str_equal(self.chromosomes, truth_alignments.chromosomes)
        position_match = np.abs(self.positions - truth_alignments.positions) <= allowed_mismatch
        match = np.where(chromosome_match & position_match)[0]

        logging.info("Number of matches: %d" % len(match))
        self.is_correct[match] = 1
        logging.info("N correct: %d" % len(match))
