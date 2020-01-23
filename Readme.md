
# Numpy alignments.
Fast way of representing mapped read positions and comparing two sets of mapped reads using numpy.
This package can make nice report pretty fast, comparing mapping accuracies of different sets of alignments. 

## Install
Clone repo and enter directory, then run:
```bash
pip3 install .
```
.. or just pip install it (Python3):
```bash
pip3 install numpy_alignments
```

## Usage
Important: Numpy alignments always needs to know the number of reads when it is run. 
This package is specifically made to work with simulated reads, and requires all reads to have names from 0 to the number of reads
(such as reads simulated by [Graph Read Simulator](https://github.com/ivargr/graph_read_simulator/)). This makes it able to efficiently 
represent the reads in numppy arrays of the given size. The following examples assume there are 100 reads.

### Example
Map reads with BWA-MEM and pipe directly to numpy alignments to avoid storing large BAM-files on disk:

```python
bwa mem ref.fa reads.fa | numpy_alignments store sam bwa 100
```

Save truth positions:
```bash
cat positions.tsv | numpy_alignments store truth truth 265154
```

Compare bwa to truth:
```bash
numpy_alignments get_correct_rates truth bwa
```


Create html report:
```bash
numpy_alignments make_report -f my-report-name --names="bwa" truth bwa purple
```

You can specify more names, alignments sets and colors by separating them with commas.

### Use as python library
```python
from numpy_alignments.comparer import Comparer
from numpy_alignments import NumpyAlignments
bwa = NumpyAlignments.from_file("bwa.npz")
truth = NumpyAlignments.from_file("truth.npz")
comparer = Comparer(truth, {"bwa": bwa})
rate = comparer.get_correct_rates()
print(rate)

```
