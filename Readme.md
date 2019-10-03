# Numpy alignments 
Fast way of representing mapped read positions and comparing two sets of mapped reads using numpy. 

## Install
Clone repo and enter directory, then run:
```bash
pip3 install .
```

## Usage
Numpy alignments always needs to know the number of reads when it is run. This makes it able to efficiently 
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


Create ROC-plots:
```bash
numpy_alignments compare -f figure.png truth bwa

```

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