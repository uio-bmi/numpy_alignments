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
bwa mem ref.fa reads.fa | numpy_alignments store sam some_name 100
```
