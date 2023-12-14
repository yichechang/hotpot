# Calculate correlation between mapped reads and feature annotation


## Summary

This pipeline takes mapped and indexed alignment files and calculate
correlation between mapped reads and Alu annotation across the 
genome binned at given window sizes.


## Detail

***What Alu sites are used?***
This is determined by the following two:
1. A FASTA file containing names of Alu subfamilies to be processed.
2. *rmsk* BED file descriving genomic intervals for repeats.

The pipeline will pull all intervals (from BED file) that 
correspond to specified Alu element names (from FASTA file).


***What is the bin size?***
You can specify one or more in the configuration file: 
- `binsizes`


## Required input files
- Sorted and indexed BAM files
- An *rmsk* BED file
- A FASTA file containing desired Alu subfamily names
- A file describing chromosome names and lengths for a given genome