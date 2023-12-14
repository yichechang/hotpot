# Calculate mapped reads profile around each Alu subfamily site


## Summary

This pipeline takes mapped and indexed alignment files and calculate
profile and heatmap around each Alu subfamily annotated sites.


## Detail

***What Alu sites are used?***
This is determined by the following two:
1. A FASTA file containing names of Alu subfamilies to be processed.
2. *rmsk* BED file descriving genomic intervals for repeats.

The pipeline will pull all intervals (from BED file) that 
correspond to specified Alu element names (from FASTA file).


***What is the bin size?***
You can specify one or more in the configuration file: 
- `bam_to_bigwig` > `binsizes`


***What is the range considered for the profile?***
You can specify in the configuration file: 
- `compute_matrix` > `params` > `before` and `after`


## Required input files
- Sorted and indexed BAM files
- An *rmsk* BED file
- A FASTA file containing desired Alu subfamily names