import pandas as pd
import pysam


def main(infile, outfile):
    sorted_bamfile_path = infile

    mapped_lengths = []

    # Iterate through each read in the BAM file
    with pysam.AlignmentFile(sorted_bamfile_path) as bamfile:
        for read in bamfile:
            if read.is_proper_pair and not read.is_secondary and not read.is_supplementary:
                # Calculate the mapped length
                if read.is_read1:
                    # Assuming read2 comes after read1 in the file,
                    # which should be the case after sorting by query names
                    read2 = next(bamfile)
                    # Do not assume read1 always is upstream (5'-) to read2
                    start = min(read.reference_start, read2.reference_start)
                    end = max(read.reference_end, read2.reference_end)
                    mapped_length = end - start
                    mapped_lengths.append(mapped_length)


    # Export to csv
    mapped_lengths = pd.Series(mapped_lengths, name='fragment_length')
    mapped_lengths.to_csv(outfile, index=False)

if __name__ == '__main__':
    main(snakemake.input[0], snakemake.output[0])