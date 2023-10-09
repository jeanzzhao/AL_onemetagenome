import os
import csv
import pysam
import argparse

# Notes:
# - potentially make sure mapped section is at least x bp long / x% of read length (instead of just total % alignment)
# - potentially split failed reads into multiple files based on the REASON they failed (e.g. map quality, alignment %, etc.)


def percent_aligned(read):
    """
    Calculate the percentage of the read that is aligned based on its CIGAR string.

    Parameters:
    - read (pysam.AlignedSegment): A single read from a BAM/SAM file, represented as an AlignedSegment in pysam.

    Returns:
    - float: The percentage of the read that is aligned to the reference.

    The CIGAR string provides detailed information about how a read aligns to the reference. In pysam, the CIGAR string is represented as
    a list of tuples (cigartuples), where each tuple consists of an operation code and its length. The operation codes are integers that
    correspond to different alignment operations:

    - 0: M (alignment match or mismatch)
    - 7: X (sequence mismatch, a more specific version of M)
    - 8: = (sequence match)

    The function calculates the total number of bases in the read that align to the reference by summing the lengths of the operations
    that represent alignment matches (M, X, =). It then divides this total by the full length of the read to compute the alignment
    percentage.

    Note:
    This function assumes that the read has been aligned (i.e., is not unmapped) and that it has a valid CIGAR string.
    """
    aligned_length = sum([length for code, length in read.cigartuples if code in [0, 7, 8]])  # M, X, and = operations
    return (aligned_length / read.query_length) * 100

def write_csv(csv_writer, read):
    """
    Write the read data into the provided CSV writer.
    """
    csv_writer.writerow([read.query_name,
                         read.flag,
                         read.cigarstring,
                         read.mapping_quality,
                         read.reference_name,
                         read.reference_start,
                         read.reference_end])

def main(args):
    """
    Split reads from input_bam based on the given criteria.
    Reads that meet the criteria are written to pass_bam and pass_csv,
    while the rest are written to fail_bam and fail_csv.
    """
    infile_basename = (args.input).rsplit(".bam", 1)[0]
    # make parameter string for output files
    # pa70_mq30_pp
    p_str = f"pa{args.alignment_threshold}"
    if args.mapq:
        p_str += f'_mq{args.mapq}'
    if args.proper_pairs:
        p_str += '_pp'

    pass_bam_file = f"{infile_basename}.{p_str}.pass.bam"
    fail_bam_file = f"{infile_basename}.{p_str}.fail.bam"
    pass_csv_file = f"{infile_basename}.{p_str}.pass.csv"
    fail_csv_file = f"{infile_basename}.{p_str}.fail.csv"

    with pysam.AlignmentFile(args.input, "rb") as infile, \
         pysam.AlignmentFile(pass_bam_file, "wb", header=infile.header) as passbam, \
         pysam.AlignmentFile(fail_bam_file, "wb", header=infile.header) as failbam, \
         open(pass_csv_file, 'w', newline='') as pass_csv, \
         open(fail_csv_file, 'w', newline='') as fail_csv:

        pass_csv_writer = csv.writer(pass_csv)
        fail_csv_writer = csv.writer(fail_csv)

        # Write the header for the CSV files
        header = ["Name", "Flag", "CIGAR", "Mapping_Quality", "Reference_Seq", "Start_Position", "End_Position"]
        pass_csv_writer.writerow(header)
        fail_csv_writer.writerow(header)

        for read in infile:
            if not read.is_unmapped and percent_aligned(read) >= args.alignment_threshold and (args.mapq is None or read.mapping_quality >= args.mapq):
                if args.proper_pairs and not read.is_proper_pair:
                    failbam.write(read)
                    write_csv(fail_csv_writer, read)
                else:
                    passbam.write(read)
                    write_csv(pass_csv_writer, read)
            else:
                failbam.write(read)
                write_csv(fail_csv_writer, read)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Split reads based on alignment percentage, proper pair criteria, and mapping quality.")
    p.add_argument("-i", "--input", required=True, help="Path to the input BAM file.")
    p.add_argument("-a", "--alignment-threshold", type=float, default=70, help="Percentage threshold for alignment. Default is 70%%.")
    p.add_argument("-p", "--proper-pairs", action="store_true", help="Only pass read alignments that are properly paired.")
    p.add_argument("-q", "--mapq", type=int, help="Only pass read alignments that meet this minimum mapping quality.")

    args = p.parse_args()
    main(args)

