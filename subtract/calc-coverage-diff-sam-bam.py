#! /usr/bin/env python
"""
TODO:
* verify start and end / do we have an off-by-one error?
* verify pileup start and end too
* think about how to speed up...
* can probably do with a generic 'cov' output, too.
* check that both BAMs have the same reference name...
"""
import sys
import argparse
import pysam
import csv


def main():
    p = argparse.ArgumentParser()
    p.add_argument('all_reads')
    p.add_argument('query_reads')
    p.add_argument('-o', '--output', help='coverage CSV',
                   required=True)
    args = p.parse_args()

    all_bam = pysam.AlignmentFile(args.all_reads, "rb")
    query_bam = pysam.AlignmentFile(args.query_reads, "rb")

    outfp = open(args.output, 'w', newline='')
    w = csv.writer(outfp)
    w.writerow(['read_name','mapping_cov'])

    # build a location cache so we only need to query once for each location!
    location_cache = {}

    # iterate over query reads
    fup = query_bam.fetch()
    for n, read in enumerate(fup):
        if n % 100 == 0:
            print('...', n)

        chr = read.reference_name
        start = read.reference_start + 1
        end = read.reference_end

        # add /1 or /2 if read is paired and is read 1 or read 2
        if read.is_paired:
            if read.is_read1:
                if '/1' not in read.qname:
                    read.qname += '/1'
            elif read.is_read2:
                if '/2' not in read.qname:
                    read.qname += '/2'
            else:
                raise ValueError(f"read {read.qname} is paired but not read 1 or read 2")
        else:
            if '/1' not in read.qname:
                read.qname += '/1'

        # for each position in query read, get coverage
        sum_cov = []
        #print('XXX', chr, start, end)
        for i in range(start, end):
            # is it in location cache?
            depth = location_cache.get(i)

            # nope, calculate and save.
            if depth is None:
                pup = all_bam.pileup(chr, i, i+1, flag_filter=0)
                #pup = list(pup)

                # get set of unique reads (by name) that map to this position
                reads = set()
                for pupcol in pup:
                    for pup_read in pupcol.pileups:
                        reads.add(pup_read.alignment.query_name)
#                        print(pup_read.is_del, pup_read.is_refskip,
#                              pup_read.alignment.query_name,
#                              pup_read.query_position,
#                              pup_read.alignment.query_sequence[pup_read.query_position])

                # depth is number of distinct reads that cover this location
                depth = len(reads)
                location_cache[i] = depth

            sum_cov.append(depth)

        #print(sum_cov)
        w.writerow([read.qname, f"{sum(sum_cov) / len(sum_cov):.2f}"])

    outfp.close()


if __name__ == '__main__':
    main()
