#! /home/zyzhao/miniconda3/bin/python
"""#! /usr/bin/env python
"""
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
        if n % 10 == 0:
            print('...', n)

        chr = read.reference_name
        start = read.reference_start
        end = read.reference_end

        # for each position in query read, get coverage
        sum_cov = []
        for i in range(start, end + 1):
            # is it in location cache?
            cov = location_cache.get(i)

            # nope, calculate and save.
            if cov is None:
                pup = all_bam.pileup(chr, i, i+1)
                cov = len(list(pup))

                location_cache[i] = cov

            sum_cov.append(cov)

        w.writerow([read.qname, f"{sum(sum_cov) / len(sum_cov):.2f}"])


if __name__ == '__main__':
    main()
