#! /usr/bin/env python
"""
TODO:
* worry about partial alignments?

note pysam docs:
https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.qname
"""
import sys
import argparse
import pysam
import csv


def main():
    p = argparse.ArgumentParser()
    p.add_argument('depth_file')
    p.add_argument('query_reads')
    p.add_argument('-o', '--output', help='coverage CSV',
                   required=True)
    args = p.parse_args()

    coverages = {}
    with open(args.depth_file, 'r', newline='') as fp:
        r = csv.reader(fp, delimiter='\t')
        for row in r:
            # unpack row
            contig, pos, depth = row

            # convert to int
            pos, depth = int(pos), int(depth)

            # find depth per contig
            d2 = coverages.get(contig, {})

            # check that it's not being overwritten - basic consistency check
            assert pos not in d2

            # write!
            d2[pos] = depth
            coverages[contig] = d2

    query_bam = pysam.AlignmentFile(args.query_reads, "rb")

    outfp = open(args.output, 'w', newline='')
    w = csv.writer(outfp)
    w.writerow(['read_name','mapping_cov', 'cigar', 'mapping_quality',
                'read_length', 'read_align_f',
                'is_proper_pair', 'is_primary_alignment'])

    # iterate over query reads
    fup = query_bam.fetch()
    for n, read in enumerate(fup):
        if n % 100 == 0:
            print('...', n)

        chr = read.reference_name
        start = read.reference_start + 1 # @CTB
        end = read.reference_end
        assert start < end

        #print('XXX', chr, start, end)

        # add /1 or /2 if read is paired and is read 1 or read 2
        if read.is_paired and not ('/1' in read.qname or '/2' in read.qname):
            if read.is_read1:
                read.qname += '/1'
            elif read.is_read2:
                read.qname += '/2'
            else:
                raise ValueError(f"read {read.qname} is paired but not read 1 or read 2")
        elif not ('/1' in read.qname or '/2' in read.qname):
            read.qname += '/1'

        # for each position in query read, get coverage
        sum_cov = []
        for pos in range(start, end + 1):
            d2 = coverages.get(chr)
            if d2:
                depth = d2.get(pos, 0)
                sum_cov.append(depth)

        if not sum_cov:
            sum_cov = [0]

        # calculate fraction aligned/mapped, fraction soft clipped
        n_match = 0
        n_softclip = 0
        for k, v in read.cigartuples:
            if k == 0:          # M
                n_match += v
            elif k == 4:
                n_softclip += v
            else:
                #print(f"unhandled cigarstring operation: {(k, v)}")
                pass

#        print(len(sum_cov), read.query_length, sum_cov)
        w.writerow([read.qname, f"{sum(sum_cov) / len(sum_cov):.2f}",
                    read.cigarstring, read.mapping_quality,
                    read.query_length, n_match / read.query_length,
                    '1' if read.is_proper_pair else '0',
                    '0' if read.is_secondary else '1'])

    outfp.close()


if __name__ == '__main__':
    main()
