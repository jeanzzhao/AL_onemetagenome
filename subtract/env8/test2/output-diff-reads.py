#! /usr/bin/env python
import screed
import argparse
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument('starting_reads')
    p.add_argument('remove_reads')
    p.add_argument('-o', '--output-remaining-reads', required=True)
    args = p.parse_args()

    remove_reads = set()
    for record in screed.open(args.remove_reads):
        remove_reads.add(record.name)

    with open(args.output_remaining_reads, 'wt') as fp:
        for record in screed.open(args.starting_reads):
            if record.name not in remove_reads:
                fp.write(f'@{record.name}\n{record.sequence}\n+\n{record.quality}\n')


if __name__ == '__main__':
    sys.exit(main())
