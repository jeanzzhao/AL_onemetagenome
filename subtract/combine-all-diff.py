#! /usr/bin/env python
import sys
import argparse
import csv
import os


def main():
    p = argparse.ArgumentParser()
    p.add_argument('all_csv')
    p.add_argument('diff_csv')
    p.add_argument('-o', '--output', required=True)
    args = p.parse_args()

    assert os.path.exists(args.all_csv)
    assert os.path.exists(args.diff_csv)

    all_d = {}
    with open(args.all_csv, 'r', newline="") as all_fp:
        r = csv.DictReader(all_fp)

        for row in r:
            name = row['read_name']
            all_d[name] = row
    print(f"read {len(all_d)} rows from 'all' file '{args.all_csv}'")

    diff_d = {}
    with open(args.diff_csv, 'r', newline="") as diff_fp:
        r = csv.DictReader(diff_fp)

        for row in r:
            name = row['read_name']
            diff_d[name] = row
    print(f"read {len(diff_d)} rows from 'diff' file '{args.diff_csv}'")

    all_s = set(all_d)
    diff_s = set(diff_d)
    assert not (diff_s - all_s)

    # check - how many paired?
    total_1 = 0
    total_pair = 0
    for name in all_d:
        if name.endswith('/1'):
            total_1 += 1
            pair = name.rstrip('1') + '2'
            if pair in all_d:
                total_pair += 1

    print(f"of {total_1} reads in 'all' ending in /1, {total_pair} have paired reads")

    # check - how many unpaired /2?
    total_2 = 0
    total_pair = 0
    for name in all_d:
        if name.endswith('/2'):
            total_2 += 1
            pair = name.rstrip('2') + '1'
            if pair in all_d:
                total_pair += 1

    print(f"of {total_2} reads in 'all' ending in /2, {total_pair} have paired reads")

    # build a set of distinct read names w/o /1 or /2
    uniq_names = set()
    for name in all_d:
        assert name.endswith('/1') or name.endswith('/2')
        name = name[:-2]
        uniq_names.add(name)

    fp = open(args.output, 'w', newline='')
    w = csv.writer(fp)

    header_row = [
        'all_name1',
        'all_1_cov',
        'all_1_cigar',
        'all_1_mapq',
        'all_1_readlen',
        'all_1_align_f',
        'all_1_is_pair',
        'all_1_is_primary',
        'all_name2',
        'all_2_cov',
        'all_2_cigar',
        'all_2_mapq',
        'all_2_readlen',
        'all_2_align_f',
        'all_2_is_pair',
        'all_2_is_primary',
        'diff_1_present',
        'diff_1_cov',
        'diff_1_cigar',
        'diff_1_mapq',
        'diff_1_readlen',
        'diff_1_align_f',
        'diff_1_is_pair',
        'diff_1_is_primary',
        'diff_2_present',
        'diff_2_cov',
        'diff_2_cigar',
        'diff_2_mapq',
        'diff_1_readlen',
        'diff_1_align_f',
        'diff_2_is_pair',
        'diff_2_is_primary',
        ]
    w.writerow(header_row)

    transfer_cols = ('mapping_cov', 'cigar', 'mapping_quality',
                     'read_length', 'read_align_f',
                     'is_proper_pair', 'is_primary_alignment')

    # make a gigantic munge of information!!
    for name in uniq_names:
        name_1 = name + '/1'
        name_2 = name + '/2'

        all_paired = 0
        if name_1 in all_d and name_2 in all_d:
            all_paired = 1

        diff_paired = 0
        if name_1 in diff_d and name_2 in diff_d:
            diff_paired = 1
   
        output_row = []
        output_row.append(name_1)
        for k in transfer_cols:
            if name_1 in all_d:
                output_row.append(all_d[name_1][k])
            else:
                output_row.append('')

        output_row.append(name_2)
        for k in transfer_cols:
            if name_2 in all_d:
                output_row.append(all_d[name_2][k])
            else:
                output_row.append('')

        if name_1 in diff_d:
            output_row.append('diff_1_PRESENT')
        else:
            output_row.append('diff_1_ABSENT')

        for k in transfer_cols:
            if name_1 in diff_d:
                output_row.append(diff_d[name_1][k])
            else:
                output_row.append('')

        if name_2 in diff_d:
            output_row.append('diff_2_PRESENT')
        else:
            output_row.append('diff_2_ABSENT')

        for k in transfer_cols:
            if name_2 in diff_d:
                output_row.append(diff_d[name_2][k])
            else:
                output_row.append('')

        w.writerow(output_row)

        assert len(output_row) == len(header_row), (len(output_row), len(header_row))


if __name__ == '__main__':
    sys.exit(main())
