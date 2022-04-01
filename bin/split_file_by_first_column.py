#! /usr/bin/env python

"""
Usage:
    python split_file_by_first_column.py [IN_FILE] -o [OUT_DIR] -s [OUT_FILE_SUFFIXES]

note:
    The IN_FILE should be sorted by the first column.

"""

import argparse
import os
import os.path
from itertools import groupby
from operator import itemgetter


def parse_data(input_file):
    for line in input_file:
        data = line.rstrip('\n').split('\t')
        yield data


def create_parser():
    parser = argparse.ArgumentParser(description="note: The IN_FILE should be sorted by the first column.")
    parser.add_argument('input_file', type=argparse.FileType('r'))
    parser.add_argument('-o', '--out_dir', default='.')
    parser.add_argument('-s', '--suffix')

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    all_data = parse_data(args.input_file)

    for key, gp in groupby(all_data, key=itemgetter(0)):
        with open(os.path.join(args.out_dir, f"{key}{args.suffix}"), 'w') as out:
            for data in gp:
                print(*data[1:], sep='\t', file=out)
