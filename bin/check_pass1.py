#! /usr/bin/env python


import argparse
from itertools import groupby
from operator import itemgetter


def parse_data(input_file):
    for line in input_file:
        data = line.rstrip('\n').split('\t')
        yield data


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=argparse.FileType('r'))

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    all_data = parse_data(args.input_file)

    for key, gp in groupby(all_data, key=itemgetter(0)):
        gp = list(gp)

        chr_d, _, _, chr_a, _, _ = key.split(':')

        count_d = len(list(filter(lambda data: data[1]==chr_d, gp)))
        count_a = len(list(filter(lambda data: data[1]==chr_a, gp)))

        if count_d and count_a:
            print(key)
