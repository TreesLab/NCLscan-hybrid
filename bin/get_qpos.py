#! /usr/bin/env python

import argparse
import re
from collections import defaultdict
from itertools import accumulate
from bisect import bisect_left

import logging

logger = logging.getLogger(__name__)


_CIGAR_PAT = re.compile(r'(([0-9]+)([MIDNSHP=X]))')


def get_query_pos(cigar_str, ref_pos):
    cigar_list = re.findall(_CIGAR_PAT, cigar_str)
    logger.debug(cigar_list)

    cigar_list_ref = list(filter(lambda c: c[2] in ['M', 'D'], cigar_list))
    logger.debug(cigar_list_ref)

    cigar_list_ref_accum = list(
        accumulate(
            [
                int(c[1])
                for c in cigar_list_ref
            ]
        )
    )
    logger.debug(cigar_list_ref_accum)

    idx_ref_pos = bisect_left(cigar_list_ref_accum, ref_pos)
    ref_pos_cigar = cigar_list_ref[idx_ref_pos]
    logger.debug(ref_pos_cigar)

    if ref_pos_cigar[2] == 'D':
        ref_pos = cigar_list_ref_accum[idx_ref_pos-1]

    idx_in_cigar_list = cigar_list.index(ref_pos_cigar)
    logger.debug(idx_in_cigar_list)

    cigar_list_before_ref_pos = cigar_list[:idx_in_cigar_list]
    logger.debug(cigar_list_before_ref_pos)

    total_D = sum(
        map(
            lambda c: int(c[1]),
            filter(
                lambda c: c[2] == 'D',
                cigar_list_before_ref_pos
            )
        )
    )
    logger.debug(total_D)

    total_I = sum(
        map(
            lambda c: int(c[1]),
            filter(
                lambda c: c[2] == 'I',
                cigar_list_before_ref_pos
            )
        )
    )
    logger.debug(total_I)

    query_pos = ref_pos + total_I - total_D

    return query_pos


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', type=argparse.FileType('r'))
    parser.add_argument('CIGAR', nargs='?', default=None)
    parser.add_argument('ref_pos', nargs='?', type=int, default=None)
    parser.add_argument('--debug', action='store_true')

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        logger.addHandler(ch)

    if args.file:
        for line in args.file:
            cigar, ref_pos = line.rstrip('\n').split('\t')
            ref_pos = int(ref_pos)

            query_pos = get_query_pos(cigar, ref_pos)
            print(query_pos)

    elif args.CIGAR and args.ref_pos:
        query_pos = get_query_pos(args.CIGAR, args.ref_pos)
        print(query_pos)

    else:
        parser.print_help()
