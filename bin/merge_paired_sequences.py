#! /usr/bin/env python

from __future__ import print_function

import sys
import itertools as it
import re


def merge_paired_sequences(fasta_data):
    sorted_fasta_data = sorted(fasta_data, key=lambda fa: [fa.ID[:-2], int(fa.ID[-1])])
    
    result_fa = []
    for k1, g1 in it.groupby(sorted_fasta_data, lambda fa: fa.ID[:-2]):
        g1_list = list(g1)
        
        if len(g1_list) != 2:
            continue
        
        fa_1, fa_2 = g1_list

        if (fa_1.ID[-1] != '1') or (fa_2.ID[-1] != '2'):
            print(k1)
            break

        result_fa.append([k1, fa_1.seq + fa_2.seq])
    
    return result_fa
    
    
class Fasta(object):
    def __init__(self, raw_fa):
        raw_fa_list = raw_fa[1:].split('\n')
        
        self.ID = raw_fa_list[0]
        self.seq = ''.join(raw_fa_list[1:])
        self.raw_fa = raw_fa
        
if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        print("""\
Usage:
    merge_paired_sequences.py [in_file] [out_file]""", file=sys.stderr)
        
        exit(1)
    
    
    with open(sys.argv[1]) as data_reader:
        all_fasta_data = map(Fasta, re.findall(r"(>[^>]*)", data_reader.read()))
    
    result_fa = merge_paired_sequences(all_fasta_data)
    
    with open(sys.argv[2], 'w') as data_writer:
        for seq_id, seq in result_fa:
            print(">{}".format(seq_id), file=data_writer)
            print("{}".format(seq.upper()), file=data_writer)