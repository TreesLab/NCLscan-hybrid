#! /usr/bin/env python
from __future__ import print_function

import sys
import argparse
import re
import gzip
from collections import namedtuple
from itertools import groupby
from operator import itemgetter

__version__ = '0.7.1'

JuncSite = namedtuple("JuncSite", ['chr', 'pos', 'strand'])

class NCLevent:
    def __init__(self, raw_data):
        self.raw_data = raw_data[:]
        self.donor = JuncSite(raw_data[0], raw_data[1], raw_data[2])
        self.accepter = JuncSite(raw_data[3], raw_data[4], raw_data[5])


class ExonExtractor:
    def __init__(self):
        self.tid_pattern = re.compile('transcript_id "([^;]*)";')
        self.exon_number_pattern = re.compile('exon_number "?([0-9]*)"?;')
        self.transcript_type_pattern = re.compile('transcript_(?:bio)?type "([^;]*)";')
        self.transcript_version_pattern = re.compile('transcript_version "([^;]*)";')

    def _basic_info_getter(self, data):
        return [data[0], data[3], data[4], data[6]]

    def extract(self, data_iter):
        self.transcripts = {}
        self.exons = []

        for line in data_iter:

            if line.startswith('#'):
                continue

            data = line.rstrip('\n').split('\t')
            region_type = data[2]

            if region_type == 'exon':
                basic_info = self._basic_info_getter(data)
                exon_tid = re.search(self.tid_pattern, data[8]).group(1)
                exon_tid_version = re.search(self.transcript_version_pattern, data[8])
                if exon_tid_version:
                    exon_tid = "{}.{}".format(exon_tid, exon_tid_version.group(1))

                exon_number = re.search(self.exon_number_pattern, data[8]).group(1)

                transcript_type = re.search(self.transcript_type_pattern, data[8]).group(1)
                is_protein_coding = int(transcript_type == "protein_coding")

                exon_len = int(data[4]) - int(data[3]) + 1

                self.exons.append(basic_info + \
                                    [exon_tid, exon_number, exon_len, is_protein_coding])

        # generate the total length of exons in one transcript
        for tid, tid_gp in groupby(self.exons, key=itemgetter(4)):
            tid_gp = list(tid_gp)
            self.transcripts[tid] = (tid_gp[0][7], sum(map(itemgetter(6), tid_gp)))


class IntronExtractor:
    def _get_intron(self, exon1, exon2):
        assert exon1[4] == exon2[4]

        strand = exon1[3]
        if strand == '+':
            intron_start = int(exon1[2]) + 1
            intron_end = int(exon2[1]) - 1

        elif strand == '-':
            intron_start = int(exon2[2]) + 1
            intron_end = int(exon1[1]) - 1

        assert intron_start <= intron_end

        intron = [exon1[0], intron_start, intron_end, strand, exon1[4], exon1[5], exon2[5], intron_end - intron_start + 1]
        return intron

    def extract(self, exons_data):
        self.introns = []
        for _, tid_gp in groupby(exons_data, key=itemgetter(4)):
            tid_gp = list(tid_gp)
            for i in range(len(tid_gp) - 1):
                self.introns.append(self._get_intron(tid_gp[i], tid_gp[i+1]))


class JunctionSitesDB:
    def __init__(self):
        self.ncl_donor = {}
        self.ncl_accepter = {}

    def _get_donor_accepter(self, exons, di, ai):
        def group_exons(exs):
            for k, gp in groupby(sorted(exs, key=itemgetter(0)), key=itemgetter(0)):
                gp = [data[1:] for data in gp]
                yield [k, gp]

        donor = map(itemgetter(di, 4, 5), exons)
        accepter = map(itemgetter(ai, 4, 5), exons)

        grouped_donor = dict(group_exons(donor))
        grouped_accepter = dict(group_exons(accepter))

        return grouped_donor, grouped_accepter

    def generate_db(self, exon_data):
        for chrm, chrm_exons in groupby(exon_data, key=itemgetter(0)):
            chrm_exons = list(chrm_exons)
            plus_part = list(filter(lambda data: data[3]=='+', chrm_exons))
            minus_part = list(filter(lambda data: data[3]=='-', chrm_exons))

            plus_donor, plus_accepter = self._get_donor_accepter(plus_part, 2, 1)
            minus_donor, minus_accepter = self._get_donor_accepter(minus_part, 1, 2)

            self.ncl_donor[chrm] = {'+': plus_donor, '-': minus_donor}
            self.ncl_accepter[chrm] = {'+': plus_accepter, '-': minus_accepter}

    def get_junc_site(self, junc_site, donor_accepter):
        if donor_accepter == "donor":
            db = self.ncl_donor
        elif donor_accepter == "accepter":
            db = self.ncl_accepter
        else:
            raise ValueError("Only 'donor' or 'accepter' is valid!")

        try:
            res = db[junc_site.chr][junc_site.strand][junc_site.pos]
        except KeyError:
            print("No such {} site: {}".format(donor_accepter, junc_site), file=sys.stderr)
            return []
        else:
            return res


class FlankingIntronDB:
    def __init__(self):
        self.db = {}

    def generate_db(self, introns_data):
        for tid, tid_gp in groupby(introns_data, key=itemgetter(4)):
            tid_gp = list(tid_gp)
            donor_introns = dict(map(itemgetter(5, 7), tid_gp))
            accepter_introns = dict(map(itemgetter(6, 7), tid_gp))
            self.db[tid] = {'donor': donor_introns, 'accepter': accepter_introns}

    def get_flanking_intron(self, tid, exon_number, donor_accepter):
        try:
            intron_len = self.db[tid][donor_accepter][exon_number]
        except KeyError:
            return "no intron"
        else:
            return intron_len


class ExonsLengthDB:
    def __init__(self):
        self.db = {}

    def generate_db(self, exons_data):
        for tid, tid_gp in groupby(exons_data, key=itemgetter(4)):
            tid_gp = list(tid_gp)
            exon_lens = list(map(itemgetter(6), tid_gp))
            self.db[tid] = exon_lens

    def get_total_len_of_intermediate_exons(self, tid, exon_number_1, exon_number_2):
        start_ind = int(exon_number_1) - 1
        stop_ind = int(exon_number_2) - 1 + 1
        return sum(self.db[tid][start_ind:stop_ind])


class TranscriptExonsDB:
    def __init__(self):
        self.db = {}

    def generate_db(self, exons_data):
        for tid, tid_gp in groupby(exons_data, key=itemgetter(4)):
            num_of_exons = len(list(tid_gp))
            self.db[tid] = num_of_exons

    def get_num_of_exons(self, tids):
        return [self.db[tid] for tid in tids]



def get_longest_tid(tids, tid_len_dict, show_all=False):
    if tids == []:
        return []

    if show_all:
        return sorted(tids, key=lambda tid: [tid_len_dict.get(tid), tid], reverse=True)

    else:
        the_longest_len = max(map(tid_len_dict.get, tids))
        the_longest = sorted([tid for tid in tids if tid_len_dict.get(tid) == the_longest_len])
        the_longest = [the_longest[0]]

        return the_longest


def get_longest_common_tid(donor_tids, accepter_tids, tid_len_dict, show_all=False):
    common = set(donor_tids) & set(accepter_tids)
    if len(common) > 0:
        the_longest = get_longest_tid(common, tid_len_dict, show_all=show_all)
        return the_longest


def get_tid_exon_number(tid_data, tids):
    tid_data_dict = dict(tid_data)
    return list(map(tid_data_dict.get, tids))


def get_flanking_introns(intron_db, tids, exon_numbers, donor_accepter):
    return list(map(lambda tid, exon_number: \
                        intron_db.get_flanking_intron(tid, exon_number, donor_accepter), \
                    tids, exon_numbers))


def get_lens_of_intermediate_exons(exons_len_db, tids, start_exons, stop_exons):
    return list(map(lambda tid, en1, en2: \
                        exons_len_db.get_total_len_of_intermediate_exons(tid, en1, en2), \
                    tids, start_exons, stop_exons))


def print_results(results):
    def list_to_str(obj):
        if type(obj) == list:
            return ','.join(map(str, obj))
        else:
            return obj

    out_res = list(map(list_to_str, results))
    print(*out_res, sep='\t')


def print_usage():
    usage_msg = \
        "\n"\
        "  cat [NCLscan_result] | ./%(prog)s [Anno_gtf] > [Output_file]"\

    return usage_msg


def open_file(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, 'rt')
    else:
        return open(filename)


def create_parser():
    parser = argparse.ArgumentParser(usage=print_usage())
    parser.add_argument('anno_gtf', help='Annotation file. (".gtf" or ".gtf.gz")')
    parser.add_argument('--show-all', dest='show_all', action='store_true', help='Show all transcripts that satisfied the criteria.')
    parser.add_argument('--expand', action='store_true', help='Show all with expand mode.')
    parser.add_argument('--detail', action='store_true', help='Show more detail infomation of transcripts.')
    parser.add_argument('-V', '--version', action='version', version='{}'.format(__version__))

    return parser


if __name__ == "__main__":
    parser = create_parser()
    args = parser.parse_args()
    show_all = args.show_all
    anno_file = args.anno_gtf


    # generate db
    ## exons
    exon_extractor = ExonExtractor()
    with open_file(anno_file) as anno_data:
        exon_extractor.extract(anno_data)

    junc_sites_db = JunctionSitesDB()
    junc_sites_db.generate_db(exon_extractor.exons)

    ## introns
    intron_extractor = IntronExtractor()
    intron_extractor.extract(exon_extractor.exons)

    intron_db = FlankingIntronDB()
    intron_db.generate_db(intron_extractor.introns)

    ## intermediate exons
    exons_len_db = ExonsLengthDB()
    exons_len_db.generate_db(exon_extractor.exons)

    ## total number of exons
    num_of_exons_db = TranscriptExonsDB()
    num_of_exons_db.generate_db(exon_extractor.exons)


    # get exon numbers
    for line in sys.stdin:
        data = line.rstrip('\n').split('\t')
        ncl_event = NCLevent(data)
        donor = junc_sites_db.get_junc_site(ncl_event.donor, 'donor')
        accepter = junc_sites_db.get_junc_site(ncl_event.accepter, 'accepter')
        donor_tids = list(map(itemgetter(0), donor))
        accepter_tids = list(map(itemgetter(0), accepter))

        res_data = []
        len_of_intermediate_exons = []

        the_longest_common_tid = get_longest_common_tid(donor_tids, accepter_tids,
                                                        exon_extractor.transcripts,
                                                        show_all=show_all)
        if the_longest_common_tid:
            exon_number_donor = get_tid_exon_number(donor, the_longest_common_tid)
            exon_number_accepter = get_tid_exon_number(accepter, the_longest_common_tid)
            res_data += [the_longest_common_tid,
                            the_longest_common_tid,
                            exon_number_donor,
                            exon_number_accepter]

            # intermediate exons
            len_of_intermediate_exons = get_lens_of_intermediate_exons(
                                            exons_len_db,
                                            the_longest_common_tid,
                                            exon_number_accepter,
                                            exon_number_donor)
        else:
            the_longest_tid_donor = get_longest_tid(donor_tids, exon_extractor.transcripts,
                                                    show_all=show_all)
            the_longest_tid_accepter = get_longest_tid(accepter_tids,
                                                        exon_extractor.transcripts,
                                                        show_all=show_all)

            exon_number_donor = get_tid_exon_number(donor, the_longest_tid_donor)
            exon_number_accepter = get_tid_exon_number(accepter, the_longest_tid_accepter)

            res_data += [the_longest_tid_donor,
                         the_longest_tid_accepter,
                         exon_number_donor,
                         exon_number_accepter]

        # lens of flanking introns
        flanking_intron_len_donor = get_flanking_introns(intron_db,
                                                         res_data[0],
                                                         res_data[2],
                                                         'donor')
        flanking_intron_len_accepter = get_flanking_introns(intron_db,
                                                            res_data[1],
                                                            res_data[3],
                                                            'accepter')
        res_data += [flanking_intron_len_donor, flanking_intron_len_accepter]


        # intermediate exons
        res_data += [len_of_intermediate_exons]

        # total number of exons
        res_data += [num_of_exons_db.get_num_of_exons(res_data[0]),
                     num_of_exons_db.get_num_of_exons(res_data[1])]

        # additional transcript info when "--detail" is toggled
        if args.detail:
            if len(res_data[0]) > 0:
                is_protein_coding_donor, transcript_len_donor = map(list, zip(*[exon_extractor.transcripts[tid] for tid in res_data[0]]))
            else:
                is_protein_coding_donor, transcript_len_donor = [], []

            if len(res_data[1]) > 0:
                is_protein_coding_acceptor, transcript_len_acceptor = map(list, zip(*[exon_extractor.transcripts[tid] for tid in res_data[1]]))
            else:
                is_protein_coding_acceptor, transcript_len_acceptor = [], []

            res_data += [is_protein_coding_donor, is_protein_coding_acceptor, transcript_len_donor, transcript_len_acceptor]

        # print results
        if show_all and args.expand:
            for res in zip(*res_data):
                res = ncl_event.raw_data + list(res)
                print_results(res)
        else:
            res_data = ncl_event.raw_data + res_data
            print_results(res_data)
