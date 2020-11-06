#!/usr/bin/env python
#coding:utf-8 -*-

import os
import sys
import logging
import argparse
import pysam

LOG = logging.getLogger(__name__)

__author__ = ("Yating Zhu",)
__email__ = "1404213815qq.com"
__version__ = "v1.0"


def stat_n50(lenths):

    sumlen = sum(lenths)
    nsum = 0

    for n in lenths:
        nsum += n
        if nsum >= (sumlen/2):
            n50 = n
            break
    return n50


def stat_result(lenths):

    reads_num = len(lenths)
    total_num = sum(lenths)
    n50 = stat_n50(lenths)
    mid_len = total_num*1.0/reads_num
    max_readslen = max(lenths)

    return total_num, reads_num, n50, mid_len, max_readslen


def pick_bam(file, qvalue, minlen, name1, name2, types="bam"):
   
    raw_data = []
    clean_data = []
    fs = pysam.AlignmentFile(file, 'rb', check_sq=False)
    if types == "bam":
        f1 = pysam.AlignmentFile("%s.bam" % name1, "wb", template=fs)
    else:
        f1 = open("%s.fq" % name1, "w")

    for line in fs:
        raw_data.append(len(line.seq))
        if len(line.seq) < minlen or line.get_tag('rq') < qvalue:
            continue
        clean_data.append(len(line.seq))

        if types != "bam":
            f1.write("@%s\n%s\n%s\n%s\n" % (line.qname, line.seq, "+", line.qual))
            continue
        f1.write(line)

    fs.close()
    f1.close()
    
    f2 = open("%s.tsv" % name2, "w")    

    total_num, reads_num, n50, mid_len, max_readslen = stat_result(raw_data)
    f2.write('#Total_base(bp)\tReads_number(bp)\tN50(bp)\tMax_lenth(bp)\tMidlenth(bp)\n')
    f2.write('Raw\t{0:,}\t{1:,}\t{2:,.2f}\t{3:,.2f}\t{4:,}\n'.format(total_num, reads_num, n50, mid_len, max_readslen))
    total_num, reads_num, n50, mid_len, max_readslen = stat_result(clean_data)
    f2.write('Clean\t{0:,}\t{1:,}\t{2:,.2f}\t{3:,.2f}\t{4:,}\n'.format(total_num, reads_num, n50, mid_len, max_readslen))

    f2.close()
    

def add_help_args(parser):

    parser.add_argument('-b', '--bam', metavar = ('file'),
        help = "Read bam file")
    parser.add_argument('-q','--qvalue',  metavar = ('number'), type = float,
        help = "tag_rg in bam file")
    parser.add_argument('-l','--lenth',  metavar = ('number'), type = float,
        help = "lenth of reads seq in bam file")
    parser.add_argument('-n1','--name1',  metavar = ('filename'),
        help = "Output the name of bamfile")
    parser.add_argument('-n2','--name2',  metavar = ('filename'),
        help = "Output the name of stat_data tsvfile")
    parser.add_argument('-t','--types',  choices = ["bam", "fastq"], default = "bam",
        help = "Output the name of stat_data tsvfile")
    return parser


def main():

    logging.basicConfig(
        stream = sys.stderr,
        level = logging.INFO,
        format = "[%(levelname)s] % (message)s"
    )
    parser = argparse.ArgumentParser(formatter_class = argparse.RawDescriptionHelpFormatter,
    description = '''
name:
    gene_seq.py Pick seq from gff.
    
attention:
    gene_seq.py --

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_help_args(parser).parse_args()

    pick_bam(args.bam, args.qvalue, args.lenth, args.name1, args.name2, args.types)


if __name__ == "__main__":

    main()
