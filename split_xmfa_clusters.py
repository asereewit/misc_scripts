#!/usr/bin/env python3

import os
import sys
import argparse
from Bio import SeqIO

def main(arguments):

    parser = argparse.ArgumentParser(
        description='Split parsnp xmfa file into unaligned fasta files grouped by cluster number',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('infile', help="Input parsnp xmfa file", type=str)

    args = parser.parse_args(arguments)

    cluster_dict={}
    for record in SeqIO.parse(args.infile,"fasta"):
        cluster_num = record.description.split(" ")[2]
        # get rid of "=" sign that is used to delimit clusters in xmfa file
        record.seq = record.seq.replace("=","")
        # get rid of "-" to produce unaligned fasta file; comment out to produce aligned fasta file
        record.seq = record.seq.replace("-","")
        if cluster_num in cluster_dict:
            cluster_dict[cluster_num].append(record)
        else:
            cluster_dict[cluster_num] = []
            cluster_dict[cluster_num].append(record)

    for key in cluster_dict:
        SeqIO.write(cluster_dict[key], f"{key}.fasta", "fasta")
    
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
