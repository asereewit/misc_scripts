#!/usr/bin/env python3

import os
import sys
import argparse
import pandas as pd
import itertools

def main(arguments):

    parser = argparse.ArgumentParser(
        description='Find the snp distances of all sequences with the same MLST value',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('mlst', help="MLST excel output file", type=str)
    parser.add_argument('snpdists', help="snp-dists excel output file", type=str)

    args = parser.parse_args(arguments)

    # group fasta sequences with the same mlst value in a dictionary
    mlst_dict = {}
    df_mlst = pd.read_excel(args.mlst, usecols=[0,2], names=['fasta','mlst_value'])
    for index, row in df_mlst.iterrows():
        if row['mlst_value'] not in mlst_dict:
            mlst_dict[row['mlst_value']] = []
        mlst_dict[row['mlst_value']].append(row['fasta'])

    # create two-element combinations of fasta sequences with the same mlst value
    fasta_combination = {}
    for key in mlst_dict:
        if len(mlst_dict[key]) > 1:
            if key not in fasta_combination:
                fasta_combination[key] = []
            fasta_combination[key].extend(list(itertools.combinations(mlst_dict[key],2)))

    # print snp-dists of each fasta combination
    df_snpdists = pd.read_excel(args.snpdists, index_col=0)
    mlst_snpdists = []
    for key in fasta_combination:    
        for pair in fasta_combination[key]:
            mlst_snpdists.append([key, pair[0], pair[1], df_snpdists.loc[pair[0],pair[1]]])

    df_out = pd.DataFrame(mlst_snpdists, columns=['mlst_value','fasta1','fasta2','snp_dists'])
    df_out.to_csv("mlst_snpdists.csv",index=False)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
