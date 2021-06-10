#!/usr/bin/env python3
 
from pprint import pprint as pp
import argparse
import numpy as np
import pandas as pd
from backbone_classification_dict import *



def read_fasta_file(file_na):
    headers_in_file, seqs_in_file = [], []
    with open(f"{file_na}", "r") as fh:
        for line in fh.readlines():
            if ">" in line:
                headers_in_file.append(line[1:].rstrip())
            else:
                seqs_in_file.append(line.strip())
    return headers_in_file, seqs_in_file




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        dest="inputf",
        help="Choosing the Input File to be Processed and Reformatted."
    )
    # parser.add_argument(
        # "-o",
        # "--output",
        # dest="outputf",
        # help="Choosing the Name of the Ouput PDB File to be Output After Processing is Complete."
    # )
    option = parser.parse_args()
    aln = Alignment()

    fasta_headers, fasta_seqs = read_fasta_file(f"{option.inputf}")
    fasta_seqs = (i.split(" ") for i in fasta_seqs)
    greek_bb = (aln.hex_to_backbone_alph(seq) for seq in fasta_seqs)

    l = []     
    for i in greek_bb:
        l.append(list(i))

    for p,i in enumerate(l):
        for j,k in enumerate(i):
            if len(k)==1:
                i[j] = f"  {k}"
            if len(k)==2:
                i[j] = f" {k}"

    l = (list(map(lambda x:" ".join(x), l)))

    comb = list(zip(fasta_headers,l))
    f = open("greek_aln.faa","w+")
    for i in comb:
        f.write(f">{i[0]}\n")
        f.write(f'{i[1]}\n')
    f.close()