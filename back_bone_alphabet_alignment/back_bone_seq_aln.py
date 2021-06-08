#!/usr/bin/bash python3

from pprint import pprint as pp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
from backbone_classification_dict import BackBone
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


def read_fasta_file(file_na):
    headers_in_file, seqs_in_file = [], []
    with open(f"{file_na}", "r") as fh:
        for line in fh.readlines():
            if ">" in line:
                headers_in_file.append(line[1:].rstrip())
            else:
                seqs_in_file.append(line.strip("-"))
    return headers_in_file, seqs_in_file


def convert_to_class_type(fasta_seq):
    temp = []
    for i in fasta_seq[1]:
        if i in BackBone.back_bone_dict.keys():
            temp.append(f"{BackBone.back_bone_dict[i]}")
    return "".join(temp)


if __name__ == "__main__":

    # Defining the arguments that will be able to be passed on the command line.
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        dest="inputf",
        help="Choosing the Input PDB File to be Processed and Reformatted.",
    )
    # parser.add_argument(
    #     "-o",
    #     "--output",
    #     dest="outputf",
    #     help="Choosing the Name of the Ouput PDB File to be Output After Processing is Complete.",
    # )

    option = parser.parse_args()

    fasta_headers, fasta_seqs = read_fasta_file(f"{option.inputf}")

    # print(convert_to_class_type(fasta_seqs))

    # f = open(f"backbone_alphabet_seqs.faa","a+")
    # f.close()
    # Example case to see if it works
    X = "-ρ-β-β-ρρρ-ρρ-αiα-ρiαααρα-ρρα-αiρ-β-β-β-ααπ-β-ρ-β-β-βρρπ-ρα-αiα-ρiαααρρ-βρα-αiα-β-β-ρ"
    Y = "-ρβ-β-αρρ-ρα-αiα-ρiρααρα-βρα-αiα-β-β-ρ-β-β-βi-β-ρ-β-β-βρρπ-ρρ-αiα-ρiαααρρ-βρα-αiα-ρ-β-β"

    alignments = pairwise2.align.globalxx(X, Y)
    for a in alignments:
        print(format_alignment(*a))
