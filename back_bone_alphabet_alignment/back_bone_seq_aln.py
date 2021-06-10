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
                seqs_in_file.append(line.strip("-\n"))
    return headers_in_file, seqs_in_file


def convert_to_class_type(fasta_seq):
    temp = []
    for i in fasta_seq:
        if i in BackBone.back_bone_dict.keys():
            temp.append(f"{BackBone.back_bone_dict[i]}")
        else:
            continue
    return "".join(temp)


def convert_to_hex(fasta_seq):
    hex_li = [] 
    for i in fasta_seq:
        l = [format(ord(f"{c}"), "x") for c in i]
        hex_li.append(" ".join(l))
    return hex_li



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
    hex_seqs = convert_to_hex(fasta_seqs)
    z = list(zip(fasta_headers, hex_seqs))

    f = open("backbone_seqs.hex","a+")
    for i in z:
        f.write(f">{i[0]}\n")
        f.write(f"{i[1]}\n")
    f.close()






    # print(convert_to_class_type(fasta_seqs))

    # f = open(f"backbone_alphabet_seqs.faa","a+")
    # f.close()
    # Example case to see if it works
    X = "-ρ-β-β-ρρρ-ρρ-αiα-ρiαααρα-ρρα-αiρ-β-β-β-ααπ-β-ρ-β-β-βρρπ-ρα-αiα-ρiαααρρ-βρα-αiα-β-β-ρ"
    Y = "-ρβ-β-αρρ-ρα-αiα-ρiρααρα-βρα-αiα-β-β-ρ-β-β-βi-β-ρ-β-β-βρρπ-ρρ-αiα-ρiαααρρ-βρα-αiα-ρ-β-β"

    # X = "ptblyypaqaryaayabyaqabbpbbfbpbbbyyipyqaraaayybyaqapbb"
    # Y = "pbbpyypyqaraaayapyaqybbblaibpbbbyyipaqaraaayybyaqabbp"

    alignments = pairwise2.align.globalxx(X, Y)
    # l1 = convert_to_class_type(alignments[0][0])
    # l2 = convert_to_class_type(alignments[0][1])

    # for a in alignments:
    #     print(format_alignment(*a))


# ptbl--yypa-qaryaa-yab-yaqab-bpbbf---bpbbbyyipy-qaraaayybyaqapbb-
# | |   |||  ||| || ||  |||   | ||    |||||||||  ||||||||||||| || 
# p-b-bpyyp-yqar-aaaya-pyaq--yb-bb-laibpbbbyyip-aqaraaayybyaqa-bbp
#   Score=42