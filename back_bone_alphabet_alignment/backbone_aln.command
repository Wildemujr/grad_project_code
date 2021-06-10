#!/bin/bash

# protSequence {whatever PDBS I end up putting here...Most likely will be fed first into biopython.} E.g., python3 obtain_pdb.py or whatever.

python3 back_bone_seq_aln.py -i bb.faa
/usr/local/libexec/mafft/hex2maffttext backbone_seqs.hex > input.ASCII
mafft --globalpair --maxiterate 1000 --text input.ASCII > output.ASCII
/usr/local/libexec/mafft/maffttext2hex output.ASCII > output.hex
awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' output.hex > lin_msa_ex.faa
mate lin_msa_ex.faa