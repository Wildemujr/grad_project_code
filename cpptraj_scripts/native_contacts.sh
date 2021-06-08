#!/bin/bash

#############################################################################

read -p $'\nWhat is the name of the system that you\'re working with?\n> ' na
read -p $'\nHow many residues are in the system that you\'re working with?\n> ' resi

topo_wa=${na}.wa.top
coord_wa=${na}.wa.dyn.mdcrd
user="$USER"
ip=$(hostname -I | awk '{print $1}')
cwd=$(pwd)
address_clean="${user}@${ip}:${cwd}"
cp ${topo_wa} x.top
cp ${coord_wa} x.crd
##############################################################################
echo -e "\n$(date)\nRunning the native contact analysis..."

cat <<eof > in
parm ${topo_wa}
trajin ${coord_wa}
reference ${na}.wa.pdb
nativecontacts name NC1 :1-${resi} byresidue out nc.all.res.dat mindist maxdist distance 3.0 reference map mapout ${na}.resmap.gnu contactpdb ${na}.native.contacts.pdb series seriesout native.dat
run
eof

cpptraj x.top < in > out
rm -f in out
echo -e "\n$(date)\nCompleted computing native contacts!\n"

##############################################################################
echo "Use this command on your local machine to secure copy the output files:"
echo -e " 'scp ${address_clean}/{nonnative.${na}.resmap.gnu,native.${na}.resmap.gnu,nc.all.res.dat,native.dat,\
${na}.native.contacts.pdb} .'\n" 