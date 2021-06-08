#!/bin/bash

##############################################################################

read -p $'\nWhat is the name of the system you are dealing with?\n> ' resi
na=${resi}.wa

# If you want to use the directory name to be the new file name
# full_path=$(realpath $0)
# dir_path=$(dirname $full_path)
# folder_name=${dir_path##*/}
# new_file="${na}.${folder_name}"
# echo "${na}"

user="$USER"
ip=$(hostname -I | awk '{print $1}')
cwd=$(pwd)
address_clean="${user}@${ip}:${cwd}"

cp ${na}.top x.top
cp ${na}.crd x.crd

##############################################################################
echo -e '\n$(date)\nRunning RMSF and B-factor calculations...'
cat <<eof> traj_test
clear all
parm ${na}.top
loadcrd ${na}.corr.mdcrd
crdaction ${na}.corr.mdcrd average ${na}.avg.pdb
crdaction ${na}.corr.mdcrd average ${na}.avg_CA.pdb @CA
parm ${na}.avg_CA.pdb
reference ${na}.avg_CA.pdb parm ${na}.avg_CA.pdb
crdaction ${na}.corr.mdcrd rms reference @CA
crdaction ${na}.corr.mdcrd atomicfluct out ${na}.avg.rmsf_fluct.dat @CA
crdaction ${na}.corr.mdcrd atomicfluct out ${na}.avg.bfactor_fluct.dat bfactor @CA

eof

# cpptraj -i traj_test; echo -e "\n"
cpptraj x.top < traj_test > out
echo -e "\n$(date)\nCompleted RMS and B-factor calculations!\n"
echo "Use this command on your local machine to secure copy ${na}_ensemble.pdb"
echo -e " 'scp ${address_clean}/{${na}.avg.pdb,${na}.avg.rmsf_fluct.dat} .' \n " 
rm traj_test

##############################################################################
