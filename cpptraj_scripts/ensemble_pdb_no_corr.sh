#!/bin/bash

##############################################################################
printf "\n"
date 
read -p $'What is the name of the system you are dealing with?\n> ' na
read -p $'\nHow many frames are there in the trajectory?\n> ' frames
read -p $'\nWhat should the frame offset be?\n> ' offset   
read -p $'\nHow many residues do you want to read into cpptraj?\n> ' resi

user="$USER"
ip=$(hostname -I | awk '{print $1}')
cwd=$(pwd)
address="${user}@${ip}:${cwd}"

#To remove trailing white spaces in the variable
address_clean=$(echo $address | sed 's/ *$//g')

cp $na.wa.top x.top
cp $na.wa.crd x.crd

##############################################################################
printf "\nCreating a protein ensemble with the details provided..."
# echo "Creating a protein ensemble with the details provided..."
cat <<eof> ensemble_file
clear all
parm ${na}.wa.top
trajin ${na}.wa.corr.mdcrd 1 ${frames} ${offset}
rms fit :1-${resi}
strip :WAT
trajout ${na}_ensemble.pdb pdb
run

eof

cpptraj x.top < ensemble_file > out
echo -e "\n"
date
echo -e "Completed creating the protein ensemble!\n"
echo "Use this command on your local machine to secure copy ${na}_ensemble.pdb"
echo -e " 'scp ${address_clean}/${na}_ensemble.pdb .' \n " 
# rm ensemble_file


##############################################################################




