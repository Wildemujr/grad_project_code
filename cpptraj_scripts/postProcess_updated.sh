#!/bin/bash

##############################################################################
# Gaining user input and setting variables
##############################################################################

read -e -p $'\nWhat is name of your topology file?\n>  ' topo
read -e -p $'\nWhat is the name of your coordinate file?\n>  ' traj

read -e -p $'\nIs there particular frame range that you would like to read in? [y/n]\n>  ' frame_q
if [[ $frame_q == y || $frame_q == yes || $frame_q == Yes ]]
then
	read -e -p $'\n\tWhat is the range? Write the starting frame followed by a space and \n\tthe final frame. E.g., 2500 3500\n\t>  ' frames 
else
	echo -e "\n\t** Using all of the frames in the trajectory **\n"
	frames="1 last"
fi 

read -e -p $'\nAre you dealing with a large system? E.g., a hexadecamer? [y/n]\n>  ' big_sys_q
if [[ $big_sys_q == y || $big_sys_q == yes || $big_sys_q == Yes ]]
then
	read -e -p $'\n\tHow many residues comprise the original molecule prior to generation of symmetry mates?\n\t>  ' res_num
	res_num=":1-${res_num}"
else
	echo -e "\n\t ** Using all residues ** \n"
	res_num="*"
fi 

# user="$USER"
# ip=$(hostname -I | awk '{print $1}')
# cwd=$(pwd)
# address_clean="${user}@${ip}:${cwd}"

##############################################################################
echo -e "\n$(date)\nBeginning post-processing analysis..."

if [[ $big_sys_q == y || $big_sys_q == yes || $big_sys_q == Yes ]]
then 
cat <<eof > in_1
parm ${topo}
trajin ${traj} ${frames}
autoimage familiar origin anchor :${res_num} 
rms fit :${res_num}&!@H= mass
trajout ${traj%%.wa.dyn.mdcrd}.corr.wa.nc parm
run
eof
else
cat <<eof > in_1
parm ${topo}
trajin ${traj} ${frames}
autoimage
rms fit :${res_num}&!@H= mass
trajout ${traj%%.wa.dyn.mdcrd}.corr.wa.nc parm
run
eof
fi 
cpptraj ${topo} < in_1 > out
rm -f in_1 out
echo -e "Completed post-processing!\n$(date)\n"

#############################################################################
echo -e "\nExtracting PDBs from AMBER trajectory & renaming files..."
cat <<eof > in_2
parm ${topo}
trajin ${traj%%.wa.dyn.mdcrd}.corr.wa.nc
trajout ${traj%%.wa.dyn.mdcrd}.corr.nowa.nc parm
strip :WAT
run
eof
cpptraj ${topo} < in_2 > out
rm -f in_2 out
echo -e "\nDone Extracting PDBs from AMBER trajectory & renaming files!\n$(date)\n"

#############################################################################
echo -e "\nStripping Waters from Topology File..."
cat <<eof > in_3
parm ${topo}
parmstrip :WAT
parmwrite out ${topo%%.wa.top}.nowa.top
run
eof
cpptraj ${topo} < in_3 > out
rm -f in_3 out
echo -e "\nDone Stripping Waters from Topology!\n$(date)\n"
#############################################################################
echo -e "\nCalculating the RMSD..."
cat <<eof > in_4
parm ${topo%%.wa.top}.nowa.top
trajin ${traj%%.wa.dyn.mdcrd}.corr.nowa.nc
rms RMSD :*&!@H= first out ${traj%%.wa.dyn.mdcrd}_rmsd.dat mass
run
eof
cpptraj "${topo%%.wa.top}.nowa.top" < in_4 > out
rm -f in_4 out
echo -e "\nDone calculating the RMSD!\n$(date)\n"
#############################################################################
echo -e "\nCalculating the RMSF..."
cat <<eof > in_5
parm ${topo%%.wa.top}.nowa.top
trajin ${traj%%.wa.dyn.mdcrd}.corr.nowa.nc 
rms first
average crdset MyAvg
run
rms ref MyAvg
atomicfluct out ${traj%%.wa.dyn.mdcrd}_rmsf.dat byres 
run
eof
cpptraj "${topo%%.wa.top}.nowa.top" < in_5 > out
rm -f in_5 out
echo -e "\nDone calculating the RMSF!\n$(date)\n"




