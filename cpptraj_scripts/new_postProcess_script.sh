#!/bin/bash


read -e -p $'\nWhat is name of your topology file?\n>  ' topo
read -e -p $'\nWhat is the name of your coordinate file?\n>  ' traj
read -e -p $'\nHow many residues are there?\n>  ' res_num
read -e -p $'\nWhat is the residue number of the first cluster?\n>  ' clust1
read -e -p $'\nWhat is the residue number of the second cluster?\n>  ' clust2
read -e -p $'\nWhat run is this?\n>  ' run_num


read -e -p $'\nIs there particular frame range that you would like to read in? [y/n]\n>  ' frame_q
if [[ $frame_q == y || $frame_q == yes || $frame_q == Yes ]]
then
	read -e -p $'\n\tWhat is the range? Write the starting frame followed by a space and \n\tthe final frame. E.g., 2500 3500\n\t>  ' frames 
else
	echo -e "\n\t** Using all of the frames in the trajectory **\n"
	frames="1 last"
fi


user="$USER"
ip=$(hostname -I | awk '{print $1}')
cwd=$(pwd)
address_clean="${user}@${ip}:${cwd}"


#############################################################################
#############################################################################
echo -e "\n$(date)\nBeginning post-processing analysis..."
cat <<eof > in_1
parm ${topo}
trajin ${traj} ${frames}
autoimage
rms fit :${res_num}&!@H= mass
trajout ${traj%%.wa.dyn.mdcrd}.wa.nc parm
run
eof
#############################################################################
cpptraj ${topo} < in_1 > out
rm -f in_1 out
echo -e "\nDone with Post Processing Files!\n$(date)\n"
#############################################################################
echo -e "\nExtracting PDBs from AMBER trajectory & renaming files..."
cat <<eof > in_2
parm ${topo}
trajin ${traj%%.wa.dyn.mdcrd}.wa.nc
trajout ${traj%%.wa.dyn.mdcrd}.nowa.nc parm
strip :WAT
run
eof
cpptraj ${topo} < in_2 > out
rm -f in_2 out
echo -e "\nDone Extracting PDBs from AMBER trajectory & renaming files!\n$(date)\n"
#############################################################################
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
#############################################################################
echo -e "\nCalculating the Cluster-Cluster Distance..."
cat <<eof > in_4
parm ${topo%%.wa.top}.nowa.top
trajin ${traj%%.wa.dyn.mdcrd}.nowa.nc  ${frames}
distance :${clust1}@* :${clust2}@* out ${traj%%.wa.dyn.mdcrd}.clustDist.${run_num}.txt
run
eof
cpptraj "${topo%%.wa.top}.nowa.top" < in_4 > out
rm -f in_4 out
echo -e "\nDone calculating the Cluster-Cluster Distance!\n$(date)\n"
#############################################################################