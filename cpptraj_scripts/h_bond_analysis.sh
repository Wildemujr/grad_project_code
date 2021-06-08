#!/bin/bash

##############################################################################

read -p $'\nWhat is the name of the system that you\'re working with?\n> ' na
read -p $'\nHow many residues are there in the system?\n> ' resi

topo_wa=${na}.wa.top
coord_wa=${na}.wa.dyn.mdcrd
user="$USER"
ip=$(hostname -I | awk '{print $1}')
cwd=$(pwd)
address_clean="${user}@${ip}:${cwd}"
cp ${topo_wa} x.top
cp ${coord_wa} x.crd
##############################################################################
echo -e "\n$(date)\nRunning the analysis..."


cat <<eof > in
clear all
parm ${topo_wa}
trajin ${coord_wa}
hbond All out All.hbvtime.dat solventdonor :WAT solventacceptor :WAT@O avgout All.UU.avg.dat solvout All.UV.avg.dat bridgeout All.bridge.avg.dat         # First line is to track solute-solute (UU) and solute-solvent (UV) H-bonds. Solute-solvent-solute (bridge) bridging interactions are also tracked.
hbond Backbone :1-${resi}@C,O,N,H avgout BB.avg.dat series uuseries bbhbond.gnu     									 # This tracks solute-solute backbone hydrogens bonds. uuseries writes the raw hydrogen bond time series to bbhbond.gnu.
create nhbvtime.agr All[UU] Backbone[UU] All[UV] All[Bridge]                        									 # Will create a file containing only the number of solute-solute hydrogen bonds, and all other aspects as well. 
rms BBrmsd :1-${resi}@C,CA,N out BBrmsd.agr                                         								  	 # Add RMSD to calculate backbone RMSD to the first frame to track unfolding.
run
lifetime Backbone[solutehb] out backbone.lifetime.dat                               									 # Runs lifetime analysis on the backbone hydrogen bonds. crv.backbone.lifetime.dat is a way to visualize the time dependent behavior of H bonds.
runanalysis
eof

cpptraj x.top < in > out
rm -f in out
echo -e "\n$(date)\nCompleted the Hydrogen Bond Analysis!\n"

##############################################################################
echo "Use this command on your local machine to secure copy the output files:"
echo -e " 'scp ${address_clean}/{All.bridge.avg.dat,All.hbvtime.dat,All.UU.avg.dat,All.UV.avg.dat,BB.avg.dat,\
BBrmsd.agr,nhbvtime.agr,backbone.lifetime.dat,crv.backbone.lifetime.dat} .'\n" 
# rm traj_test
