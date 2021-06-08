#1/bin/bash

##############################################################################

date
read -p "What is the name of the coordinate file (...wa.dyn.mdcrd) you are dealing with?   " na
read -p "How many residues do you want to read into cpptraj?    " resi

user="$USER"
ip=$(hostname -I | awk '{print $1}')
cwd=$(pwd)
address="${user}@${ip}:${cwd}"

#To remove trailing white spaces in the variable
address_clean=$(echo $address | sed 's/ *$//g')

cp $na.top x.top
cp $na.crd x.crd

##############################################################################
echo "Converting "${na}".wa.dyn.mdcrd to "${na}".nowa.dyn.dcd..."
cat <<eof> mdcrd_to_dcd_nowa
parm ${na}.wa.top
trajin ${na}.wa.dyn.mdcrd
rms fit :1-${resi}
strip :WAT
trajout ${na}.nowa.dyn.dcd dcd
run
eof

cpptraj -i mdcrd_to_dcd_nowa
date
echo -e "Completed the conversion from mdcrd to dcd!\n"
echo "Use this command on your mac to transfer it: "
echo -e " 'scp ${address_clean}/${na}_nowa.dyn.dcd .' \n " 
rm mdcrd_to_dcd_nowa
##############################################################################


