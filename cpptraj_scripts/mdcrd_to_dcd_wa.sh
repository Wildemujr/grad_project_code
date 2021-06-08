#!/bin/bash

##############################################################################

date
read -p $'\nWhat is the name of the coordinate file (...wa.dyn.mdcrd) you are dealing with?\n> ' na
read -p $'\nHow many residues do you want to read into cpptraj?\n> ' resi
read -p $'\nWould you like to specify a specific frame range? [Y/n]\n> ' choice

# To easier SCP contents to a remote system.
user="$USER"
ip=$(hostname -I | awk '{print $1}')
cwd=$(pwd)
address="${user}@${ip}:${cwd}"

#To remove trailing white spaces in the variable.
address_clean=$(echo $address | sed 's/ *$//g')

# Creating temporary backups just in case.
cp $na.top x.top
cp $na.crd x.crd

if [ $choice == "Y" ]
then
	read -p $'\nWhat frame would you like to start with?\n> ' frm_start
	read -p $'\nWhat frame would you like to end with?\n> ' frm_end
	echo -e "\n$(date)\n Converting ${na}.wa.dyn.mdcrd to ${na}.nowa.dyn.dcd..."
cat <<eof > mdcrd_to_dcd_wa
parm ${na}.wa.top
trajin ${na}.wa.dyn.mdcrd ${frm_start} ${frm_end}
rms fit :1-${resi}
trajout ${na}.wa.dyn.dcd dcd
run
eof
	cpptraj x.top < in > out
	rm -f in out
	echo -e "\n$(date)\nCompleted the conversion!\n"
	echo "Use this command on your mac to transfer it: "
	echo -e " 'scp ${address_clean}/${na}.wa.dyn.dcd .' \n "
	# rm mdcrd_to_dcd_wa
	
else
	read -p $'How many frames make up the trajectory?\n> ' frm
	echo -e "\n$(date)\n Converting ${na}.wa.dyn.mdcrd to ${na}.nowa.dyn.dcd..."
cat <<eof > mdcrd_to_dcd_wa
parm ${na}.wa.top
trajin ${na}.wa.dyn.mdcrd
rms fit :1-${resi}
trajout ${na}.wa.dyn.dcd dcd
run
eof
	cpptraj x.top < in > out
	rm -f in out
	echo -e "\n$(date)\nCompleted the conversion!\n"
	echo "Use this command on your mac to transfer it: "
	echo -e " 'scp ${address_clean}/${na}.wa.dyn.dcd .' \n "
	# rm mdcrd_to_dcd_wa
fi 