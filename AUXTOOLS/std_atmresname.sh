#!/bin/bash

NARG=$#
declare -a ARG=""
ARG=($(echo $*))

if [ $NARG -lt 2  ] 
then
    echo "Usage: std_atmresname.sh [ -atm | -res | -atmres ]  pdb1 pdb2 ... "
    exit
fi
OPT=${ARG[0]}

IATM=0
IRES=0
if [ $OPT == '-atm' ]; then IATM=1; fi
if [ $OPT == '-res' ]; then IRES=1; fi
if [ $OPT == '-atmres' ]; then IATM=1; IRES=1; fi
if [ $OPT == '-resatm' ]; then IATM=1; IRES=1; fi

if [ $IATM -eq 0 ] && [ $IRES -eq 0 ]
then
    echo "First argument: $OPT is not correct!"
    echo "Usage: std_atmresname.sh [ -atm | -res | -atmres ]  pdb1 pdb2 ... "
    exit
fi

#
#  AMBER residue/atom names are fixed here 
#
if [ $IRES -eq 1 ]
then 
	for ((I=2;I<=NARG;I++))
	do
		let "J=$I-1"
		pdb=${ARG["$J"]} 
		sed  's/HIE/HIS/' $pdb |sed  's/HID/HIS/' |sed  's/GLH/GLU/' |sed  's/ASH/ASP/' |sed  's/HIP/HIS/' \
		|sed  's/ARN/ARG/' |sed  's/LYN/LYS/' |sed  's/CYM/CYS/' |sed  's/CYX/CYS/' \
		|sed  's/DA3/ DA/' |sed  's/DC3/ DC/' |sed  's/DG3/ DG/' |sed  's/DT3/ DT/' |sed  's/DA5/ DA/' \
		|sed  's/DC5/ DC/' |sed  's/DG5/ DG/' |sed  's/DT5/ DT/' >tmpstdresname_${pdb}
		mv -f tmpstdresname_${pdb} $pdb
		sed -i 's/ W.. / WAT /' $pdb
	done
fi

if [ $IATM -eq 1 ]
then
	for ((I=2;I<=NARG;I++))
	do
		let "J=$I-1"
		pdb=${ARG["$J"]} 
		sed  's/ 1HB / HB1 /'  $pdb | sed  's/ 2HB / HB2 /'  | sed  's/ 3HB / HB3 /'  | sed  's/ 2HA / HA2 /'  \
		| sed  's/ 3HA / HA3 /'   | sed  's/ 2HG / HG2 /'   | sed  's/ 3HG / HG3 /'   | sed  's/ 2HD / HD2 /'  \
		| sed  's/ 3HD / HD3 /'   | sed  's/ 1HG / HG1 /'   | sed  's/ 2HE / HE2 /'   | sed  's/ 3HE / HE3 /'  \
		| sed  's/ 1HZ / HZ1 /'   | sed  's/ 2HZ / HZ2 /'   | sed  's/ 3HZ / HZ3 /'   | sed  's/ 1HE / HE1 /'  \
		| sed  's/ 2HE / HE2 /'   | sed  's/ 3HE / HE3 /'   | sed  's/ 2HD / HD2 /'   | sed  's/ 3HD / HD3 /'  \
		| sed  's/ 1HH3 / HH31 /' | sed  's/ 2HH3 / HH32 /' | sed  's/ 3HH3 / HH33 /' | sed  's/ 1HD1 / HD11 /' \
		| sed  's/ 2HD1 / HD12 /' | sed  's/ 3HD1 / HD13 /' | sed  's/ 1HD2 / HD21 /' | sed  's/ 2HD2 / HD22 /' \
		| sed  's/ 3HD2 / HD23 /' | sed  's/ 1HH1 / HH11 /' | sed  's/ 2HH1 / HH12 /'  | sed  's/ 1HH2 / HH21 /'  \
		| sed  's/ 2HH2 / HH22 /'  | sed  's/ 1HE2 / HE21 /'  | sed  's/ 2HE2 / HE22 /'  | sed  's/ 1HG1 / HG11 /'  \
		| sed  's/ 2HG1 / HG12 /'  | sed  's/ 3HG1 / HG13 /'  | sed  's/ 1HG2 / HG21 /'  | sed  's/ 2HG2 / HG22 /'  \
		| sed  's/ 3HG2 / HG23 /'  | sed  's/ 1HD1 / HD11 /'   > tmpstdresname_${pdb}
		mv -f tmpstdresname_${pdb} $pdb
	done
fi

