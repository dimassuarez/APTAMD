#!/bin/bash 

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi

if [ "$#" -eq 0 ]; then more $APTAMD/DOC/do_system_edition.txt ; exit; fi

# do_system_edition: Starting from PDB file, this scripts performs the following tasks
#                    tleap edition adding counterions to satisfy the required ionic force

source $APTAMD/ENV/aptamd_env.sh

extension="${1##*.}"
if [ $extension == "src" ]
then
    echo "Sourcing $1"
    source $1
else
    INITIAL=$1
fi
if [ -z "$INITIAL" ]; then more $APTAMD/DOC/do_system_edition.txt ; exit; fi
if [ -z "$IONIC" ]; then IONIC="0.150"; echo 'Assuming IONIC STRENGTH=0.150 M' ; fi
if [ -z "$BUFFER_SOLV" ]; then BUFFER_SOLV=16.0; echo 'Assuming BUFFER SOLV= 16.0 A ' ; fi
if [ -z "$FF" ]; then FF=19; echo 'Assuming FF19SB / OPC ' ; fi
# GAFF
if [ -z $GAFF ]
then
   GAFF="gaff2"
fi
if [ $GAFF != "gaff2" ] && [ $GAFF != "gaff" ] 
then
   echo "GAFF=$GAFF, but it can only be gaff2 or gaff"
   exit
fi


if [ ! -e $INITIAL ]; then echo "$INITIAL does not exist in the current location"; exit; fi

if [ "$FF" -ne 19 ] &&  [ "$FF" -ne 14 ] ; then echo "FF=14 or 19, but FF=$FF"; exit ; fi

# Filenames
MOL=${INITIAL/_initial.pdb/}

if [ $MOL == $INITIAL ]; then echo "$INITIAL does not conform MOL_initial.pdb filename"; exit; fi

if [ -e ${MOL}.off ] 
then
   echo "An AMBER OFF library file ${MOL}.off is present in the edition directory"
   echo "It will be incorporated in the edition scripts"
   OFFLIB="YES"
else
   OFFLIB="NO"
fi
if [ -e ${MOL}.frcmod ] 
then
   echo "An AMBER Params File ${MOL}.frcmod is present in the edition directory"
   echo "It will be incorporated in the edition scripts"
   FRCMOD="YES"
else
   FRCMOD="NO"
fi
if [ -e ${MOL}_extra.src ] 
then
   echo "A LEaP extra src file ${MOL}_extra.src is present in the edition directory"
   echo "It will be incorporated into the edition scripts after loading coordinates"
   echo "Typically, this is needed to fix some missing detail in the parameterization"
   echo "Prepare carefully ${MOL}_extra.src and use only mol as unit name"
   EXTRA="YES"
else
   EXTRA="NO"
fi
if [ -z ${BELLYMASK} ]
then
   BELLYMASK="NONE"
else
   MAXCYC=500 
   echo "BELLYMASK=${BELLYMASK} --> Selected atoms will be relaxed " 
fi

INITMOL=$INITIAL
TOP_SOLUTE=${MOL}_solute.top
CRD_SOLUTE=${MOL}_solute.crd
PDB_SOLUTE=${MOL}_solute.pdb
TOP=${MOL}.top
CRD=${MOL}.crd
PDB=${MOL}.pdb

# Edition with tleap 
if [ $FF -eq 19 ]
then 

echo '# Force Field data' >edit_leap_solute.src 
echo 'source leaprc.DNA.bsc1'  >>edit_leap_solute.src
echo 'source leaprc.protein.ff19SB' >> edit_leap_solute.src
echo "source leaprc.${GAFF}" >> edit_leap_solute.src
if [ $OFFLIB == "YES" ]; then echo "loadoff ${MOL}.off" >> edit_leap_solute.src ; fi
if [ $FRCMOD == "YES" ]; then echo "loadAmberParams ${MOL}.frcmod" >> edit_leap_solute.src ; fi
echo 'WAT=OPC' >> edit_leap_solute.src
echo 'loadamberparams frcmod.ionslm_iod_opc' >>edit_leap_solute.src 
echo '# Build System' >>edit_leap_solute.src 
echo 'mol=loadpdb' $INITMOL >>edit_leap_solute.src 
if [ $EXTRA == "YES" ]; then cat ${MOL}_extra.src  >> edit_leap_solute.src ; fi
echo 'check mol ' >>edit_leap_solute.src 
echo 'charge mol ' >>edit_leap_solute.src 
echo 'saveamberparm mol ' $TOP_SOLUTE  $CRD_SOLUTE >>edit_leap_solute.src 
echo 'savepdb mol ' $PDB_SOLUTE >>edit_leap_solute.src 
echo 'quit' >>edit_leap_solute.src 

elif [ $FF -eq 14 ]
then

echo '# Force Field data' >edit_leap_solute.src 
echo 'source leaprc.DNA.bsc1'  >>edit_leap_solute.src
echo 'source leaprc.protein.ff14SB' >> edit_leap_solute.src
echo 'loadamberparams frcmod.ionsjc_tip3p' >>edit_leap_solute.src 
echo "source leaprc.${GAFF}" >> edit_leap_solute.src
if [ $OFFLIB == "YES" ]; then echo "loadoff ${MOL}.off" >> edit_leap_solute.src ; fi
if [ $FRCMOD == "YES" ]; then echo "loadAmberParams ${MOL}.frcmod" >> edit_leap_solute.src ; fi
echo '# Build System ' >>edit_leap_solute.src 
echo 'mol=loadpdb' $INITMOL >>edit_leap_solute.src 
if [ $EXTRA == "YES" ]; then cat ${MOL}_extra.src  >> edit_leap_solute.src ; fi
echo 'check mol ' >>edit_leap_solute.src 
echo 'charge mol ' >>edit_leap_solute.src 
echo 'saveamberparm mol ' $TOP_SOLUTE  $CRD_SOLUTE >>edit_leap_solute.src 
echo 'savepdb mol ' $PDB_SOLUTE >>edit_leap_solute.src 
echo 'quit' >>edit_leap_solute.src 

else

echo "FF not defined!"

fi

rm -f leap.log
echo  Editing solute coordinates
$AMBERHOME/bin/tleap -f edit_leap_solute.src > mlog 
mv leap.log  edit_leap_solute.log
tail -1 edit_leap_solute.log

#  Relaxation of selected residues 
if [ ${BELLYMASK}  != "NONE" ]
then

cat  <<EOF  > sander_relax_solute.inp 
Partial relaxation (PO4 backbone atom fixed)  Distance dependent eps
&cntrl
 imin=1, ncyc=100,  maxcyc=${MAXCYC},  ntmin=2, drms=0.02
 ntb=0, ntf=1, ntc=1, ntpr=500,
 cut=1000.0, nsnb=5000,  igb=0,
 ibelly=1,
 bellymask="${BELLYMASK}"
 /
 &ewald
  eedmeth=5,
 /
EOF

echo "Running sander to relax BELLYMASK atoms"

$AMBERHOME/bin/sander -O -i sander_relax_solute.inp -p $TOP_SOLUTE -c $CRD_SOLUTE -r sander_relax_solute.crd -o sander_relax_solute.out

echo "$CRD_SOLUTE contains the partially relaxed coordinates"
echo "$PDB_SOLUTE contains the partially relaxed coordinates"
mv sander_relax_solute.crd $CRD_SOLUTE

$AMBERHOME/bin/cpptraj $TOP_SOLUTE<<EOF
trajin $CRD_SOLUTE
autoimage
trajout ${PDB_SOLUTE} dumpq include_ep
go
EOF

fi

# Edition with tleap to add solvent box 
# First we edit without adding counter ions

if [ "$FF" -eq 14 ]
then 

echo '# Force Field data' >edit_leap.src 
echo 'source leaprc.DNA.bsc1'  >>edit_leap.src
echo 'source leaprc.protein.ff14SB'  >>edit_leap.src
echo "source leaprc.${GAFF}" >> edit_leap.src
echo 'source leaprc.water.tip3p'  >>edit_leap.src
echo 'loadamberparams frcmod.ionsjc_tip3p' >>edit_leap.src 
if [ $OFFLIB == "YES" ]; then echo "loadoff ${MOL}.off" >> edit_leap.src ; fi
if [ $FRCMOD == "YES" ]; then echo "loadAmberParams ${MOL}.frcmod" >> edit_leap.src ; fi
echo '# Build system' >>edit_leap.src 
echo 'mol=loadpdb' $PDB_SOLUTE >>edit_leap.src 
if [ $EXTRA == "YES" ]; then cat ${MOL}_extra.src  >> edit_leap.src ; fi
echo 'alignaxes mol ' >>edit_leap.src
echo 'solvateOct mol TIP3PBOX ' $BUFFER_SOLV >>edit_leap.src
echo 'saveamberparm mol ' $TOP  $CRD >>edit_leap.src 
echo 'savepdb mol ' $PDB >>edit_leap.src 
echo 'quit' >>edit_leap.src 

elif [ $FF -eq 19 ]
then

echo '# Force Field data' >edit_leap.src
echo 'source leaprc.DNA.bsc1'  >>edit_leap.src
echo 'source leaprc.protein.ff19SB' >> edit_leap.src
echo "source leaprc.${GAFF}" >> edit_leap.src
echo 'WAT=OPC' >> edit_leap.src
echo 'source leaprc.water.opc' >> edit_leap.src
echo 'loadamberparams frcmod.ionslm_iod_opc' >>edit_leap.src 
if [ $OFFLIB == "YES" ]; then echo "loadoff ${MOL}.off" >> edit_leap.src ; fi
if [ $FRCMOD == "YES" ]; then echo "loadAmberParams ${MOL}.frcmod" >> edit_leap.src ; fi
echo '# Build system' >>edit_leap.src 
echo 'mol=loadpdb' $PDB_SOLUTE >>edit_leap.src 
if [ $EXTRA == "YES" ]; then cat ${MOL}_extra.src  >> edit_leap.src ; fi
echo 'alignaxes mol ' >>edit_leap.src
echo "solvateOct mol OPCBOX  $BUFFER_SOLV" >>edit_leap.src
echo 'saveamberparm mol ' $TOP  $CRD >>edit_leap.src 
echo 'savepdb mol ' $PDB >>edit_leap.src 
echo 'quit' >>edit_leap.src 

fi

echo "Adding solvent box with BUFFER SOLV= $BUFFER_SOLV" 
rm -f leap.log
$AMBERHOME/bin/tleap -f edit_leap.src > mlog 
mv leap.log  edit_leap.log
tail -1 edit_leap.log
rm -f mlog

NWAT=$(grep -c 'O   WAT' $PDB)
Q=$(grep 'Total unperturbed charge:' edit_leap_solute.log | awk '{print $4}')

echo "System contains $NWAT waters. Solute charge = $Q. Ionic strenght= $IONIC (M) " 
$OCTAVE -q  <<EOF  > mlog
q=round( ${Q} );
if q < 0 
  num_Na = round( ( ${IONIC} / 55 )* ${NWAT});
  num_Cl = num_Na - abs(q) ; 
else
  num_Cl = round( ( ${IONIC} / 55 )* ${NWAT});
  num_Na = num_Cl - abs(q) ; 
end
disp([ 'Na= ',num2str(num_Na)])
disp([ 'Cl= ',num2str(num_Cl)])
EOF

cat mlog 

NUM_NA=$(head -1 mlog  | awk '{print $2}')
NUM_CL=$(tail -1 mlog  | awk '{print $2}')

if [ $NUM_NA -gt 0 ]  ||  [ $NUM_CL -gt  0 ]    # Only if counterions are needed !
then 

if [ "$FF" -eq 14 ]
then 

echo '# Force Field data' >edit_leap.src 
echo 'source leaprc.DNA.bsc1'  >>edit_leap.src
echo 'source leaprc.protein.ff14SB'  >>edit_leap.src
echo "source leaprc.${GAFF}" >> edit_leap.src
echo 'source leaprc.water.tip3p'  >>edit_leap.src
echo 'loadamberparams frcmod.ionsjc_tip3p' >>edit_leap.src 
if [ $OFFLIB == "YES" ]; then echo "loadoff ${MOL}.off" >> edit_leap.src ; fi 
if [ $FRCMOD == "YES" ]; then echo "loadAmberParams ${MOL}.frcmod" >> edit_leap.src ; fi
echo '# Build System' >>edit_leap.src 
echo 'mol=loadpdb' $PDB_SOLUTE >>edit_leap.src 
if [ $EXTRA == "YES" ]; then cat ${MOL}_extra.src  >> edit_leap.src ; fi
echo 'alignaxes mol ' >>edit_leap.src
echo 'solvateOct mol TIP3PBOX ' $BUFFER_SOLV >>edit_leap.src
echo "addionsrand mol  Na+ $NUM_NA  Cl- $NUM_CL " >>edit_leap.src
echo 'saveamberparm mol ' $TOP  $CRD >>edit_leap.src 
echo 'savepdb mol ' $PDB >>edit_leap.src 
echo 'quit' >>edit_leap.src 

elif [ $FF -eq 19 ]
then

echo '# Force Field data' >edit_leap.src
echo 'source leaprc.DNA.bsc1'  >>edit_leap.src
echo 'source leaprc.protein.ff19SB' >> edit_leap.src
echo "source leaprc.${GAFF}" >> edit_leap.src
echo 'WAT=OPC' >> edit_leap.src
echo 'source leaprc.water.opc' >> edit_leap.src
echo 'loadamberparams frcmod.ionslm_iod_opc' >>edit_leap.src 
if [ $OFFLIB == "YES" ]; then echo "loadoff ${MOL}.off" >> edit_leap.src ; fi
if [ $FRCMOD == "YES" ]; then echo "loadAmberParams ${MOL}.frcmod" >> edit_leap.src ; fi
echo '# Build system' >>edit_leap.src 
echo 'mol=loadpdb' $PDB_SOLUTE >>edit_leap.src 
if [ $EXTRA == "YES" ]; then cat ${MOL}_extra.src  >> edit_leap.src ; fi
echo 'alignaxes mol ' >>edit_leap.src
echo "solvateOct mol OPCBOX  $BUFFER_SOLV" >>edit_leap.src
echo "addionsrand mol  Na+ $NUM_NA  Cl- $NUM_CL " >>edit_leap.src
echo 'saveamberparm mol ' $TOP  $CRD >>edit_leap.src 
echo 'savepdb mol ' $PDB >>edit_leap.src 
echo 'quit' >>edit_leap.src 

fi

echo "Adding again solvent box, but with $NUM_NA sodiums and $NUM_CL chlorides" 
rm -f leap.log
$AMBERHOME/bin/tleap -f edit_leap.src > mlog 
mv leap.log  edit_leap.log
tail -1 edit_leap.log

echo 'Randomizing ion positions'
NRES=$(grep 'ATOM  ' $PDB_SOLUTE | tail  -1  | awk '{print $5}')
$AMBERHOME/bin/cpptraj $TOP <<EOF   > mlog 
trajin  $CRD 
autoimage
randomizeions @Na+,Cl- around :1-${NRES} by 6.0  overlap 4.0 
trajout ${CRD}_rand.crd restrt 
go
EOF
echo "$CRD contains now the randomized ion positions"
mv ${CRD}_rand.crd  $CRD

fi 

NATOM_SOLUTE=$(grep -c 'ATOM  ' $PDB_SOLUTE ) 
NATOM=$(grep -c 'ATOM  ' $PDB ) 
echo "Total solute atoms = $NATOM_SOLUTE" 
echo "Total number of atoms in final system= $NATOM" 

rm -f mlog mdinfo
