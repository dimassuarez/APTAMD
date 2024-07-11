#!/bin/bash 

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi

if [ "$#" -eq 0 ]; then more $APTAMD/DOC/do_aptamer_edition.txt ; exit; fi

# do_aptamer_edition: Starting from aptamer pdb file, this scripts performs the following tasks:
#    1) partial relaxation of sugar and base moieties"
#    2) tleap edition adding counterions to satisfy the required ionic force"

source $APTAMD/ENV/aptamd_env.sh 

extension="${1##*.}"
if [ $extension == "src" ]
then 
    echo "Sourcing $1"
    source $1
else
    INITIAL=$1
fi
if [ -z "$INITIAL" ]; then more $APTAMD/DOC/do_aptamer_edition.txt ; exit; fi
if [ -z "$IONIC" ]; then IONIC="0.150"; echo 'Assuming IONIC STRENGTH=0.150 M' ; fi
if [ -z "$BUFFER_SOLV" ]; then BUFFER_SOLV=16.0; echo 'Assuming BUFFER SOLV= 16.0 A ' ; fi
if [ -z "$MAXCYC" ]; then MAXCYC=500; echo 'Assuming MAXCYC=500 ' ; fi

# 
if [ ! -e $INITIAL ]; then echo "$INITIAL does not exist in the current location"; exit; fi

# Filenames
MOL=${INITIAL/_initial.pdb/}
if [ $MOL == $INITIAL ]; then echo "$INITIAL does not conform MOL_initial.pdb filename"; exit; fi

TOP_SOLUTE=${MOL}_solute.top
CRD_SOLUTE=${MOL}_solute.crd
PDB_SOLUTE=${MOL}_solute.pdb
TOP=${MOL}.top
CRD=${MOL}.crd
PDB=${MOL}.pdb

# Some cleaning of the INITIAL PDB file is required 
# when generated by RNAComposer
sed -i '/MODEL/d'  $INITIAL
sed -i '/ENDMDL/d' $INITIAL
# Remove the Uracyl H5 atoms that it is replaced by a CH3 in Thymine. 
sed -i '/H5    U/d' $INITIAL
# Remove the Ribose O2' and HO2' atoms to get the Deoxiribose sugar of DNA.
sed -i "/O2' /d" $INITIAL
sed -i "/HO2' /d" $INITIAL
# Rename residues Uracyl to Thymine.
sed -i 's/U/T/' $INITIAL
#Rename all residues according to AMBER names
sed -i 's/ T A/DT A/' $INITIAL
sed -i 's/ A A/DA A/' $INITIAL
sed -i 's/ G A/DG A/' $INITIAL
sed -i 's/ C A/DC A/' $INITIAL

# Edition with tleap to build CH3 in Thymine.
echo '# Force Field data' >edit_leap_solute.src 
echo 'source leaprc.DNA.bsc1'  >>edit_leap_solute.src
echo 'loadamberparams frcmod.ionsjc_tip3p' >>edit_leap_solute.src 
echo 'source leaprc.RNA.OL3'  >>edit_leap_solute.src
echo '# Build Aptamer' >>edit_leap_solute.src 
echo 'apt=loadpdb' $INITIAL >>edit_leap_solute.src 
echo 'check apt' >>edit_leap_solute.src 
echo 'charge apt' >>edit_leap_solute.src 
echo 'saveamberparm apt ' $TOP_SOLUTE  $CRD_SOLUTE >>edit_leap_solute.src 
echo 'savepdb apt ' $PDB_SOLUTE >>edit_leap_solute.src 
echo 'quit' >>edit_leap_solute.src 

rm -f leap.log
echo  Editing solute coordinates
$AMBERHOME/bin/tleap -f edit_leap_solute.src > mlog 
mv leap.log  edit_leap_solute.log

# Minimization of the nucleotide bases PO4 backbone fixed 
cat  <<EOF  > sander_relax_solute.inp 
Partial relaxation (PO4 backbone atom fixed)  Distance dependent eps
&cntrl
 imin=1, ncyc=100,  maxcyc=${MAXCYC},  ntmin=2, drms=0.02
 ntb=0, ntf=1, ntc=1, ntpr=500,
 cut=1000.0, nsnb=5000,  igb=0,
 ibelly=1,
 bellymask="!@P,OP1,OP2,O3',O5'"
 /
 &ewald
  eedmeth=5,
 /
EOF

echo "Running sander to relax sugar and base moieties (PO4 backbone fixed)"

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


# Edition with tleap to add solvent box 
# First we edit without adding counter ions
echo '# Force Field data' >edit_leap.src 
echo 'source leaprc.DNA.bsc1'  >>edit_leap.src
echo 'source leaprc.water.tip3p'  >>edit_leap.src
echo 'loadamberparams frcmod.ionsjc_tip3p' >>edit_leap.src 
echo 'source leaprc.RNA.OL3'  >>edit_leap.src
echo '# Build Aptamer' >>edit_leap.src 
echo 'apt=loadpdb' $PDB_SOLUTE >>edit_leap.src 
echo 'alignaxes apt' >>edit_leap.src
echo 'solvateOct apt TIP3PBOX ' $BUFFER_SOLV >>edit_leap.src
echo 'saveamberparm apt ' $TOP  $CRD >>edit_leap.src 
echo 'savepdb apt ' $PDB >>edit_leap.src 
echo 'quit' >>edit_leap.src 

echo "Adding solvent box with BUFFER SOLV= $BUFFER_SOLV" 
rm -f leap.log
$AMBERHOME/bin/tleap -f edit_leap.src > mlog 
mv leap.log  edit_leap.log
rm -f mlog

NWAT=$(grep -c 'O   WAT' $PDB)
Q=$(grep 'Total unperturbed charge:' edit_leap_solute.log | awk '{print $4}')

echo "System contains $NWAT waters. Solute charge = $Q. Ionic strenght= $IONIC (M) " 
$OCTAVE -q  <<EOF  > mlog
num_Na = round( ( ${IONIC} / 55 )* ${NWAT});
num_Cl = num_Na - abs(${Q}) ; 
disp([ 'Na= ',num2str(num_Na)])
disp([ 'Cl= ',num2str(num_Cl)])
EOF

cat mlog 

NUM_NA=$(head -1 mlog  | awk '{print $2}')
NUM_CL=$(tail -1 mlog  | awk '{print $2}')

echo '# Force Field data' >edit_leap.src 
echo 'source leaprc.DNA.bsc1'  >>edit_leap.src
echo 'source leaprc.water.tip3p'  >>edit_leap.src
echo 'loadamberparams frcmod.ionsjc_tip3p' >>edit_leap.src 
echo 'source leaprc.RNA.OL3'  >>edit_leap.src
echo '# Build Aptamer' >>edit_leap.src 
echo 'apt=loadpdb' $PDB_SOLUTE >>edit_leap.src 
echo 'alignaxes apt' >>edit_leap.src
echo 'solvateOct apt TIP3PBOX ' $BUFFER_SOLV >>edit_leap.src
echo "addionsrand apt  Na+ $NUM_NA  Cl- $NUM_CL " >>edit_leap.src
echo 'saveamberparm apt ' $TOP  $CRD >>edit_leap.src 
echo 'savepdb apt ' $PDB >>edit_leap.src 
echo 'quit' >>edit_leap.src 

echo "Adding again solvent box, but with $NUM_NA sodiums and $NUM_CL chlorides" 
rm -f leap.log
$AMBERHOME/bin/tleap -f edit_leap.src > mlog 
mv leap.log  edit_leap.log

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

NATOM_SOLUTE=$(grep -c 'ATOM  ' $PDB_SOLUTE ) 
NATOM=$(grep -c 'ATOM  ' $PDB ) 
echo "Total solute atoms = $NATOM_SOLUTE" 
echo "Total number of atoms in final system= $NATOM" 


rm -f mlog mdinfo
