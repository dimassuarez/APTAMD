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
if [ ! -z "$IONIC" ] && [ -z "$NA_IONIC" ] && [ -z "$MG_IONIC" ]
then
    echo "IONIC=$IONIC is interpretad as NA_IONIC=$IONIC MG_IONIC=0"
    NA_IONIC=$IONIC
    MG_IONIC=0
fi
if [ -z "$NA_IONIC" ]; then NA_IONIC="0.150"; echo 'Assuming NA_IONIC STRENGTH=0.150 M' ; fi
if [ -z "$MG_IONIC" ]; then MG_IONIC="0.000"; echo 'Assuming MG_IONIC STRENGTH=0.000 M' ; fi

if [ -z "$BUFFER_SOLV" ]; then BUFFER_SOLV=16.0; echo 'Assuming BUFFER SOLV= 16.0 A ' ; fi
if [ -z "$FF" ]; then FF=19; echo 'Assuming FF19SB / OPC ' ; fi
if [ -z "$DNAFF" ]; then DNAFF="bsc1"; echo "Assuming DNAFF=$DNAFF" ; fi

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
if [ -z $BOX ]
then
   BOX="OCT"
fi
if [ ${BOX} == "OCT" ]
then
   SOLVATE="solvateOct"
   echo "Truncated Octahedral BOX will be used"
elif [ ${BOX} == "CUBOID" ]
then
   SOLVATE="solvateBox"
   echo "Cuboid BOX will be used"
else
   SOLVATE="solvateOct"
   echo "Truncated Octahedral BOX will be used"
fi

# Checking FF 
if [ "$FF" -ne 19 ] &&  [ "$FF" -ne 14 ] ; then echo "FF=14 or 19, but FF=$FF"; exit ; fi
# Checking DNA FF 
if [ $DNAFF != "bsc1" ] && [  $DNAFF != "OL24" ]  && [  $DNAFF != "OL21" ]  && [  $DNAFF != "OL15" ]
then
        echo "DNA=$DNAFF, but only bsc1 or OL15/OL21/OL24 can be selected" 
        exit
fi
# FF selects WATMODEL and IONFF
if [ $FF -eq 14 ]
then
   WATMODEL="tip3p"
   WATBOX="TIP3PBOX"
else
   WATMODEL="opc"
   WATBOX="OPCBOX"
fi
# Select ion parameters
if [ $WATMODEL == "tip3p" ]
then
        if [ "${MG_IONIC}" == "0.000" ] || [ "${MG_IONIC}" == "0.00" ] || [ "${MG_IONIC}" == "0.0" ] || [ "${MG_IONIC}" == "0" ]
        then
           echo "Using Joung–Cheatham ion parameters for TIP3P"
           IONFF_NA="ionsjc_tip3p"
           IONFF_MG="none"
        else
           echo "Using Li-Merz ion parameters for TIP3P"
           IONFF_NA="ions1lm_126_tip3p"
           IONFF_MG="ions234lm_1264_tip3p"
        fi
else
        if [ "${MG_IONIC}" == "0.000" ] || [ "${MG_IONIC}" == "0.00" ] || [ "${MG_IONIC}" == "0.0" ] || [ "${MG_IONIC}" == "0" ]
        then
                IONFF_NA="ionslm_126_opc"
                IONFF_MG="none"
        else
                IONFF_NA="ionslm_126_opc"
                IONFF_MG="ionslm_1264_opc"
        fi
        echo "Using Li-Merz ion parameters for OPC"
fi

# ALIGN shuould be deactivated for SLAB systems 
if [ -z $ALIGN ]
then
      ALIGN="YES"
fi
if [ ! -e $INITIAL ]; then echo "$INITIAL does not exist in the current location"; exit; fi

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

echo '# Force Field data' >edit_leap_solute.src 
echo "source leaprc.DNA.${DNAFF}"  >>edit_leap_solute.src
echo "source leaprc.protein.ff${FF}SB" >> edit_leap_solute.src
echo "source leaprc.${GAFF}" >> edit_leap_solute.src
echo '# GLYCAM FF '>>edit_leap_solute.src
echo 'source leaprc.GLYCAM_06j-1'>>edit_leap_solute.src
echo 'loadOff GLYCAM_amino_06j_12SB.lib'>>edit_leap_solute.src
echo 'loadOff GLYCAM_aminont_06j_12SB.lib'>>edit_leap_solute.src
echo 'loadOff GLYCAM_aminoct_06j_12SB.lib'>>edit_leap_solute.src
if [ $OFFLIB == "YES" ]; then echo "loadoff ${MOL}.off" >> edit_leap_solute.src ; fi
if [ $FRCMOD == "YES" ]; then echo "loadAmberParams ${MOL}.frcmod" >> edit_leap_solute.src ; fi
if [ $WATMODEL == "opc" ]; then echo 'WAT=OPC' >> edit_leap_solute.src ; fi 
echo "source leaprc.water.${WATMODEL}" >> edit_leap_solute.src
echo "loadamberparams frcmod.${IONFF_NA}" >>edit_leap_solute.src
if [ $IONFF_NA != $IONFF_MG ]; then echo "loadamberparams frcmod.${IONFF_MG}" >>edit_leap_solute.src ; fi
echo '# Build System' >>edit_leap_solute.src 
echo 'mol=loadpdb' $INITMOL >>edit_leap_solute.src 
if [ $EXTRA == "YES" ]; then cat ${MOL}_extra.src  >> edit_leap_solute.src ; fi
echo 'check mol ' >>edit_leap_solute.src 
echo 'charge mol ' >>edit_leap_solute.src 
echo 'saveamberparm mol ' $TOP_SOLUTE  $CRD_SOLUTE >>edit_leap_solute.src 
echo 'savepdb mol ' $PDB_SOLUTE >>edit_leap_solute.src 
echo 'quit' >>edit_leap_solute.src 

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

echo '# Force Field data' >edit_leap.src
echo "source leaprc.DNA.${DNAFF}"  >>edit_leap.src
echo "source leaprc.protein.ff${FF}SB" >> edit_leap.src
echo "source leaprc.${GAFF}" >> edit_leap.src
echo '# GLYCAM FF '>>edit_leap.src
echo 'source leaprc.GLYCAM_06j-1'>>edit_leap.src
echo 'loadOff GLYCAM_amino_06j_12SB.lib'>>edit_leap.src
echo 'loadOff GLYCAM_aminont_06j_12SB.lib'>>edit_leap.src
echo 'loadOff GLYCAM_aminoct_06j_12SB.lib'>>edit_leap.src
if [ $WATMODEL == "opc" ]; then echo 'WAT=OPC' >> edit_leap.src ; fi
echo "source leaprc.water.${WATMODEL}" >> edit_leap.src
echo "loadamberparams frcmod.${IONFF_NA}" >>edit_leap.src
if [ $IONFF_NA != $IONFF_MG ]; then echo "loadamberparams frcmod.${IONFF_MG}" >>edit_leap.src ; fi
if [ $OFFLIB == "YES" ]; then echo "loadoff ${MOL}.off" >> edit_leap.src ; fi
if [ $FRCMOD == "YES" ]; then echo "loadAmberParams ${MOL}.frcmod" >> edit_leap.src ; fi
echo '# Build system' >>edit_leap.src 
echo 'mol=loadpdb' $PDB_SOLUTE >>edit_leap.src 
if [ $EXTRA == "YES" ]; then cat ${MOL}_extra.src  >> edit_leap.src ; fi
if [ ${ALIGN} == "YES" ]; then echo 'alignaxes mol ' >>edit_leap.src; fi 
echo "${SOLVATE} mol $WATBOX $BUFFER_SOLV" >>edit_leap.src
echo 'saveamberparm mol ' $TOP  $CRD >>edit_leap.src 
echo 'savepdb mol ' $PDB >>edit_leap.src 
echo 'quit' >>edit_leap.src 


echo "Adding solvent box with BUFFER SOLV= $BUFFER_SOLV" 
rm -f leap.log
$AMBERHOME/bin/tleap -f edit_leap.src > mlog 
mv leap.log  edit_leap.log
tail -1 edit_leap.log
rm -f mlog

NWAT=$(grep -c 'O   WAT' $PDB)
Q=$(grep 'Total unperturbed charge:' edit_leap_solute.log | awk '{print $4}')

echo "System contains $NWAT waters. Solute charge = $Q. Ionic conc NA = $NA_IONIC MG $MG_IONIC (M) "
$OCTAVE -q  <<EOF  > mlog
num_Na = round( ( ${NA_IONIC} / 55 )* ${NWAT});
num_Mg = round( ( ${MG_IONIC} / 55 )* ${NWAT});
num_Cl = num_Na + 2 * num_Mg  - abs(${Q}) ;
disp([ 'Na= ',num2str(num_Na)])
disp([ 'Cl= ',num2str(num_Cl)])
disp([ 'Mg= ',num2str(num_Mg)])
EOF

NUM_NA=$(head -1 mlog  | awk '{print $2}')
NUM_CL=$(head -2 mlog  | tail -1  | awk '{print $2}')
NUM_MG=$(tail -1 mlog  | awk '{print $2}')

if [ $NUM_CL -lt 0 ] || [ $NUM_NA -lt 0 ]  || [ $NUM_MG -lt 0 ] 
then
        echo "Solute has $Q charge"
        echo "NUM_NA=$NUM_NA  NUM_CL=$NUM_CL   NUM_MG=$NUM_MG"
        echo "Solvent box with BUFFER_SOLV=${BUFFER_SOLV} is too small!!"
        echo "Increase BUFFER_SOLV and run again do_system_edition"
        rm -f *.top *.crd
        exit
fi

if [ $NUM_NA -gt 0 ]  ||  [ $NUM_CL -gt  0 ] || [ $NUM_MG -gt 0 ]    # Only if counterions are needed !
then 

echo '# Force Field data' >edit_leap.src
echo "source leaprc.DNA.${DNAFF}"  >>edit_leap.src
echo "source leaprc.protein.ff${FF}SB" >> edit_leap.src
echo "source leaprc.${GAFF}" >> edit_leap.src
echo '# GLYCAM FF '>>edit_leap.src
echo 'source leaprc.GLYCAM_06j-1'>>edit_leap.src
echo 'loadOff GLYCAM_amino_06j_12SB.lib'>>edit_leap.src
echo 'loadOff GLYCAM_aminont_06j_12SB.lib'>>edit_leap.src
echo 'loadOff GLYCAM_aminoct_06j_12SB.lib'>>edit_leap.src
if [ $WATMODEL == "opc" ]; then echo 'WAT=OPC' >> edit_leap.src ; fi
echo "source leaprc.water.${WATMODEL}" >> edit_leap.src
echo "loadamberparams frcmod.${IONFF_NA}" >>edit_leap.src
if [ $IONFF_NA != $IONFF_MG ]; then echo "loadamberparams frcmod.${IONFF_MG}" >>edit_leap.src ; fi
if [ $OFFLIB == "YES" ]; then echo "loadoff ${MOL}.off" >> edit_leap.src ; fi
if [ $FRCMOD == "YES" ]; then echo "loadAmberParams ${MOL}.frcmod" >> edit_leap.src ; fi
echo '# Build system' >>edit_leap.src 
echo 'mol=loadpdb' $PDB_SOLUTE >>edit_leap.src 
if [ $EXTRA == "YES" ]; then cat ${MOL}_extra.src  >> edit_leap.src ; fi
if [ ${ALIGN} == "YES" ]; then echo 'alignaxes mol ' >>edit_leap.src; fi 
echo "${SOLVATE} mol $WATBOX  $BUFFER_SOLV" >>edit_leap.src
if [ $NUM_NA -gt 0 ]  &&  [ $NUM_CL -gt  0 ]
then
   echo "addionsrand mol  Na+ $NUM_NA  Cl- $NUM_CL " >>edit_leap.src
elif [ $NUM_NA -gt 0 ]
then
   echo "addionsrand mol  Na+ $NUM_NA  " >>edit_leap.src
else
   echo "addionsrand mol  Cl- $NUM_CL " >>edit_leap.src
fi
if [ $NUM_MG -gt 0 ]
then
   echo "addions mol  MG $NUM_MG  " >>edit_leap.src
fi
echo 'saveamberparm mol ' $TOP  $CRD >>edit_leap.src 
echo 'quit' >>edit_leap.src 

echo "Adding again solvent box, but with $NUM_NA sodiums, $NUM_CL chlorides and $NUM_MG magnesiums" 
rm -f leap.log
$AMBERHOME/bin/tleap -f edit_leap.src > mlog 
mv leap.log  edit_leap.log
tail -1 edit_leap.log

echo 'Randomizing ion positions'
NRES=$(grep 'ATOM  ' $PDB_SOLUTE | grep -v 'WAT\|HOH\|Cl-\|Na+|MG\|Mg' | tail  -1  | awk '{print $5}')
$AMBERHOME/bin/cpptraj $TOP <<EOF   > mlog 
trajin  $CRD 
autoimage
randomizeions @Na+,Cl-,MG around :1-${NRES} by 6.0  overlap 4.0 
trajout ${CRD}_rand.crd restrt 
go
EOF
echo "$CRD contains now the randomized ion positions"
mv ${CRD}_rand.crd  $CRD

fi 

# Maybe some crystallographic waters or counterions were added to the 
# initial PDB file. In such a case, it is necessary to filter the _solute. files

NSOLV_IN_SOLUTE=$( sed '1,/FLAG RESIDUE_LABEL/d' $TOP_SOLUTE | sed '/FLAG/,$d'  |  grep 'WAT\|HOH\|Cl-\|Na+\|MG\|Mg' | tr " \t" "\n" | grep -c 'WAT\|HOH\|Na+\|Cl-\|MG\|Mg' )

if [ ${NSOLV_IN_SOLUTE} -gt 0 ]
then

echo "WAT/NA+/Cl-/Mg2+' detected in initial(solute) structure" 
echo "Removing them from $TOP_SOLUTE, $CRD_SOLUTE and $PDB_SOLUTE"
echo "WARNING: Counterion handling may be not adequate for your system!"

$AMBERHOME/bin/parmed -n $TOP_SOLUTE <<EOF
strip :WAT,Na+,Cl-,MG
parmout tmp.top 
go
EOF
mv tmp.top $TOP_SOLUTE
grep -v 'WAT\|HOH\|Na+\|Cl-|MG' $PDB_SOLUTE > tmp.pdb
mv tmp.pdb $PDB_SOLUTE 
$TOOLS/pdbcrd < $PDB_SOLUTE > $CRD_SOLUTE

fi

# For large systems, tLeap PDB file is not readable by Rasmol or other programs
$AMBERHOME/bin/cpptraj $TOP<<EOF  > mlog
trajin $CRD
autoimage
trajout ${PDB} dumpq include_ep
go
EOF

# Adding CCOEF to topology files for 12-6-4 LJ ion potentials
if [ ${IONFF_MG} != "none" ]
then

for topology in ${TOP} ${TOP_SOLUTE}
do

rm -f temp_C4.top
echo "Adding C4 parameters to ${topology} using parmed"
$AMBERHOME/bin/parmed -n  ${topology} <<EOF
setOverwrite True
add12_6_4 :MG
outparm temp_C4.top
EOF
mv -f temp_C4.top ${topology}
done

fi

# Size of the system.....
NATOM_SOLUTE=$(grep -c 'ATOM  ' $PDB_SOLUTE ) 
NATOM=$(grep -c 'ATOM  ' $PDB ) 
echo "Total solute atoms = $NATOM_SOLUTE" 
echo "Total number of atoms in final system= $NATOM" 

rm -f mlog mdinfo
