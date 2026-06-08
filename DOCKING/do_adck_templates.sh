#!/bin/bash

# This scripts prepares the Autodock template files
# that should be inspected before running Docking calcs

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
if [ "$#" -eq 0 ]; then more $APTAMD/DOC/do_adck_templates.txt ; exit; fi

source $APTAMD/ENV/aptamd_env.sh

# Specific MGLTOOLS to be used
AUTODOCKTOOLS="$MGLTOOLS/MGLToolsPckgs/AutoDockTools"
PYTHONSH="$MGLTOOLS/bin/pythonsh"
if [ ! -e $PYTHONSH ] || [ ! -e $AUTODOCKTOOLS ]
then
   echo "$PYTHONSH and/or $AUTODOCKTOOLS do no exist, but required"
   exit
fi

# Checking options
extension="${1##*.}"
if [ $extension == "src" ]
then
    echo "Sourcing $1"
    source $1
else
 # LIST Files containing the PDB files to be processed 
 RECEPTOR_PDB_LIST=$1 
 LIGAND_PDB_LIST=$2
 if [  -z "$RECEPTOR_PDB_LIST" ]; then echo 'Usage: do_adck_templates.sh RECEPTOR_PDB_LIST LIGAND_PDB_LIST' ; exit; fi
 if [  -z "$LIGAND_PDB_LIST" ];   then echo 'Usage: do_adck_templates.sh RECEPTOR_PDB_LIST LIGAND_PDB_LIST'; exit; fi
fi

if [ -z "$NPROCS" ]
then
   NPROCS=$(cat /proc/cpuinfo | grep -c processor)
   if [ $NPROCS -gt 12 ]; then NPROCS=8; fi 
   echo "Using NPROCS=$NPROCS available"
else
   echo "Using NPROCS=$NPROCS as predefined "
fi

# AUTODOCK OPTIONS
if [ -z "$GRID" ]
then
   GRID="0.333"
   echo "Using GRID=${GRID}"
else
   echo "Using GRID=${GRID} as predefined "
fi
if [ -z "$GA_RUN" ]
then
   GA_RUN="50"
   echo "Using GA_RUN=${GA_RUN}"
else
   echo "Using GA_RUN=${GA_RUN} as predefined "
fi
if [ -z "$LIGFORMAT" ]
then
   LIGFORMAT="RIGID"
else
   if [ $LIGFORMAT != "RIGID" ] && [ $LIGFORMAT != "FLEX" ] && [ $LIGFORMAT != "ALLATOM" ] && [ $LIGFORMAT != "TORS" ]
   then
      echo "LIGFORMAT=$LIGFORMAT NOT VALID, SELECTING LIGFORMAT=RIGID"
      LIGFORMAT="RIGID"
   else
      echo "Using LIGFORMAT=$LIGFORMAT as predefined"
   fi 
fi
if [ -z "$LIGCHARGE" ]
then
   LIGCHARGE="GASTEIGER"
else
   if [ $LIGCHARGE != "GASTEIGER" ] && [ $LIGCHARGE != "RESP" ] 
   then
      echo "LIGCHARGE=$LIGCHARG NOT VALID, SELECTING LIGCHARGE=GASTEIGER"
      LIGCHARGE="GASTEIGER"
   else
      echo "Using LIGCHARGE=$LIGCHARGE as predefined"
   fi
fi
if [ -z "$RECCHARGE" ]
then
   RECCHARGE="GASTEIGER"
else
   if [ $RECCHARGE != "GASTEIGER" ] && [ $RECCHARGE != "RESP" ] 
   then
      echo "RECCHARG=$RECCHARG NOT VALID, SELECTING RECCHARG=GASTEIGER"
      RECCHARG="GASTEIGER"
   else
      echo "Using RECCHARGE=$RECCHARGE as predefined"
   fi
fi
if [ -z $KEEPSOLV_REC ]
then
    KEEPSOLV_REC="NO"
elif [ $KEEPSOLV_REC == "YES" ]
then 
    echo 'KEEPSOLV_REC=YES'
    echo 'Na+ ions and H2O molecules preserved in receptor coordinates'
fi 
if [ $KEEPSOLV_REC == "NO" ] ; then echo 'Counterions and H2O molecules are removed from receptor coordinates' ; fi

if [ -z $IATCEN ]
then 
   IATCEN=-1
else
  echo "Selecting ATOM ${IATCEN} in receptor PDB files as grid center"
  echo "Grid size is then determined by ligand size" 
fi 

ZPOT=0
ZAMPL=10.0  # kcal/mol
ZWIDTH=3.0  # Angstrom ( taken from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4815316/ )
if [ -z "${ZATM_REC}" ]
then 
   ZATM_REC=-1
else
  echo "Selecting Z-ATOMS ${ZATM_REC} in receptor PDB for covalent-like constraint"
  let "ZPOT=$ZPOT+1"
fi 
if [ -z "${ZATM_LIG}" ]
then 
   ZATM_LIG=-1
else
  echo "Selecting Z-ATOMS ${ZATM_LIG} in ligand PDB for covalent-like constraint"
  let "ZPOT=$ZPOT+1"
fi 

if [ -z "$RELAX" ]
then
   RELAX="YES"
fi
if [ $RELAX == "YES" ]; then echo "Relaxing structures before ADCK atom typing"; fi 
# Force filed options 

if [ -z "$FF" ]; then FF=19; echo 'Assuming FF19SB  ' ; fi
if [ "$FF" -ne 19 ] &&  [ "$FF" -ne 14 ] ; then echo "FF=14 or 19, but FF=$FF"; exit ; fi

if [ ! -z "$TOPOLOGY_REC" ]
then
   echo "Using TOPOLOGY_REC=${TOPOLOGY_REC}"
   if [ ${TOPOLOGY_REC} != "receptor.top" ]; then rm -f receptor.top;  cp ${TOPOLOGY_REC}  receptor.top ; fi
fi
if [ ! -z "$TOPOLOGY_LIG" ]
then
   echo "Using TOPOLOGY_LIG=${TOPOLOGY_LIG}"
   if [ ${TOPOLOGY_LIG} != "ligand.top" ]; then rm -f ligand.top; cp  ${TOPOLOGY_LIG}  ligand.top ; fi
fi
if [ ! -z "$TOPOLOGY_CMPLX" ]
then
   echo "Using TOPOLOGY_CMPLX=${TOPOLOGY_CMPLX}"
   if [ ${TOPOLOGY_CMPLX} != "complex.top" ]; then rm -f complex.top;  cp ${TOPOLOGY_CMPLX}  complex.top ; fi
fi

WORKDIR=$PWD
echo "WORKDIR=$PWD"
echo "No other working directory is made!"

# Structures with larger size are taken as references
rm -r -f receptor.rmax
for file in $(cat  $RECEPTOR_PDB_LIST)
do
   if [ $KEEPSOLV_REC == "NO" ]
   then
       grep -v 'WAT\|HOH' $file | grep -v 'Na+' | grep -v 'K+' | grep -v 'Cl-' > temp.pdb
   else
       grep -v 'K+\|Cl-\|EPW' $file  > temp.pdb
   fi
   $TOOLS/center_geom -noatm -i temp.pdb -o temp.size  > mlog
   R=$(grep 'COMMENT  RADMAX' temp.size  |  awk '{printf("%f\n",$4)}')
   echo $file $R >> receptor.rmax
   rm -f temp.pdb
done
R_LARGE=$(sort -n -k 2 receptor.rmax | tail -1 | awk '{print $1}')
cp $R_LARGE receptor.pdb
if [ $KEEPSOLV_REC == "NO" ]
then
       grep -v 'WAT\|HOH' receptor.pdb | grep -v 'Na+' | grep -v 'K+' | grep -v 'Cl-' > temp.pdb
else
       grep -v 'K+\|Cl-\|EPW' receptor.pdb  > temp.pdb
fi
mv -f temp.pdb receptor.pdb
rm -f temp.size
   
rm -r -f ligand.rmax
for file in $(cat  $LIGAND_PDB_LIST)
do
   grep -v 'WAT\|HOH' $file  | grep -v 'Na+' | grep -v 'K+' | grep -v 'Cl-' > temp.pdb
   $TOOLS/center_geom -noatm -i temp.pdb  -o temp.size  > mlog
   R=$(grep 'COMMENT  RADMAX' temp.size  |  awk '{printf("%f\n",$4)}')
   echo $file $R >> ligand.rmax
   rm -f temp.pdb 
done
L_LARGE=$(sort -n -k 2 ligand.rmax | tail -1 | awk '{print $1}')
cp $L_LARGE ligand.pdb
rm -f temp.size

echo "Selected receptor = $R_LARGE "
echo "Selected ligand = $L_LARGE "

# Recentering coordinates
cp receptor.pdb temp.pdb
$TOOLS/center_geom -i temp.pdb -o receptor.pdb
grep COMMENT receptor.pdb > receptor.size
sed -i '/^COMMENT/d' receptor.pdb
rm -f temp.pdb

cp ligand.pdb  temp.pdb
$TOOLS/center_geom -i temp.pdb -o ligand.pdb
grep COMMENT ligand.pdb > ligand.size 
sed -i '/^COMMENT/d'  ligand.pdb 
rm -f temp.pdb

# Determining Receptor and Ligand dimensions 
DIM_REC=$(grep 'XYZ DIM' receptor.size |  awk '{printf("%f %f %f\n",$5,$6,$7)}')
R_REC=$(grep 'COMMENT  RADMAX' receptor.size |  awk '{printf("%f\n",$4)}')
DIM_LIG=$(grep 'XYZ DIM' ligand.size |  awk '{printf("%f %f %f\n",$5,$6,$7)}')
R_LIG=$(grep 'COMMENT  RADMAX' ligand.size  |  awk '{printf("%f\n",$4)}')
echo "DIM_REC=$DIM_REC"
echo "R_REC=$R_REC"
echo "DIM_LIG=$DIM_LIG"
echo "R_LIG=$R_LIG"


# It is important to determine wheter the ligand is a peptide or an aptamer
LIG_TYPE='UNK'
for molecule in ligand receptor
do

  grep 'ATOM ' ${molecule}.pdb | awk '{printf(" %s %3i\n",$4,$5)}' | uniq > tmp.res
  NRES=$(grep 'ATOM ' ${molecule}.pdb | grep -v 'Na+' | grep -v 'WAT'  | tail -1 | awk '{print $5}')

  nd=0
  for res in DA DT DC DG 
  do 
     n=$(grep -c $res tmp.res)
     let "nd=$nd+$n"
  done
  naa=0 
  for res in ACE ALA ARG ASH ASN ASP CYM CYS CYX GLH GLN GLU GLY \
  HID HIE HIP HYP ILE LEU LYN LYS MET PHE PRO SER THR TRP TYR VAL NME 
  do 
     n=$(grep -c $res tmp.res)
     let "naa=$naa+$n"
  done
  rm -f tmp.res
  if [ $NRES -eq $nd ]
  then 
     echo "Good $molecule seems to be a nice aptamer"
     MOL_TYPE="APT"
     if [ ${molecule} == "ligand" ]; then LIG_TYPE="APT"; fi
  elif [ $NRES -eq $naa ]
  then
     echo "Good $molecule seems to be a nice peptide"
     MOL_TYPE="PEP"
     if [ ${molecule} == "ligand" ]; then LIG_TYPE="PEP"; fi
  else
     echo "Sorry, no luck in guessing $molecule type " 
     if [ ! -e ${molecule}.top ]
     then 
        echo "${molecule}.top does not exist nor is defined"
        echo "Provide TOPOLOGY_REC or TOPOLOGY_LIG option "
        echo "and rerun do_adck_templates"
        exit
     fi
     MOL_TYPE='UNKNOWN'
  fi
  ISOLV=$(grep 'ATOM  ' ${molecule}.pdb | grep -c 'Na+\|WAT')
  if  [ ! -e ${molecule}.top ]
  then 
    echo "Creating ${molecule}.top"
    echo '# Force Field data' >  edition_${molecule}.src
    echo 'source leaprc.DNA.bsc1'  >> edition_${molecule}.src
    echo "source leaprc.protein.ff${FF}SB" >> edition_${molecule}.src
    if  [ $ISOLV -gt 0 ] && [ "$FF" -eq 19 ]
    then
       echo 'WAT=OPC' >> edition_${molecule}.src
       echo 'source leaprc.water.opc' >> edition_${molecule}.src
       echo 'loadamberparams frcmod.ionslm_iod_opc' >>edition_${molecule}.src 
    elif [ $ISOLV -gt 0 ]
    then 
       echo 'source leaprc.water.tip3p' >>edition_${molecule}.src
       echo 'loadamberparams frcmod.ionsjc_tip3p'>>edition_${molecule}.src
    fi
    echo "mol=loadpdb  ${molecule}.pdb "  >> edition_${molecule}.src
    echo 'check mol ' >>   edition_${molecule}.src
    echo 'charge mol ' >>   edition_${molecule}.src
    echo "saveamberparm mol  ${molecule}.top dummy.crd "  >>  edition_${molecule}.src
    echo "saveoff mol   ${molecule}.off   "  >>  edition_${molecule}.src
    echo 'quit' >>   edition_${molecule}.src
    $AMBERHOME/bin/tleap -f  edition_${molecule}.src > edition_${molecule}.log 
    rm -f leap.log 
    if [ $molecule == "receptor" ] ; then TOPOLOGY_REC="receptor.top"; fi
    if [ $molecule == "ligand" ] ; then TOPOLOGY_LIG="ligand.top"; fi
  fi
done
if [ -z $TOPOLOGY_REC ]; then TOPOLOGY_REC="receptor.top"; fi 
if [ -z $TOPOLOGY_LIG ]; then TOPOLOGY_LIG="ligand.top"; fi 

# Relaxing structures prior to atom typing with Autodock
if [ "$RELAX" == "YES" ]
then 

# Minimization of the nucleotide bases PO4 backbone fixed 
echo "Minimal relaxation of RECEPTOR and LIGAND geometries for ADCK atom typing"

cat  <<EOF  >  sander_min.inp 
Partial relaxation  Distance dependent eps
&cntrl
 imin=1, ncyc=100,  maxcyc=100, ntmin=2, drms=0.02
 ntb=0, ntf=1, ntc=1, ntpr=500,
 cut=1000.0, nsnb=5000,  igb=0,
 ibelly=1,
 bellymask="DUMMY_BELLY"
 /
 &ewald
  eedmeth=5,
 /
EOF

for molecule in ligand receptor 
do
   $AMBERHOME/bin/cpptraj ${molecule}.top <<EOF > mlog
trajin ${molecule}.pdb
trajout ${molecule}.crd restart
go
EOF
   NAT=$(grep -c 'ATOM  ' ${molecule}.pdb)
   NRES=$(grep 'ATOM ' ${molecule}.pdb | grep -v 'Na+' | grep -v 'WAT'  | tail -1 | awk '{print $5}')
   sed "s/DUMMY_BELLY/:1-${NRES}/" sander_min.inp  > sander_min_${molecule}.inp 
   echo " Relaxing ${molecule}"
   if [ $NAT -lt 500 ]
   then 
     $AMBERHOME/bin/sander -O -i sander_min_${molecule}.inp -p ${molecule}.top -c ${molecule}.crd -r ${molecule}.rst -o ${molecule}.out
   else
     $MPI_HOME/bin/mpirun -q -np $NPROCS $AMBERHOME/bin/sander.MPI -O -i sander_min_${molecule}.inp -p ${molecule}.top -c ${molecule}.crd -r ${molecule}.rst -o ${molecule}.out
   fi
   $AMBERHOME/bin/cpptraj ${molecule}.top <<EOF > mlog
trajin  ${molecule}.rst 
trajout ${molecule}.pdb pdb pdbatom 
go
EOF
# Preparing a pqr Amber file 
  $AMBERHOME/bin/cpptraj ${molecule}.top  <<EOF > mlog
  trajin  ${molecule}.rst 
  trajout temp.pdb dumpq 
  go
EOF
  mv temp.pdb ${molecule}.pqr

done

fi

# Build up the PDBQT file for the RECEPTOR.

echo "prepare_receptor4.py -r receptor.pdb"

# Convert AMBER atom/residue names to std names before using
# prepare_receptor4.py / prepare_ligand4.py 
# But we cannot forget AMBER residue names for ligand reconstruction!

for molecule in receptor ligand
do
	cp ${molecule}.pdb ${molecule}_amber.pdb
	$TOOLS/std_atmresname.sh -atmres ${molecule}.pdb
done

if [ $KEEPSOLV_REC == "YES" ]
then

grep 'Na+\|WAT' receptor.pdb > receptor_ionsolv.pdb
grep -v 'Na+\|WAT' receptor.pdb > receptor_solute.pdb 
$PYTHONSH $AUTODOCKTOOLS/Utilities24/prepare_receptor4.py -r receptor_solute.pdb
cat receptor_solute.pdbqt receptor_ionsolv.pdb > receptor.pdbqt
# fixq_pdbqt assigns Na+/OW, H1,H2 types and charges
$TOOLS/fixq_pdbqt receptor.pdbqt tmp_receptor.pdbqt receptor_amber.pdb ; mv -f tmp_receptor.pdbqt receptor.pdbqt
$TOOLS/reorder < receptor.pdbqt > tmp_receptor.pdbqt ; mv -f tmp_receptor.pdbqt receptor.pdbqt
   
else 

$PYTHONSH $AUTODOCKTOOLS/Utilities24/prepare_receptor4.py -r receptor.pdb
$TOOLS/fixq_pdbqt receptor.pdbqt tmp_receptor.pdbqt receptor_amber.pdb ; mv -f tmp_receptor.pdbqt receptor.pdbqt

fi

#Fix atom names in pdbqt
$TOOLS/std_atmresname.sh -atm receptor.pdbqt
# Replace the atom type Ho (in terminal residues of aptamers) by HD
sed -i 's/Ho/HD/' receptor.pdbqt

# Build up the PDBQT file for the LIGAND. Option -Z means that all torsions are inactive.
# -U lps means that all atoms are preserved 
echo "prepare_ligand4.py ligand.pdb "

# Check if ligand has more than one TER record meaning different molecules 
NTER=$(grep -c TER ligand.pdb)
if [ $NTER -gt  1 ]
then
      echo "ligand.pdb has more than one TER record"
      echo "prepare_ligand4.py will be applied sequentially on every fragment."
      if [ $LIGFORMAT == "FLEX" ] || [ $LIGFORMAT == "TORS" ]
      then
            echo "Flexible docking with multimolecule ligand is NOT valid."
            exit
      fi
      csplit -n 3 -s -f ligand_frag_  ligand.pdb  '/TER/' '{*}'
      ifrag=0
      echo 'ROOT' > ligand.pdbqt
      rm -f ligand.frag
      for lfrag in $(ls ligand_frag_*)
      do
          NATFRAG=$(grep -c 'ATOM  ' $lfrag )
          if [ $NATFRAG -eq 0 ]; then continue ; fi
          let "ifrag=$ifrag+1"
          rm -f ligand_frag.pdbqt 
          grep 'ATOM  ' $lfrag > ligand_temp.pdb 
	  if [ $LIGFORMAT == "ALLATOM" ]
          then
             echo "Selecting ligand_allatom_notors.pdbqt for fragment $ifrag"
             $PYTHONSH $AUTODOCKTOOLS/Utilities24/prepare_ligand4.py -U lps -Z -l ligand_temp.pdb -o  ligand_frag.pdbqt
          else
             echo "Selecting ligand_notors.pdbqt for fragment $ifrag"
             $PYTHONSH $AUTODOCKTOOLS/Utilities24/prepare_ligand4.py -Z -l ligand_temp.pdb -o ligand_frag.pdbqt
          fi
          grep 'ATOM  \|HETATM' ligand_frag.pdbqt >> ligand.pdbqt
          NATFRAG=$(grep -c 'ATOM  \|HETATM' ligand_frag.pdbqt) 
          echo $NATFRAG >> ligand.frag
      done
      rm -f ligand_frag* ligand_temp.pdb 
      $TOOLS/fixq_pdbqt ligand.pdbqt tmp_ligand.pdbqt ligand_amber.pdb; mv -f tmp_ligand.pdbqt  ligand.pdbqt
      echo 'ENDROOT' >> ligand.pdbqt
      echo 'TORSDOF   0 ' >> ligand.pdbqt

else

if [ $LIGFORMAT == "ALLATOM" ]
then
     echo "Selecting ligand_allatom_notors.pdbqt"
     $PYTHONSH $AUTODOCKTOOLS/Utilities24/prepare_ligand4.py -U lps -Z -l ligand.pdb -o ligand_allatom_notors.pdbqt
     $TOOLS/fixq_pdbqt ligand_allatom_notors.pdbqt tmp_ligand.pdbqt ligand_amber.pdb; mv -f tmp_ligand.pdbqt  ligand_allatom_notors.pdbqt
     echo 'ROOT' > ligand.pdbqt
     grep 'ATOM ' ligand_allatom_notors.pdbqt >> ligand.pdbqt
     echo 'ENDROOT' >> ligand.pdbqt
     echo 'TORSDOF   0 ' >> ligand.pdbqt
elif [ $LIGFORMAT == "RIGID" ]
then
     echo "Selecting ligand_notors.pdbqt"
     $PYTHONSH $AUTODOCKTOOLS/Utilities24/prepare_ligand4.py -Z -l ligand.pdb -o ligand_notors.pdbqt
     $TOOLS/fixq_pdbqt ligand_notors.pdbqt tmp_ligand.pdbqt ligand_amber.pdb ; mv -f tmp_ligand.pdbqt  ligand_notors.pdbqt
     echo 'ROOT' > ligand.pdbqt
     grep 'ATOM ' ligand_notors.pdbqt >> ligand.pdbqt
     echo 'ENDROOT' >> ligand.pdbqt
     echo 'TORSDOF   0 ' >> ligand.pdbqt
elif [ $LIGFORMAT == "FLEX" ] || [ $LIGFORMAT == "TORS" ]
then
     echo "Selecting ligand_tors.pdbqt"
     $PYTHONSH $AUTODOCKTOOLS/Utilities24/prepare_ligand4.py -l ligand.pdb -o ligand_tors.pdbqt
     $TOOLS/fixq_pdbqt ligand_tors.pdbqt tmp_ligand.pdbqt ligand_amber.pdb ; mv -f tmp_ligand.pdbqt  ligand_tors.pdbqt
     cp ligand_tors.pdbqt ligand.pdbqt
else
     echo "Selecting ligand_notors.pdbqt"
     $PYTHONSH $AUTODOCKTOOLS/Utilities24/prepare_ligand4.py -Z -l ligand.pdb -o ligand_notors.pdbqt
     $TOOLS/fixq_pdbqt ligand_notors.pdbqt tmp_ligand.pdbqt ligand_amber.pdb; mv -f tmp_ligand.pdbqt  ligand_notors.pdbqt
     echo 'ROOT' > ligand.pdbqt
     grep 'ATOM ' ligand_notors.pdbqt >> ligand.pdbqt
     echo 'ENDROOT' >> ligand.pdbqt
     echo 'TORSDOF   0 ' >> ligand.pdbqt
fi
NATFRAG=$(grep -c 'ATOM  \|HETATM' ligand.pdbqt) 
echo $NATFRAG > ligand.frag

fi

#Fix atom names in pdbqt
$TOOLS/std_atmresname.sh -atm ligand.pdbqt
# Replace the atom type Ho (in terminal residues of aptamers) by HD
sed -i 's/Ho/HD/' ligand*pdbqt

# Changing GASTEIGER charges by RESP charges if required
if [ $LIGCHARGE == "RESP" ]
then 
     $APTAMD/AUXTOOLS/addZ_pqr_pdbqt  ligand.pqr ligand.pdbqt temp.pdbqt
     mv temp.pdbqt ligand.pdbqt
fi
if [ $LIGCHARGE == "GASTEIGER" ] && [ $LIGFORMAT == "FLEX" ]
then 
     $APTAMD/AUXTOOLS/addZ_pqr_pdbqt ligand.pqr ligand.pdbqt temp.pdbqt keep_charge
     mv temp.pdbqt ligand.pdbqt
fi
if [ $RECCHARGE == "RESP" ]
then 
     $APTAMD/AUXTOOLS/addZ_pqr_pdbqt  receptor.pqr receptor.pdbqt temp.pdbqt
     mv temp.pdbqt receptor.pdbqt
fi


# Incoporating Z-ATM if defined
if [ $ZPOT -eq  2 ]
then
   rm -f receptor.zatm ligand.zatm
   for IDATM in ${ZATM_REC}
   do
       echo "ZATM    ${IDATM} " >> receptor.zatm 
   done
   for IDATM in ${ZATM_LIG}
   do
       echo "ZATM    ${IDATM} " >> ligand.zatm 
   done
   cat receptor.zatm >> receptor.pqr
   grep -v PDBQT_PQR_MATCH receptor.pdbqt > tmp; mv tmp receptor.pdbqt
   $APTAMD/AUXTOOLS/addZ_pqr_pdbqt  receptor.pqr receptor.pdbqt temp.pdbqt keep_charge
   mv temp.pdbqt receptor.pdbqt
   cat ligand.zatm >> ligand.pqr
   grep -v PDBQT_PQR_MATCH ligand.pdbqt > tmp; mv tmp ligand.pdbqt
   $APTAMD/AUXTOOLS/addZ_pqr_pdbqt  ligand.pqr ligand.pdbqt temp.pdbqt keep_charge
   mv temp.pdbqt ligand.pdbqt
fi

# Determining Grid dimensions (We expand receptor dimensions by ligand max radius)
if [ $IATCEN -le 0 ]
then 
   SIZE=$(echo $DIM_REC $R_LIG $GRID | awk '{ix=1.5*($1+$4)/$5; iy=1.5*($2+$4)/$5; iz=1.5*($3+$4)/$5; printf("%i,%i,%i \n",ix,iy,iz)}')
else 
   SIZE=$(echo $DIM_LIG $GRID | awk '{ix=1.5*$1/$4; iy=1.5*$2/$4; iz=1.5*$3/$4; printf("%i,%i,%i \n",ix,iy,iz)}')
fi
echo "SIZE=$SIZE"
# Preparing input for AUTOGRID 
echo prepare_gpf4.py -l ligand.pdbqt -r receptor.pdbqt  -o receptor_ligand.gpf -p npts=${SIZE} -p spacing="${GRID}" -p gridcenter="0.,0.,0." 
$PYTHONSH $AUTODOCKTOOLS/Utilities24/prepare_gpf4.py -p parameter_file=$APTAMD/DOCKING/AD4_parameters.dat  \
    -l ligand.pdbqt -r receptor.pdbqt  -o receptor.gpf \
    -p npts="${SIZE}" -p spacing="${GRID}" -p gridcenter="0.,0.,0."

# Preparing input for AUTODOCK
echo prepare_dpf4.py -l ligand.pdbqt -r receptor.pdbqt -o receptor_ligand.dpf -p ga_run="${GA_RUN}"
$PYTHONSH $AUTODOCKTOOLS/Utilities24/prepare_dpf4.py -p parameter_file=$APTAMD/DOCKING/AD4_parameters.dat   \
   -l ligand.pdbqt -r receptor.pdbqt -o receptor_ligand.dpf -p ga_run="${GA_RUN}"

# Grid center to be particularized = coordinates of the IATCEN atom of receptor) 
if [ $IATCEN -gt 0 ] 
then
   sed -i '/gridcenter/d' receptor.gpf 
   echo "#GRIDCENTER COORDINATES: TO BE TAKEN FROM ATOM IATCEN IN RECEPTOR PDB FILE"  >> receptor.gpf
   echo "#IATCEN= ${IATCEN}"  >> receptor.gpf
   echo 'gridcenter DUMMY_IATCEN_COORD   # xyz-coordinates ' >> receptor.gpf
else
   echo '#IATCEN= -1' >> receptor.gpf
fi

# Zatm tuning of receptor.gpf
if [ ${ZPOT} -eq  2 ]
then 
   echo "#COVALENTMAP  : ZATM potential " >> receptor.gpf
   echo "#IZATM= 1"  >> receptor.gpf
   echo "covalentmap  ${ZWIDTH}  ${ZAMPL}  DUMMY_ZPOT_CENTER  " >> receptor.gpf
else
   echo '#IZATM= -1' >> receptor.gpf
fi

# Preparing a ligand library  (required for the reconstruction of the full ligand atoms
# after docking  
if [ ! -e ligand.off ]
  then
     echo "Preparing ligand library (OFF file) from ligand.top "
     $AMBERHOME/bin/parmed  ligand.top  <<EOF > mlog
     loadcoordinates ligand.pdb
     writeOFF ligand.off
     go
EOF
fi
# Leap input for ligand reconstruction 
echo '# Force Field data' >  input_leap_ligand.src
echo 'source leaprc.DNA.bsc1'  >> input_leap_ligand.src
echo "source leaprc.protein.ff${FF}SB" >> input_leap_ligand.src
echo "loadoff ${WORKDIR}/ligand.off" >> input_leap_ligand.src
# Make sure using addpdbresmap that terminal residues are correctly mapped
RESINI=$(grep 'ATOM  ' ligand.pdbqt | head -1 | awk '{print $4}')
RESEND=$(grep 'ATOM  ' ligand.pdbqt | tail -1 | awk '{print $4}')
if [ ${LIG_TYPE} == "PEP" ]
then 
    IFIXED=$(echo $RESINI | grep -c  "N...")
    if [ $IFIXED -eq 0  ]; then echo " addpdbresmap { { 0 \"${RESINI}\" \"N${RESINI}\" } } " >> input_leap_ligand.src ;fi
    IFIXED=$(echo $RESEND | grep -c  "C...")
    if [ $IFIXED -eq 0  ]; then echo " addpdbresmap { { 1 \"${RESEND}\" \"C${RESEND}\" } } " >> input_leap_ligand.src ;fi
elif [ ${LIG_TYPE} == "APT" ] 
then
    IFIXED=$(echo $RESINI | grep -c  "D.5")
    if [ $IFIXED -eq 0 ]; then echo " addpdbresmap { { 0 \"${RESINI}\" \"${RESINI}5\" } }   " >> input_leap_ligand.src ; fi 
    IFIXED=$(echo $RESEND | grep -c  "D.3")
    if [ $IFIXED -eq 0 ]; then echo " addpdbresmap { { 1 \"${RESEND}\" \"${RESEND}3\" } }   " >> input_leap_ligand.src ; fi 
fi
echo "mol = loadpdb ligand.pdb" >> input_leap_ligand.src
echo "savepdb mol ligand_edited.pdb" >> input_leap_ligand.src
echo 'quit' >> input_leap_ligand.src


# Preparing complex.top 
if [ ! -e complex.top ]
then 

CMAP_REC=$(grep -i -c 'FLAG CMAP_PARAMETER' receptor.top  ) 
CMAP_LIG=$(grep -i -c 'FLAG CMAP_PARAMETER' ligand.top    ) 

if [ $CMAP_REC -gt 0 ] || [ $CMAP_LIG -gt 0 ]
then
  echo "CMAP parameters present in receptor.top and/or ligand.top"
  echo "Merging of topologies using parmed is not feasible"
fi
  
if [ $CMAP_REC -eq 0 ] && [ $CMAP_LIG -eq 0 ]
then
  echo "Merging ${TOPOLOGY_REC} and ${TOPOLOGY_LIG} into complex.top" 

# If topologies of the receptor and ligand are available
# then we merge them using python and parmed library
   $AMBERHOME/miniconda/bin/python <<EOF
import parmed as pmd
parm1 = pmd.load_file('${TOPOLOGY_REC}','receptor_amber.pdb' )
parm2 = pmd.load_file('${TOPOLOGY_LIG}','ligand_amber.pdb')
joined = parm1 + parm2
joined.save('complex.parm7')
joined.save('complex.rst7')
EOF

  mv complex.parm7  complex.top
  mv complex.rst7   complex.crd

else
# Inserting Tleap  edition details may be necessary here....
  echo "Trying tLEaP edition to build complex.top" 
  echo "...maybe additional edition commands are required!"
  echo "Check the edition.src file."
  echo '# Force Field data' >edition.src
  echo '# DNA FF ' >>edition.src
  echo 'source leaprc.DNA.bsc1' >>edition.src
  echo '#  Protein and Solvent FFs ' >>edition.src
  if [ "$FF" -eq 19 ]
  then
    echo 'source leaprc.protein.ff19SB' >> edition.src
    echo 'WAT=OPC' >> edition.src
    echo 'source leaprc.water.opc' >> edition.src
    echo 'loadamberparams frcmod.ionslm_iod_opc' >>edition.src 
  else 
    echo 'source leaprc.protein.ff14SB'>>edition.src
    echo 'source leaprc.water.tip3p' >>edition.src
    echo 'loadamberparams frcmod.ionsjc_tip3p'>>edition.src
  fi
  echo '# GLYCAM FF '>>edition.src
  echo '# source leaprc.GLYCAM_06j-1'>>edition.src
  echo '#loadOff GLYCAM_amino_06j_12SB.lib'>>edition.src
  echo '#loadOff GLYCAM_aminont_06j_12SB.lib'>>edition.src
  echo '#loadOff GLYCAM_aminoct_06j_12SB.lib'>>edition.src
  echo 'receptor=loadpdb receptor_amber.pdb' >> edition.src
  echo 'ligand=loadpdb ligand_amber.pdb' >>edition.src
  echo 'complex=combine{receptor ligand}'>> edition.src
  echo 'saveamberparm complex complex.top complex.crd' >>edition.src
  echo 'quit'>>edition.src
  $AMBERHOME/bin/tleap -f edition.src > mlog
  tail -1 mlog

  echo "Check carefully the results of the tleap edition for the complex!"
  echo "Adapt edition.src if necessary and rerun tleap"

fi 


fi

rm -f mlog dummy.crd mdinfo 

echo "Check out all template files before launching do_adck.sh"
