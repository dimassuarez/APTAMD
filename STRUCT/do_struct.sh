#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
if [ "$#" -eq 0 ]; then more $APTAMD/DOC/do_struct.txt ; exit; fi

source $APTAMD/ENV/aptamd_env.sh
if [ -z "$X3DNAHOME" ]; then echo 'X3DNAHOME is not available, but needed'; exit; fi

# Checking options
extension="${1##*.}"
if [ $extension == "src" ]
then
    echo "Sourcing $1"
    source $1
else
 # Put into double quotes the list of trajectories to be processed
 MD_TRAJ=$1
 #  GAMD or MD trajectory type should be specified 
 MD_TYPE=$2
fi

if [ ! -z "$MOL" ] && [  -z "$MD_TRAJ" ]; then MD_TRAJ="${MOL}"; echo $MD_TRAJ; unset MOL; fi

if [ -z "$MD_TRAJ" ]; then echo 'Usage: do_struct.sh "MOL1   MOL2 .."   [ MD | GAMD ] '; exit ; fi
if [ -z "$MD_TYPE" ]; then echo 'Usage: do_struct.sh "MOL1   MOL2 .."   [ MD | GAMD ] '; exit ; fi
if [ -z "$MD_PROD" ]; then MD_PROD="5.PRODUCTION"; else echo "Assuming MD_PROD=$MD_PROD"; fi
if [ -z "$SOLUTE_TYPE" ]; then echo "Guessing SOLUTE_TYPE";  SOLUTE_TYPE="UNKNOWN"; else echo "Considering SOLUTE_TYPE=$SOLUTE_TYPE"; fi
if [ -z "$SIEVE" ]; then echo 'Considering SIEVE=1';  SIEVE=1; else echo "Using SIEVE=$SIEVE"; fi
if [ -z "$DIST" ]; then echo 'Considering DIST=4.0';  DIST="4.0"; else echo "Using DIST=$DIST"; fi
if [ -z "$ANGLE" ]; then echo 'Considering ANGLE=120.0';  ANGLE="120.0"; else echo "Using ANGLE=$ANGLE"; fi
if [ -z "$SUFFIX_MDCRD" ]; then echo 'Considering SUFFIX_MDCRD=_solute.mdcrd';  SUFFIX_MDCRD="_solute.mdcrd"; else echo "Using SUFFIX_MDCRD=$SUFFIX_MDCRD"; fi
if [ -z "$PREFIX_MDCRD" ]; then echo 'Considering PREFIX_MDCRD=md_';  PREFIX_MDCRD="md_"; else echo "Using PREFIX_MDCRD=$PREFIX_MDCRD"; fi
if [ -z "$DO_INF" ]; then DO_INF="YES"; else echo "Considering DO_INF=$DO_INF"; fi
if [ -z "$DO_SURF" ]; then DO_SURF="NO"; else echo "Considering DO_SURF=$DO_SURF"; fi
if [ -z "$STRUCT_DIR" ]; then echo 'Considering STRUCT_DIR=STRUCT';  STRUCT_DIR="STRUCT"; else echo "Using STRUCT_DIR=$STRUCT_DIR"; fi
if [ -z "$MASK" ]; then MASK='NONE'; MASK_READ="NO"; echo "Guessing MASK for RMSD";  else echo "Using MASK=$MASK"; MASK_READ="YES";   fi

if [ -n "$PBS_ENVIRONMENT" ] ; then
  NPROCS=$(cat $PBS_NODEFILE | wc -l)
fi

if [ -z "$NPROCS" ]
then
   NPROCS=$(cat /proc/cpuinfo | grep -c processor)
   echo "Using $NPROCS available processors"
else
   echo "Using $NPROCS processors as predefined "
fi

if [ ${SOLUTE_TYPE} == "APT" ] && [ ${MASK} == "NONE" ]  
then
   MASK="P,OP1,OP2,O5',C5',C4',C3',C2',C1',O4',O3'"
   echo "SOLUTE_TYPE=''APT''   MASK=$MASK"
elif [ ${SOLUTE_TYPE} == "PEP" ] && [ ${MASK} == "NONE" ]
then
   MASK="C,O,CA,N"
   echo "SOLUTE_TYPE=''PEP''   MASK=$MASK"
fi


WORKDIR_TRJ=$PWD

for MOL in $MD_TRAJ 
do
  if [ -e ${MOL}_${MD_TYPE} ]
  then
    cd ${MOL}_${MD_TYPE} 
  else
    echo "${MOL}_${MD_TYPE} does not exist. Exiting!"
    exit
  fi
  if [ ${SUFFIX_MDCRD} == ".mdcrd" ]
  then
     TOPOLOGY=${MOL}.top 
  elif [ ${SUFFIX_MDCRD} == "_solute.mdcrd" ] 
  then
     TOPOLOGY=${MOL}_solute.top 
   
  elif [ ${SUFFIX_MDCRD} == "_solutewat.mdcrd" ] 
  then
     TOPOLOGY=${MOL}_solutewat.top 
  else
     echo "Unexpected SUFFIX_MDCRD=${SUFFIX_MDCRD}"
     if [ -z $TOPOLOGY ]
     then 
       echo "Trying TOPOLOGY=$TOPOLOGY, but maybe it does not work!"
     fi
  fi
  echo "Selecting TOPOLOGY=$TOPOLOGY for SUFFIX_MDCRD=${SUFFIX_MDCRD}"
  if [ ! -e $TOPOLOGY ]
  then 
      echo "$TOPOLOGY not found in ${MOL}_MD. Exiting"
      exit
  fi
  if [ -z $REFPDB ]
  then 
     REFPDB=${WORKDIR_TRJ}/${MOL}_${MD_TYPE}/1.EDITION/${MOL}_solute.pdb 
  fi 
  if [ ! -e $REFPDB ]
  then 
     echo "$REFPDB not found. Exiting."
     exit
  else
     echo "Using REFPDB=$REFPDB"
  fi
  NRES=$(grep 'ATOM  '  $REFPDB | grep -v 'WAT\|Na+\|Cl-' | tail -1 | awk '{print $5}')
  NAT=$(grep -v 'WAT\|Na+\|Cl-' $REFPDB  | grep -c  'ATOM  ' )
  if [ $NRES -gt 10 ]
  then
     IRES=3
     let "JRES=$NRES-2"
  else 
     IRES=1
     JRES=$NRES
  fi
# Getting information from REFPDB to select the proper mask!
  if [ ${MASK_READ} == "NO" ] &&   [ ${SOLUTE_TYPE} == "PEP" ]
  then
     MOL_MASK=":${IRES}-${JRES}@${MASK}"
  elif [ ${MASK_READ} == "NO" ] &&   [ ${SOLUTE_TYPE} == "APT" ]
  then
     MOL_MASK=":${IRES}-${JRES}@${MASK}"
  elif [ ${MASK_READ} == "YES" ] &&   [ ${SOLUTE_TYPE} == "PEP" ]
  then
     MOL_MASK="${MASK}"
  elif [ ${MASK_READ} == "YES" ] &&   [ ${SOLUTE_TYPE} == "APT" ]
  then
     MOL_MASK="${MASK}"
  elif [  ${SOLUTE_TYPE} == "UNKNOWN" ]
  then
     if [ $MASK_READ == "YES" ]
     then
        echo "Guessing SOLUTE_TYPE "
     else
        echo "Guessing both SOLUTE_TYPE and MASK "
     fi
     grep 'ATOM ' $REFPDB | awk '{printf(" %s %3i\n",$4,$5)}' | uniq > tmp.res
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
        echo "Good. This seems to be a nice aptamer"
        SOLUTE_TYPE="APT"
        if [ "$MASK_READ" == "NO" ]
        then
           MASK="P,OP1,OP2,O5',C5',C4',C3',C2',C1',O4',O3'"
           MOL_MASK=":${IRES}-${JRES}@${MASK}"
        else
           MOL_MASK=${MASK} 
        fi
     elif [ $NRES -eq $naa ]
     then
        echo "Good. This seems to be a nice peptide"
        SOLUTE_TYPE="PEP"
        if [ "$MASK_READ" == "NO" ]
        then
          MASK="C,O,CA,N"
          MOL_MASK=":${IRES}-${JRES}@${MASK}"
        else
           MOL_MASK=${MASK} 
        fi
     elif  [ "$MASK_READ" == "YES" ]
     then
        echo "Unknown solute type"
        SOLUTE_TYPE="UNKNOWN"
        MOL_MASK=${MASK}
     else
        echo "Sorry, no luck in guessing solute type "
        echo "Better definne your own MASK for clustering"
        exit
     fi
  else
     echo "No MASK for RMSD/RGYR, no way"
     exit
  fi
  # Auxiliary octave script file 
  MFILE_RMSD_RGYR="$APTAMD/STRUCT/rmsd_rgyr_metrics_DUMMY.m"
  MFILE_RMSD_SURF="$APTAMD/STRUCT/rmsd_surf_metrics_DUMMY.m"
  if [ ${SOLUTE_TYPE} == "PEP" ]
  then
     MFILE_RMSD_INF="$APTAMD/STRUCT/rmsd_inf_peptide_metrics_DUMMY.m"
  elif [ ${SOLUTE_TYPE} == "APT" ]
  then
     MFILE_RMSD_INF="$APTAMD/STRUCT/rmsd_inf_aptamers_metrics_DUMMY.m"
  fi
  if [ ! -e 6.ANALYSIS ]
  then 
     echo "6.ANALYSIS directory not found in ${MOL}_MD. Making it"
     mkdir 6.ANALYSIS 
  fi
  cd 6.ANALYSIS
  if [ ! -e $STRUCT_DIR ]
  then 
      echo "Creating $STRUCT_DIR directory in ${MOL}_MD/6.ANALYSIS"
      mkdir $STRUCT_DIR
  else
      echo "Found $STRUCT_DIR directory in ${MOL}_MD/6.ANALYSIS"
  fi
  cd $STRUCT_DIR
  WORKDIR=$PWD 

  # Printing out options
  cat <<EOF > options.txt
MOL=$MOL     
TOPOLOGY=$TOPOLOGY
NRES=$NRES
SOLUTE_TYPE=$SOLUTE_TYPE

STRUCT_DIR=$STRUCT_DIR
MASK=$MASK
MOL_MASK=$MOL_MASK
REFPDB=$REFPDB
SIEVE=$SIEVE
DIST=$DIST
ANGLE=$ANGLE
DO_INF=$DO_INF
DO_SURF=$DO_SURF

MD_TYPE=$MD_TYPE
MD_PROD=$MD_PROD
SUFFIX_MDCRD=$SUFFIX_MDCRD
PREFIX_MDCRD=$PREFIX_MDCRD

NPROCS=$NPROCS
EOF

# Obtaining HB reference file for the initial structure

  echo "Running cpptraj to generate ${MOL}_ref.dat "

  if [ "${SOLUTE_TYPE}" == "APT" ]
  then 

  echo "Running x3dna to generate ${MOL}_ref.dat and ${MOL}_ref.out_x3dna "

  # Running x3-dna 
  echo 'MODEL        1' >temp_ref.pdb
  grep 'ATOM '  ${REFPDB} >> temp_ref.pdb;
  echo 'ENDMDL'        >> temp_ref.pdb
  $X3DNAHOME/bin/x3dna-dssr --nmr --no-pair --auxfile=no -i=temp_ref.pdb  -o=${MOL}_ref.out_x3dna >& ${MOL}_ref_dssr.log

  #Getting data from base pair interactions 
  grep -n 'List of ' ${MOL}_ref.out_x3dna | grep 'base pairs' | sed 's/:/  /' |awk '{L=$1+$4+1;printf("%i,%ip\n",$1,L)}' > ${MOL}_ref_bp.sed
  sed -n -f ${MOL}_ref_bp.sed ${MOL}_ref.out_x3dna > ${MOL}_ref_bp.dat
  sed -i 's/:D/   D/g' ${MOL}_ref_bp.dat
  grep List ${MOL}_ref_bp.dat  | awk '{print $3}' > ${MOL}_ref_bp_n.dat
  grep -v List ${MOL}_ref_bp.dat | grep -v DSSR   | awk '{print $3, $5}' | sed 's/3\///' | sed 's/5\///' | sed 's/D[A-T]/   /g' > ${MOL}_ref_bp_ij.dat

  #Getting data for non-pairing interactions, either direct base contacts or base stacking 
  grep -n 'List of ' ${MOL}_ref.out_x3dna | grep 'non-pairing interactions' | sed 's/:/  /' |awk '{l=$1+$4;printf("%i,%ip\n",$1,l)}' > ${MOL}_ref_np.sed
  sed -n -f ${MOL}_ref_np.sed ${MOL}_ref.out_x3dna > ${MOL}_ref_np.dat
  sed -i 's/:D/   D/g' ${MOL}_ref_np.dat
  grep 'List of' ${MOL}_ref_np.dat  | awk '{print $3}' > ${MOL}_ref_np_n.dat

  sed -i 's/interBase/ 0 interBase/g' ${MOL}_ref_np.dat
  sed -i 's/stacking:/ 1 stacking/g' ${MOL}_ref_np.dat
  grep -v List ${MOL}_ref_np.dat | awk '{print $3, $5, $6}' | sed 's/3\///' | sed 's/5\///' | sed 's/D[A-T]/   /g' > ${MOL}_ref_np_ij.dat

  cat ${MOL}_ref_bp_n.dat ${MOL}_ref_bp_ij.dat  ${MOL}_ref_np_n.dat ${MOL}_ref_np_ij.dat >  ${MOL}_ref.dat

  rm -f ${MOL}_ref_bp_n.dat ${MOL}_ref_bp_ij.dat  ${MOL}_ref_np_n.dat ${MOL}_ref_np_ij.dat  ${MOL}_ref_np.sed ${MOL}_ref_bp.sed ${MOL}_ref_dssr.log  ${MOL}_ref_bp.dat ${MOL}_ref_np.dat temp_ref.pdb

  fi

  # Preparing SCRIPT file for running MULTIPDB cpptraj calculations 
  # Environment file
  cat <<EOF >  environment.sh
export  MOL="${MOL}"
export  NRES="${NRES}"
export  NAT="${NAT}"
export  WORKDIR="${WORKDIR}"
export  AMBERHOME="${AMBERHOME}"
export  OCTAVE="${OCTAVE}"
export  X3DNAHOME="${X3DNAHOME}"
export  MFILE_RMSD_RGYR="${MFILE_RMSD_RGYR}"
export  MFILE_RMSD_INF="${MFILE_RMSD_INF}"
export  SCRATCH="${SCRATCH}"
export  TOPOLOGY="${TOPOLOGY}"
export  REFPDB="${REFPDB}"
export  MOL_MASK="${MOL_MASK}"
export  DIST="${DIST}"
export  ANGLE="${ANGLE}"
export  SIEVE="${SIEVE}"
export  MD_PROD="${MD_PROD}"
export  DO_INF="${DO_INF}"
export  DO_SURF="${DO_SURF}"
EOF

  rm -f TASK.sh 
 
  I=0
  for file in $(ls ../../${MD_PROD}/${PREFIX_MDCRD}???${SUFFIX_MDCRD}  )
  do 
    let "I=$I+1"
    if [ $I -lt 10 ]
    then
      txt="00${I}"
    elif [ $I -lt 100 ]
    then
      txt="0${I}"
    else 
      txt=${I}
    fi
    if [ ${SOLUTE_TYPE} == "PEP" ]
    then
       echo "cd $WORKDIR; $APTAMD/STRUCT/aux_struct_peptide.sh $file md_${txt}  > $SCRATCH/tmp_${txt}.log; rm -f  $SCRATCH/tmp_${txt}.log" >> TASK.sh
    elif [ ${SOLUTE_TYPE} == "APT" ]
    then
       echo "cd $WORKDIR; $APTAMD/STRUCT/aux_struct_aptamers.sh $file md_${txt}  > $SCRATCH/tmp_${txt}.log; rm -f  $SCRATCH/tmp_${txt}.log" >> TASK.sh
    fi
  done
  echo ${SOLUTE_TYPE}
  echo "Running parallel task in $PWD"
  cat TASK.sh  | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS  

  # Care is needed here. X3DNA multipdb calcs may be incomplete.
  # Hence we make sure that RMSD and RGYR data are paired with INF data 
  if [ $DO_INF == "YES" ]
  then 
     rm -f RMSD.dat RGYR.dat SURF.dat INF.dat 
     for file in $(ls md_*.inf)
     do
       id=${file%%.*} 
       nl_inf=$(cat $file | grep -v '#' | wc -l)
       nl_rmsd=$(cat ${id}.rmsd | grep -v '#' | wc -l) 
       if [ ${nl_inf} -ne ${nl_rmsd} ]
       then
           echo "${id}.inf and ${id}.rmsd data files have different size."
          echo "we get only ${nl_inf} lines out of ${nl_rmsd} in ${id}.rmsd" 
        fi
       grep -v '#' ${id}.inf  | awk -F ',' '{print $NF}' >> INF.dat 
       grep -v '#' ${id}.rmsd | head -${nl_inf} | awk '{print $2}' >> RMSD.dat
       grep -v '#' ${id}.rgyr | head -${nl_inf} | awk '{print $2}' >> RGYR.dat
       if [ $DO_SURF == "YES" ]; then grep -v '#' ${id}.surf | head -${nl_inf} | awk '{print $2}' >> SURF.dat; fi
     done
     $OCTAVE --no-gui  -q $APTAMD/STRUCT/rmsd_rgyr_surf_plot.m
     $OCTAVE --no-gui  -q $APTAMD/STRUCT/rmsd_inf_plot.m
  else
     cat md_*.rmsd | grep -v '#' | awk '{print $2}' > RMSD.dat
     cat md_*.rgyr | grep -v '#' | awk '{print $2}' > RGYR.dat
     if [ $DO_SURF == "YES" ]; then grep -v '#' ${id}.surf | head -${nl_inf} | awk '{print $2}' >> SURF.dat; fi
     $OCTAVE --no-gui  -q $APTAMD/STRUCT/rmsd_rgyr_surf_plot.m
  fi

  cd $WORKDIR_TRJ 

done

