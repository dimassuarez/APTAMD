#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
if [ "$#" -eq 0 ]; then more $APTAMD/DOC/do_sconform.txt ; exit; fi

source $APTAMD/ENV/aptamd_env.sh

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
# IFRAG=0  full solute
    IFRAG=$3
# DO_CC_MLA
    DO_CC_MLA=$4
fi

if [ ! -z "$MOL" ] && [  -z "$MD_TRAJ" ]; then MD_TRAJ="${MOL}"; echo $MD_TRAJ; unset MOL; fi

if [ -n "$PBS_ENVIRONMENT" ] ; then
  NPROCS=$(cat $PBS_NODEFILE | wc -l)
fi

if [ -z "$NPROCS" ]
then
   export NPROCS=$(cat /proc/cpuinfo | grep -c processor)
   echo "Using $NPROCS available processors"
else
   echo "Using $NPROCS processors as predefined "
fi

#                       done and then only CC_MLA is performed

PRSTOP="NO"
if [ -z "$MD_TRAJ" ]; then PRSTOP="YES";  fi
if [ -z "$MD_TYPE" ]; then PRSTOP="YES";  fi
if [ -z "$IFRAG" ]; then IFRAG=0;  fi
if [ -z "$DO_CC_MLA" ]; then DO_CC_MLA=0 ; fi 
if [ $PRSTOP == "YES"  ]; then more $APTAMD/DOC/do_sconform.txt ; exit; fi 
#GETTOR options
if [ -z "$GETTOR_OPT" ]
then
     GETTOR_OPT="-puck -noMet -noArom"
     echo "Using default GETTOR_OPT=$GETTOR_OPT"
else
     echo "Using GETTOR_OPT=$GETTOR_OPT"
fi
# CCMLA calcs
if [ $DO_CC_MLA -eq 0 ] || [ $DO_CC_MLA == "NO" ]  
then
  DO_CC_MLA="0"
  echo 'CC_MLA calcs are not requested. Only first-order entropies will be computed.'
elif [ $DO_CC_MLA -eq 1 ] || [ $DO_CC_MLA == "YES" ]
then
  DO_CC_MLA="1"
  echo ' Both CC_MLA and first-order entropies will be computed.'
elif [ $DO_CC_MLA -eq 2 ] || [ $DO_CC_MLA == "ONLY" ] 
then
  DO_CC_MLA="2"
  echo 'CC_MLA entropies will be computed using an already existing MATRIX.dat.'
fi
if [ -z "$CUTOFF" ]; then echo 'Choosing CUTOFF=-1'; CUTOFF="-1"; else echo "Using CUTOFF=$CUTOFF"; fi
if [ -z "$COMPOSITE" ]; then echo 'Choosing COMPOSITE=NO'; COMPOSITE="0";  \
   elif [ $COMPOSITE == "YES" ] || [ $COMPOSITE -eq 1 ];  then COMPOSITE=1; echo "Using CC_MLA COMPOSITE method"; fi
if [ -z "$PERCEN" ]; then echo 'Using PERCEN=0 (whole data set)'; PERCEN="0"; else echo "Using PERCEN=$PERCEN"; fi
if [ -z "$SUFFIX_MDCRD" ]; then echo 'Considering SUFFIX_MDCRD=_solute.mdcrd';  SUFFIX_MDCRD="_solute.mdcrd"; else echo "Using SUFFIX_MDCRD=$SUFFIX_MDCRD"; fi
if [ -z "$PREFIX_MDCRD" ]; then echo 'Considering PREFIX_MDCRD=md_';  PREFIX_MDCRD="md_"; else echo "Using PREFIX_MDCRD=$PREFIX_MDCRD"; fi
if [ -z "$MD_PROD" ]; then MD_PROD="5.PRODUCTION"; else echo "Assuming MD_PROD=$MD_PROD"; fi
if [ -z "$SCONFORM_DIR" ]; then echo 'Considering SCONFORM_DIR=SCONFORM'; SCONFORM_DIR="SCONFORM"; fi
if [ -z "$SOLUTE_ALIAS" ]; then echo 'Considering SOLUTE_ALIAS=solute'; SOLUTE_ALIAS="solute" ; fi 

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
  TOPOLOGY=${MOL}_${SOLUTE_ALIAS}.top 
  if [ ! -e $TOPOLOGY ]
  then 
      echo "$TOPOLOGY not found in ${MOL}_MD. Exiting"
      exit
  fi
  SOLUTE_CRD=${MOL}_${SOLUTE_ALIAS}.crd 
  SOLUTE_TOP=${MOL}_${SOLUTE_ALIAS}.top
  SOLUTE_PDB=${MOL}_${SOLUTE_ALIAS}.pdb
  if [ ! -e 1.EDITION/$SOLUTE_CRD ] || [ ! -e 1.EDITION/$SOLUTE_TOP ]
  then 
     echo "$SOLUTE_CRD or $SOLUTE_TOP not found in ${MOL}_MD/1.EDITION. Exiting"
     exit
  else
     cd 1.EDITION
     $AMBERHOME/bin/cpptraj $SOLUTE_TOP <<EOF > mlog 
trajin $SOLUTE_CRD
trajout $SOLUTE_PDB  pdb pdbatom
go
EOF
     NRES=$(grep 'ATOM  ' $SOLUTE_PDB   | tail -1 | awk '{print $5}')
     NFRAG=$(grep -c 'TER' $SOLUTE_PDB)
     declare -a IRES=""
     declare -a JRES=""
     JRES=($(grep TER $SOLUTE_PDB | awk '{print $4}'))
     IRES[0]="1"
     for ((i=1;i<=NFRAG-1;i++))
     do  
       let "j=$i-1"
       let "resid=${JRES["$j"]}+1"
       IRES["$i"]=$resid
     done
     echo "Detected $NFRAG fragments"
     for ((i=0; i<=NFRAG-1;i++))
     do
          echo  "Fragment :${IRES["$i"]}-${JRES["$i"]}"
     done
     cd ../
  fi
  if [ ! -e 6.ANALYSIS ]
  then 
     echo "6.ANALYSIS directory not found in ${MOL}_MD."
     exit             
  fi
  cd 6.ANALYSIS
  if [ $IFRAG -lt 0 ] || [ $IFRAG -gt $NFRAG ]
  then 
     echo "IFRAG=$IFRAG wrong!"
     exit
  fi
  if [ $IFRAG -gt 0 ]
  then
     SCONFORM_DIR_WORK=${SCONFORM_DIR}_FRAG_${IFRAG}
  else
     SCONFORM_DIR_WORK=${SCONFORM_DIR}
  fi
  if [ -e  $SCONFORM_DIR_WORK ]
  then 
      echo "$SCONFORM_DIR_WORK found in ${MOL}_MD/6.ANALYSIS"
  else
      mkdir $SCONFORM_DIR_WORK 
  fi 
  if [ $IFRAG -eq 0 ]
  then 
    echo "Computing Sconform for full solute atoms"
    MOL_MASK=":1-${NRES}"
  else
    echo "Computing Sconform for fragment=$IFRAG"
    let "I=$IFRAG-1"
    MOL_MASK=":${IRES["$I"]}-${JRES["$I"]}"
  fi
  echo $MOL_MASK
  cd $SCONFORM_DIR_WORK
  if [ $DO_CC_MLA -eq 2 ] && [ ! -e MATRIX.dat ]
  then 
      echo "DO_CC_MLA=2 requested, but MATRIX.dat not found"
      echo "Probably, DO_CC_MLA=1 needed."
      exit
  fi 

  WORKDIR=$PWD 
  cp $APTAMD/ENTROPY/run_sconform.sh .
  chmod 755 run_sconform.sh

cat <<EOF > job_sconform.sh
  env NPROCS=$NPROCS  MOL_MASK="$MOL_MASK" REFTOP="${WORKDIR_TRJ}/${MOL}_${MD_TYPE}/${TOPOLOGY}" \
  TRAJDIR="${WORKDIR_TRJ}/${MOL}_${MD_TYPE}/${MD_PROD}/" PREFIX_MDCRD="${PREFIX_MDCRD}" SUFFIX_MDCRD="${SUFFIX_MDCRD}" \
  DO_CC_MLA=$DO_CC_MLA DO_COMP_CC_MLA=$COMPOSITE  CUTOFF="$CUTOFF" PERCEN="$PERCEN"  \
  SCRATCH="$SCRATCH" GETTOR_OPT="$GETTOR_OPT"  ./run_sconform.sh > mlog
EOF
  chmod 755 job_sconform.sh

  ./job_sconform.sh 

  cd $WORKDIR_TRJ 

done


