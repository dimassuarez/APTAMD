#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
source $APTAMD/ENV/aptamd_env.sh

# PDB files with a water shell around selected solute atoms 

# Topology file 
TOP=$1
# INPUT Coord set: it can be either a coord filename or list filename
COORD=$2
# Thickness of the water shell
RPEEL=$3
# Solute Mask: Default is !:WAT,Na+,Cl-
MASK=$4
# ID for the PDB files to be printed ouf
IDMOL=$5
# Initial numbering for PDBs
INIT=$6 

if [ -z "$TOP" ] ||  [ -z "$COORD" ] ||  [ -z "$RPEEL" ]
then
   echo "Usage: peel_cpptraj topology_file coord_file rpeel [solute_mask]  [idmol]  [init]  "
   echo "      or  "
   echo "Usage: peel_cpptraj topology_file coord_filelist rpeel [solute_mask]  [idmol]  [init] "
   echo " "
   echo "Example: "
   echo " peel_cpptraj.sh  apt.top  LISTA 12.0   ":1-40,Na+"  apt 1000 "
   echo " "
   echo "    where  LISTA is a filelist like"
   echo " "
   echo "../../md_1.mdcrd  1 20000 20"
   echo "../../md_2.mdcrd  1 20000 20"
   echo "etc"
   exit
fi

if [ -z "$NPROCS" ]
then
   NPROCS=$(cat /proc/cpuinfo | grep -c processor)
   echo "Using $NPROCS available processors"
else
   echo "Using NPROCS=$NPROCS processors as predefined "
fi
export OMP_NUM_THREADS=$NPROCS

TEST_COORD_FILE=$(echo $COORD | awk '{print $1}')

ILIST=$(file $TEST_COORD_FILE | grep -c ASCII)
if [ $ILIST -gt 0 ]
then
   echo "$TEST_COORD_FILE assumed to be a list of COORD sets"
else
  IDATA=$(file $TEST_COORD_FILE | awk '{print $NF}' | grep -c data)
  if [ $IDATA -gt 0 ]
  then
     echo "$TEST_COORD_FILE seems to be a COORD set"
  else
     echo "$TEST_COORD_FILE unknown type"
     exit
  fi
fi

if [ -z $MASK ]
then
   MASK="!:WAT,Na+,Cl-"
fi
if [ -z $IDMOL ]
then
   IDMOL="md_peel"
fi
if [ -z $INIT ]
then
   INIT=0
fi

# TEMP Directories
if [ ! -n "$PBS_ENVIRONMENT" ] ; then
   SCRATCH=/scratch
fi
TT=$(date +%N)
TMPDIR=${SCRATCH}/TMPDIR_${TT}

WORKDIR=$PWD
mkdir $TMPDIR

if [ $ILIST -gt 0 ] 
then 
  cp $TEST_COORD_FILE $TMPDIR/LISTA_COORD
else
  echo "$COORD" > $TMPDIR/LISTA_COORD
fi

cd $TMPDIR

NSET=$(cat LISTA_COORD | wc -l)

rm -f input.cpptraj
for ((ISET=1;ISET<=$NSET;ISET++))
do
  COORD=$(sed -n "${ISET},${ISET}p" LISTA_COORD)
  echo  "trajin $WORKDIR/$COORD " >> input.cpptraj
done

IPEEL=${RPEEL/.*}

if [ $IPEEL -eq 0 ]; then echo "PDBs will contain only solute atoms"; fi

cat <<EOF >> input.cpptraj
trajout tmp ncrestart
watershell $MASK lower $RPEEL upper $RPEEL out tmp_WS.dat
go
EOF

$AMBERHOME/bin/cpptraj.OMP $WORKDIR/$TOP < input.cpptraj > output.cpptraj 

if [ $IPEEL -gt 0 ]
then 
  sed -i '1,1d' tmp_WS.dat
  ncoord=$(cat tmp_WS.dat | wc -l)
  maxWS=$(awk '{print $2}' tmp_WS.dat | sort -k 1 -n| tail -1)
  minWS=$(awk '{print $2}' tmp_WS.dat | sort -k 1 -n| head -1)
  echo minWS=$minWS
  echo maxWS=$maxWS
else
  ncoord=$(ls tmp* | wc -l)
fi

rm -f TASK.sh
I=$INIT
for ((icoord=1;icoord<=ncoord;icoord++))
do
   let "I=$I+1"
   if [ $I -lt 10 ]
   then
      txt="0000${I}"
   elif [ $I -lt 100 ]
   then
      txt="000${I}"
   elif [ $I -lt 1000 ]
   then
      txt="00${I}"
   elif [ $I -lt 10000 ]
   then
      txt="0${I}"
   else
      txt="${I}"
   fi
   if [ $ncoord -eq 1 ] ; then mv tmp tmp.1;fi

   if  [ $IPEEL -gt  0 ]
   then

   NWAT=$(sed -n "${icoord},${icoord}p" tmp_WS.dat | awk '{print $2}')
cat <<EOF > input.${icoord}
trajin  tmp.${icoord}
trajout $WORKDIR/${IDMOL}_${txt}.pdb pdb vdw include_ep
autoimage origin anchor $MASK
closest $NWAT $MASK solventmask :WAT@O,Na+,Cl- 
go
EOF

   else 

cat <<EOF > input.${icoord}
trajin  tmp.${icoord}
trajout $WORKDIR/${IDMOL}_${txt}.pdb pdb vdw include_ep
strip   :WAT,Na+,Cl- 
autoimage origin anchor $MASK
go
EOF

   fi

   echo "$AMBERHOME/bin/cpptraj $WORKDIR/$TOP < input.${icoord} > output.${icoord}; rm -f tmp.${icoord}" >> TASK.sh
   
done

echo "Processing all files in parallel" 
cat TASK.sh | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS   >/dev/null 2>&1

cd $WORKDIR

rm -r -f $TMPDIR
