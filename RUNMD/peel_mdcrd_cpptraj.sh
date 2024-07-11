#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
source $APTAMD/ENV/geuo_env.sh

# parmed requires AMBER PYTHON environment
export PYTHON=$AMBERHOME/miniconda/bin
export PYTHONPATH=$AMBERHOME/lib/python3.9/site-packages/

# Topology file 
TOP=$1
# INPUT Coord set: it can be either a coord filename or list filename
COORD=$2
# Thickness of the water shell
RPEEL=$3
# SIEVE for snapshot processing 
SIEVE=$4
# NWAT_INP (not required, but optional )
NWAT_INP=$5

if [ -z "$TOP" ] ||  [ -z "$COORD" ] ||  [ -z "$RPEEL" ] ||  [ -z "$SIEVE" ]
then
   echo " "
   echo "Usage: peel_mdcrd_cpptraj.sh topology_file coord_file rpeel sieve [nwat] "
   echo "      or  "
   echo "Usage: peel_mdcrd_cpptraj.sh topology_file coord_filelist rpeel sieve [nwat] "
   echo " "
   echo " Warning: nwat can be selected automatically basing on the rpeel value."
   echo " If nwat is given then rpeel is ignored. Be careful then."
   echo " "
   echo "Examples: "
   echo " " 
   echo " peel_mdcrd_cpptraj.sh  apt.top  md_001.mdcrd  12.0  50 "
   echo " " 
   echo " peel_mdcrd_cpptraj.sh  apt.top  LISTA 12.0 1 "
   echo " " 
   echo "    where  LISTA can be a filelist like"
   echo " "
   echo "../../md_1.mdcrd "
   echo "../../md_2.mdcrd "
   echo "etc"
   exit
fi

if [ -z "$NWAT_INP" ]
then
   NWAT_INP=0
else
   echo "Selecting $NWAT_INP water molecules around the solute molecule"
fi

IDMOL=$(basename $TOP)
IDMOL=${IDMOL%%.*}_solutewat

if [ -z "$NPROCS" ]
then
   NPROCS=$(cat /proc/cpuinfo | grep -c processor)
   echo "Using $NPROCS available processors"
else
   echo "Using $NPROCS processors as predefined "
fi
export OMP_NUM_THREADS=$NPROCS


if [ -z "$DATAMASH" ]; then echo 'datamash is not available, but needed'; exit; fi

TEST_COORD_FILE=$(echo $COORD | awk '{print $1}')

ILIST=$(file $TEST_COORD_FILE | grep -i -c ASCII)
if [ $ILIST -gt 0 ]
then
   echo "$TEST_COORD_FILE assumed to be a list of COORD sets"
else
  IDATA=$(file $TEST_COORD_FILE |  grep -i -c data)
  if [ $IDATA -gt 0 ]
  then
     echo "$TEST_COORD_FILE seems to be a COORD set"
  else
     echo "$TEST_COORD_FILE unknown type"
     exit
  fi
fi

# Mask of atoms to search for closest waters around
# Probably, this shouldn't be changed
if [ -z $MASK ]
then
   MASK="!:WAT,Na+,Cl-"
fi

# Temp directory
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

if [ $NWAT_INP -eq 0 ]
then

I=0

rm -f input.trajin
for ((ISET=1;ISET<=$NSET;ISET++))
do
  COORD=$(sed -n "${ISET},${ISET}p" LISTA_COORD)
  FILE=$(echo $COORD | awk '{print $1}')
  NSNAP=$(ncdump -h $WORKDIR/${FILE} |grep frame | grep UNLI | awk '{print $6}' | sed 's/(//')
  echo  "trajin $WORKDIR/$FILE 1 $NSNAP $SIEVE " >> input.cpptraj
done

cat <<EOF >> input.cpptraj
trajout tmp ncrestart
watershell $MASK lower $RPEEL upper $RPEEL out tmp_WS.dat
go
EOF
$AMBERHOME/bin/cpptraj.OMP $WORKDIR/$TOP < input.cpptraj > output.cpptraj 

sed -i '1,1d' tmp_WS.dat
ncoord=$(cat tmp_WS.dat | wc -l)
maxWS=$(awk '{print $2}' tmp_WS.dat | sort -k 1 -n| tail -1)
minWS=$(awk '{print $2}' tmp_WS.dat | sort -k 1 -n| head -1)

echo minWS=$minWS
echo maxWS=$maxWS

medianWS=$(awk '{print $2}' tmp_WS.dat | $DATAMASH median 1)
medianWS=${medianWS%%.*} 

echo medianWS=$medianWS

#Maximum number of water molecules to be extracted
echo  "nw=${medianWS} ;" >  temp.m
echo  'mw=[1E6 1E5 1E4 1E3 1E2 1E1] ;' >> temp.m
echo  'lw=length(mw) ;' >> temp.m
echo  'for i=1:lw' >> temp.m
echo  '    if  nw >  mw(i)  ' >> temp.m
echo  '        uw=floor(nw/mw(i))*mw(i)*1.20;' >> temp.m
echo  '        printf(" %d ",uw);' >> temp.m
echo  '        exit' >> temp.m
echo  '    end ' >> temp.m
echo  'end ' >> temp.m
echo  'uw=nw ;' >> temp.m
echo  'printf(" %d ",uw);' >> temp.m
echo  'exit' >> temp.m

NWAT=$($OCTAVE -q temp.m)
echo "NWAT to be extracted = $NWAT"

else

NWAT=$NWAT_INP

fi

# Preparing TOPOLOGY FILE 
REFCRD=$(head -1 LISTA_COORD | awk '{print $1}')
REFPDB=temp.pdb 
$AMBERHOME/bin/cpptraj $WORKDIR/$TOP <<EOF
trajin  ${WORKDIR}/${REFCRD} 1 1 1
trajout ${REFPDB} pdb 
autoimage origin 
go
EOF

IPROT=$(grep 'ATOM  ' $REFPDB | grep -v 'WAT' | grep -v 'HOH' | grep -v 'Na+' | grep -v 'Cl-' | tail -1 | awk '{print $5}')
ISOLV=$(grep 'ATOM  ' $REFPDB | grep  'WAT\|HOH'  | head -1 | awk '{print $5}')
NTOPWAT=$(grep -c 'O   WAT' $REFPDB)

let " ARES = $ISOLV - 1  + $NWAT + 1 "
let " BRES = $ISOLV - 1  + $NTOPWAT "

$AMBERHOME/bin/parmed -n  $WORKDIR/$TOP <<EOF
strip :${ARES}-${BRES}
parmout solutewat.top
go
EOF
$AMBERHOME/bin/cpptraj  solutewat.top <<EOF
parmbox nobox  
parmwrite out solutewat_nobox.top
go
EOF
mv solutewat_nobox.top $WORKDIR/${IDMOL}.top 
rm -f solutewat.top

echo "# IPROT ISOLV  NWAT  RPEEL " > PEEL.dat
echo " $IPROT $ISOLV $NWAT $RPEEL" >> PEEL.dat
mv PEEL.dat $WORKDIR/${IDMOL}.info

# Loop over COORD sets
rm -f TASK.sh
I=0
for ((ISET=1;ISET<=$NSET;ISET++))
do
  COORD=$(sed -n "${ISET},${ISET}p" LISTA_COORD)
  FILE=$(echo $COORD | awk '{print $1}')
  FILE=${FILE%%.*} 
  echo  "trajin $WORKDIR/$COORD " > input.cpptraj

  let "I=$I+1"

# Note that the closest command selects water around solute
# atoms only, but preserves the Na+/Cl- counterions
cat <<EOF >> input.cpptraj
trajout $WORKDIR/${FILE}_solutewat.mdcrd netcdf 
autoimage origin anchor $MASK
closest $NWAT $MASK solventmask :WAT oxygen 
go
EOF
  mv input.cpptraj input.${I}

  echo "$AMBERHOME/bin/cpptraj $WORKDIR/$TOP < input.${I} > output.${I}; rm -f tmp.${I}" >> TASK.sh

done

echo "Processing all files in parallel" 
cat TASK.sh | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS   >/dev/null 2>&1

cd $WORKDIR

rm -r -f $TMPDIR

