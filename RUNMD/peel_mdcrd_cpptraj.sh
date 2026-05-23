#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
source $APTAMD/ENV/aptamd_env.sh

# Topology file 
TOPOLOGY=$1
# INPUT Coord set: it can be either a coord filename or list filename
COORD=$2
# Thickness of the water shell
RPEEL=$3
# SIEVE for snapshot processing 
SIEVE=$4
# NWAT_INP (not required, but optional )
NWAT_INP=$5

if [ -z "$TOPOLOGY" ] ||  [ -z "$COORD" ] ||  [ -z "$RPEEL" ] ||  [ -z "$SIEVE" ]
then
   echo " "
   echo "Usage: peel_mdcrd_cpptraj.sh topology_file coord_file rpeel sieve [nwat] "
   echo "      or  "
   echo "Usage: peel_mdcrd_cpptraj.sh topology_file coord_filelist rpeel sieve [nwat] "
   echo " "
   echo " Warning: nwat can be selected automatically basing on the rpeel value."
   echo " If nwat is given then rpeel is ignored. Be careful then."
   echo " "
   echo "         If nwat = 0 then only solute atoms are preserved" 
   echo "         and sieve is ignored"
   echo " "
   echo "         If nwat = N then only counterions and N waters are preserved" 
   echo "         and sieve is take into account"
   echo " "
   echo "         If nwat = '0 N' then two files are generated without and with " 
   echo "         N-waters/counterions"
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

SOLUTE=0       # Do not extract 
SOLUTE_ONLY=0
if [ -z "${NWAT_INP}" ]
then
   NWAT=-1
   echo "Selecting waters to be extracted using peel=$RPEEL"
else
  icheck=0
  for IWAT in ${NWAT_INP}
  do
       let "icheck=$icheck+1"
       if [ ${IWAT} -eq 0 ]; then SOLUTE=1; fi 
       if [ ${SOLUTE} -eq 0 ] ; then NWAT=$IWAT; fi
       if [ ${SOLUTE} -eq 1 ] && [ $icheck -eq 2 ] ; then NWAT=$IWAT; fi
  done
  if [ $icheck -gt 2 ]
  then
     echo "Option= ${NWAT_INP} cannot be processed!" 
     exit
  fi
  if [ $icheck -eq 2 ] && [ $SOLUTE -eq 0 ] 
  then
     echo "Option= ${NWAT_INP} cannot be processed!"
     exit
  fi
  if [ $icheck -eq  1 ] && [ $SOLUTE -eq 1 ]
  then 
     echo "Only solute atoms are preserved "
     SOLUTE_ONLY=1
  elif [ $icheck -eq  1 ] && [ $SOLUTE -eq 0 ]
  then 
     echo "Selecting $NWAT water molecules around the solute molecule"
  elif [ $icheck -eq  2 ] && [ $SOLUTE -eq 1 ]
  then 
     echo "Selecting solute atoms and selecting $NWAT water molecules"
     echo "around the solute molecule in separate files"
  else
     echo "No idea what to do"
     exit
  fi
fi

IDMOL=$(basename $TOPOLOGY)
IDMOL=${IDMOL%%.*}

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

if [ $NWAT -lt 0 ] && [ $SOLUTE_ONLY -eq 0 ]
then

I=0

rm -f input.trajin
for ((ISET=1;ISET<=$NSET;ISET++))
do
  COORD=$(sed -n "${ISET},${ISET}p" LISTA_COORD)
  FILE=$(echo $COORD | awk '{print $1}')
  echo  "trajin $WORKDIR/$FILE 1 last $SIEVE " >> input.cpptraj
done

cat <<EOF >> input.cpptraj
trajout tmp ncrestart
watershell $MASK lower $RPEEL upper $RPEEL out tmp_WS.dat
go
EOF
$AMBERHOME/bin/cpptraj.OMP $WORKDIR/$TOPOLOGY < input.cpptraj > output.cpptraj 

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

echo "NWAT to be extracted (defined in input) =$NWAT"

fi

if [ $SOLUTE_ONLY -eq 0 ]
then 

# Preparing TOPOLOGY FILE 
REFCRD=$(head -1 LISTA_COORD | awk '{print $1}')
REFPDB=temp.pdb 
$AMBERHOME/bin/cpptraj $WORKDIR/$TOPOLOGY <<EOF
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

$AMBERHOME/bin/parmed -n  $WORKDIR/$TOPOLOGY <<EOF
strip :${ARES}-${BRES}
parmout solutewat.top
go
EOF
$AMBERHOME/bin/cpptraj  solutewat.top <<EOF
parmbox nobox  
parmwrite out solutewat_nobox.top
go
EOF
mv solutewat_nobox.top $WORKDIR/${IDMOL}_solutewat.top 
rm -f solutewat.top

echo "# IPROT ISOLV  NWAT  RPEEL " > PEEL.info
echo " $IPROT $ISOLV $NWAT $RPEEL" >> PEEL.info
mv PEEL.info $WORKDIR/${IDMOL}_solutewat.info

# Loop over COORD sets
rm -f TASK.sh

I=0
for ((ISET=1;ISET<=$NSET;ISET++))
do
  COORD=$(sed -n "${ISET},${ISET}p" LISTA_COORD)
  FILE=$(echo $COORD | awk '{print $1}')
  FILE=${FILE%%.*} 
  echo  "trajin $WORKDIR/$COORD 1 last ${SIEVE} " > input.cpptraj

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

  echo "$AMBERHOME/bin/cpptraj $WORKDIR/$TOPOLOGY < input.${I} > output.${I}; rm -f tmp.${I}" >> TASK.sh

done

echo "Processing all MDCRD files in parallel" 
cat TASK.sh | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS   >/dev/null 2>&1

fi # ENDIF of SOLUTE_ONLY condition


if [ $SOLUTE -eq  1 ]
then 

rm -f TASK.sh

$AMBERHOME/bin/parmed -n  $WORKDIR/$TOPOLOGY <<EOF
strip :WAT,Na+,Cl-
parmout solute.top
go
EOF
mv solute.top ${IDMOL}_solute.top 
cp ${IDMOL}_solute.top $WORKDIR/$OUTDIR/

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
trajout $WORKDIR/${FILE}_solute.mdcrd netcdf 
strip :WAT,Na+,Cl-
go
EOF
  mv input.cpptraj input.${I}

  echo "$AMBERHOME/bin/cpptraj $WORKDIR/$TOPOLOGY < input.${I} > output.${I}; rm -f tmp.${I}" >> TASK.sh

done

echo "Processing all MDCRD files in parallel" 
cat TASK.sh | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS   >/dev/null 2>&1

fi  # ENDIF of SOLUTE condition 

cd $WORKDIR

rm -r -f $TMPDIR

