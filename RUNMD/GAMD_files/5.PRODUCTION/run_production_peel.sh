#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
source $APTAMD/ENV/geuo_env.sh

#Local directory
WORKDIR=$PWD

# TOPOLOGY pathfile which is the script argument should be relative to WORKDIR
TOPOLOGY=$1
INI=$2
NUM=$3                  

if [ -z $TOPOLOGY ]; then  echo "TOPOLOGY not defined"; exit; fi
if [ -z $INI ]; then echo "ID of initial rst file not defined"; exit; fi
if [ -z $NUM ]; then echo "Number of MD production runs not defined"; exit; fi

if [ -z $SIEVE  ]; then SIEVE=20 ; fi
if [ -z $PEEL ]; then PEEL=14 ; fi
if [ -z $USE_GPU ]; then USE_GPU='YES' ; fi

if [ -z "$NPROCS" ]
then
     NPROCS=$(cat /proc/cpuinfo | grep -c processor)
     echo "Using $NPROCS available processors"
else
     echo "Using $NPROCS processors as predefined "
fi
export OMP_NUM_THREADS=$NPROCS

# Probably, this should not be changed 
if [ -z $SOLUTE_MASK ]; then SOLUTE_MASK='!:WAT,Na+,Cl-' ; fi

MD="md"               # alias of MD files
INPUT="md_npt_gamd_restart.inp"    # GAMD input file 

#Local directory
WORKDIR=$PWD

# Peeling info
if [ $PEEL -gt 0 ] && [ ! -e PEEL.info ] 
then
   echo "PEEL=$PEEL, but PEEL.info does not exist"
   exit
fi

if [ -e PEEL.info ] 
then
   MAXWAT=$(tail -1 $WORKDIR/PEEL.info | awk '{print $3}')
fi

# TEMP Directories
TT=$(date +%N)
TMPDIR=${SCRATCH}/TMPDIR_${TT}

# All work is done in a temporal directory
mkdir $TMPDIR

# GAMD equilibration data
if [ -e gamd-restart.dat ]; then cp gamd-restart.dat $TMPDIR/; \
else echo 'gamd-restart missing, but needed'; exit ; fi

for ((JOB = INI ; JOB <= NUM-1; JOB++))
do

if [ -e "STOP" ]
then
   exit
fi 

OLD=${JOB}
NEW=$(expr $OLD + 1)

if [ $NEW -lt 10 ]
then
      txt="00${NEW}"
elif [ $NEW -lt 100 ]
then
      txt="0${NEW}"
else
      txt="${NEW}"
fi

if [ $OLD -lt 10 ]
then
      txt0="00${OLD}"
elif [ $OLD -lt 100 ]
then
      txt0="0${OLD}"
else
      txt0="${OLD}"
fi

CRD="${MD}_${txt0}.rst"
RST="${MD}_${txt}.rst"
OUT="${MD}_${txt}.out"
MDEN="${MD}_${txt}.mden"
MDCRD="${MD}_${txt}.mdcrd"
GAMDLOG="${MD}_${txt}.gamd_log"

# MD is done on the TMPDIR
rm -f mdinfo mdout restrt 
ln -s $TMPDIR/$OUT   mdout
ln -s $TMPDIR/mdinfo mdinfo
ln -s $TMPDIR/$RST   restrt


cd $TMPDIR

touch $WORKDIR/../running.5_PRODUCTION_${JOB}

if [ $USE_GPU == "YES" ] 
then 
   $AMBERHOME/bin/pmemd.cuda -O -i $WORKDIR/$INPUT -p $WORKDIR/$TOPOLOGY -c $WORKDIR/$CRD -r $RST -x $MDCRD -e $MDEN -o $OUT -inf mdinfo
else
  $MPI_HOME/bin/mpirun -np $NPROCS $AMBERHOME/bin/pmemd.MPI -O -i $WORKDIR/$INPUT -p $WORKDIR/$TOPOLOGY -c $WORKDIR/$CRD -r $RST -x $MDCRD -e $MDEN -o $OUT -inf mdinfo
fi
mv gamd.log $WORKDIR/$GAMDLOG

# Textracting solute coordinates
$AMBERHOME/bin/cpptraj.OMP $WORKDIR/$TOPOLOGY <<EOF
trajin ${MDCRD}
autoimage 
strip :WAT,Cl-,Na+
trajout ${MD}_${txt}_solute.mdcrd netcdf nobox
go
EOF

# Sieving operation  .....
if [ $SIEVE -gt 1  ]
then 
  NSNAP='last'
  $AMBERHOME/bin/cpptraj.OMP $WORKDIR/$TOPOLOGY <<EOF
trajin ${MDCRD} 1  ${NSNAP}  ${SIEVE}
autoimage 
trajout temp.mdcrd netcdf 
go
EOF
  mv temp.mdcrd ${MDCRD}
fi

# Solvent peeling 
if [ $PEEL -gt 0 ]
then 
   $AMBERHOME/bin/cpptraj.OMP $WORKDIR/$TOPOLOGY <<EOF
trajin ${MDCRD}
trajout ${MD}_${txt}_solutewat.mdcrd netcdf nobox
autoimage origin anchor $SOLUTE_MASK
closest $MAXWAT $SOLUTE_MASK solventmask :WAT  oxygen
go
EOF
   mv ${MD}_${txt}_solutewat.mdcrd  $WORKDIR/
fi

# Moving data files
mv $RST $MDCRD $OUT ${MD}_${txt}_solute.mdcrd $WORKDIR/

rm -f $WORKDIR/../running.5_PRODUCTION_${JOB}

cd $WORKDIR/

done 

rm -r -f $TMPDIR 
rm -f restrt mdinfo mdout


