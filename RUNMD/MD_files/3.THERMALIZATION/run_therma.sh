#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
source $APTAMD/ENV/aptamd_env.sh

TOPOLOGY=$1

mpirun=$MPI_HOME/bin/mpirun
if [ -z "$NPROCS" ]
then
   NPROCS=$(cat /proc/cpuinfo | grep -c processor)
   echo "Using $NPROCS available processors"
else
   echo "Using $NPROCS processors as predefined "
fi

$mpirun -np $NPROCS $AMBERHOME/bin/pmemd.MPI -O -i therma_50.inp -o therma_50.out -p $TOPOLOGY -c min_all.rst -r therma_50.rst -e therma_50.mden -x therma_50.mdcrd

gzip therma_50.mden therma_50.mdcrd

$mpirun -np $NPROCS $AMBERHOME/bin/pmemd.MPI -O -i therma_100.inp -o therma_100.out -p $TOPOLOGY -c therma_50.rst -r therma_100.rst -e therma_100.mden -x therma_100.mdcrd

gzip therma_100.mden therma_100.mdcrd

$mpirun -np $NPROCS $AMBERHOME/bin/pmemd.MPI -O -i therma_150.inp -o therma_150.out -p $TOPOLOGY -c therma_100.rst -r therma_150.rst -e therma_150.mden -x therma_150.mdcrd

gzip therma_150.mden therma_150.mdcrd

$mpirun -np $NPROCS $AMBERHOME/bin/pmemd.MPI -O -i therma_200.inp -o therma_200.out -p $TOPOLOGY -c therma_150.rst -r therma_200.rst -e therma_200.mden -x therma_200.mdcrd


gzip therma_200.mden therma_200.mdcrd

$mpirun -np $NPROCS $AMBERHOME/bin/pmemd.MPI -O -i therma_250.inp -o therma_250.out -p $TOPOLOGY -c therma_200.rst -r therma_250.rst -e therma_250.mden -x therma_250.mdcrd

gzip therma_250.mden therma_250.mdcrd

$mpirun -np $NPROCS $AMBERHOME/bin/pmemd.MPI -O -i therma_300.inp -o therma_300.out -p $TOPOLOGY -c therma_250.rst -r therma_300.rst -e therma_300.mden -x therma_300.mdcrd

gzip therma_300.mden therma_300.mdcrd

