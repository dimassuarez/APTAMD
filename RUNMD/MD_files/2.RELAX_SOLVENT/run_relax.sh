#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
source $APTAMD/ENV/geuo_env.sh

TOPOLOGY=$1

mpirun=$MPI_HOME/bin/mpirun 
if [ -z "$NPROCS" ]
then
   NPROCS=$(cat /proc/cpuinfo | grep -c processor)
   echo "Using $NPROCS available processors"
else
   echo "Using $NPROCS processors as predefined "
fi

$mpirun -np $NPROCS $AMBERHOME/bin/pmemd.MPI -O -i min_solv.inp -o min_solv.out -p $TOPOLOGY -c ../initial.crd -r min_solv.rst

$mpirun -np $NPROCS $AMBERHOME/bin/pmemd.MPI -O -i md_solv.inp -o md_solv.out -p $TOPOLOGY -c min_solv.rst -r md_solv.rst -e md_solv.mden -x md_solv.mdcrd

$mpirun -np $NPROCS $AMBERHOME/bin/pmemd.MPI -O -i min_solv_postmd.inp -o min_solv_postmd.out -p $TOPOLOGY -c md_solv.rst -r min_solv_postmd.rst

$mpirun -np $NPROCS $AMBERHOME/bin/pmemd.MPI -O -i min_all.inp -o min_all.out -p $TOPOLOGY -c min_solv_postmd.rst -r min_all.rst


