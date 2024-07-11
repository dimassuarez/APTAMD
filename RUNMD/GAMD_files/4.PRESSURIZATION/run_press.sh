#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
source $APTAMD/ENV/geuo_env.sh

TOPOLOGY=$1
if [ -z $USE_GPU ]; then USE_GPU="YES" ; fi

mpirun=$MPI_HOME/bin/mpirun
if [ -z "$NPROCS" ]
then
   NPROCS=$(cat /proc/cpuinfo | grep -c processor)
   echo "Using $NPROCS available processors"
else
   echo "Using $NPROCS processors as predefined "
fi

NVIDIA=$(whereis -b nvidia-smi | awk '{print NF}') 

if [ "$NVIDIA" -gt "1" ]  &&  [ "$USE_GPU" == "YES" ] 
then
     $AMBERHOME/bin/pmemd.cuda -O -i md_npt_mc.inp -c therma_300.rst -o press.out -p $TOPOLOGY -x press.mdcrd -r press.rst 
else
     $mpirun -np $NPROCS  $AMBERHOME/bin/pmemd.MPI -O -i md_npt_mc.inp -c therma_300.rst -o press.out -p $TOPOLOGY -x press.mdcrd -r press.rst 
fi

