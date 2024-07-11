#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
source $APTAMD/ENV/geuo_env.sh

#Local directory
WORKDIR=$PWD

# TOPOLOGY pathfile which is the script argument should be relative to WORKDIR
TOPOLOGY=$1

if [ -z $TOPOLOGY ]; then  echo "TOPOLOGY not defined"; exit; fi
if [ -z $USE_GPU ]; then USE_GPU='YES' ; fi

# Number of atoms in topology
NAT=$(grep -A3 'FLAG POINTERS' $TOPOLOGY | grep -v '%' | head -1 | awk '{print $1}')
# GAMD settings are selected on the basis of the NAT value
let "NTAVE=4*$NAT"
NTCMDPREP="0"
let "NTCMD=14*$NTAVE"
let "NTEBPREP=$NTCMD"
let "NTEB=10*$NTCMD"
let "NTSLIM=$NTAVE+$NTCMD+$NTEBPREP+$NTEB+$NTCMDPREP"

cp md_npt_gamd_eq_DUMMY.inp md_npt_gamd_eq.inp

sed -i "s/DUMMY_NAT/${NAT}/" md_npt_gamd_eq.inp
sed -i "s/DUMMY_NTAVE/${NTAVE}/" md_npt_gamd_eq.inp
sed -i "s/DUMMY_NTCMDPREP/${NTCMDPREP}/" md_npt_gamd_eq.inp
sed -i "s/DUMMY_NTCMD/${NTCMD}/" md_npt_gamd_eq.inp
sed -i "s/DUMMY_NTEBPREP/${NTEBPREP}/" md_npt_gamd_eq.inp
sed -i "s/DUMMY_NTEB/${NTEB}/" md_npt_gamd_eq.inp
sed -i "s/DUMMY_NTSLIM/${NTSLIM}/" md_npt_gamd_eq.inp

if [ $USE_GPU == "YES" ]
then
   $AMBERHOME/bin/pmemd.cuda -O -i md_npt_gamd_eq.inp -c gamd.rst -o gamd_eq.out -p $TOPOLOGY -x gamd_eq.mdcrd -r gamd_eq.rst -amd gamd_eq.log 
else
  if [ -z "$NPROCS" ]
  then
     NPROCS=$(cat /proc/cpuinfo | grep -c processor)
     echo "Using $NPROCS available processors"
  else
     echo "Using $NPROCS processors as predefined "
  fi
  $MPI_HOME/bin/mpirun -np $NPROCS $AMBERHOME/bin/pmemd.MPI -O -i md_npt_gamd_eq.inp -c gamd.rst -o gamd_eq.out -p $TOPOLOGY -x gamd_eq.mdcrd -r gamd_eq.rst -amd gamd_eq.log 
fi

mv gamd_eq.rst md_000.rst 

let "NTAVE_PROD=2*NTAVE"

cp md_npt_gamd_restart_DUMMY.inp md_npt_gamd_restart.inp
sed -i "s/DUMMY_NTAVE/${NTAVE_PROD}/"  md_npt_gamd_restart.inp

