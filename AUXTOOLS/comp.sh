#!/bin/bash

# gfortran compilation of auxiliary FORTRAN codes

DEBUG="0"
EXTRA=""

ARG=$1

if [ -z $ARG ]
then
   echo "Usage: comp.sh [ all | nmode | program1.f program2.f ...]"
   exit
fi

if [ $ARG == "all" ]
then
    list_f=$(ls *.f)
elif [ $ARG == "nmode" ]
then
     list_f=""
else
    list_f="$*"
fi


for file_f in  $list_f
do


code=${file_f%%.*}
file_exe=${code}

echo "Compiling $code ..."

if [ $code == "fullcontact" ] 
then 
    EXTRA0=$EXTRA
    EXTRA="${EXTRA0} -fconvert=big-endian -I/usr/include/ -L/usr/lib64/ -lnetcdf -lnetcdff"
fi
if [ $DEBUG -eq 1 ]
then 
   gfortran $file_f -fbacktrace -finit-real=nan -g -pg -fcheck=all -fdump-core  \
   -fmax-errors=1 -ffpe-trap=invalid -Wconversion-extra  ${EXTRA}  -o $file_exe
else
   gfortran $file_f ${EXTRA}  -O2 -o $file_exe
fi
if [ $code == "fullcontact" ] ; then EXTRA=$EXTRA0; fi 

done

if [ $ARG == "all" ] ||  [ $ARG == "nmode" ]
then 

if [ -e nmode_standalone ]
then

cd nmode_standalone
make
mv nmode ../
mv ../

fi

fi

