#!/bin/bash

# Statistical analysis of base pair contacts (bp) and
# non-pairing (np) interactions (either base stacking or other ) 

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi

source $APTAMD/ENV/aptamd_env.sh

if [ ! -e options.txt ]
then
   echo "options.txt file not found"
   exit
fi

NRES=$(grep NRES options.txt | sed 's/NRES=//')
NFILES=$(ls md*bp_n.dat | wc -l) 

if [ -z $NFILES ]
then 
   echo "data files md_???_bp_n.dat, md_???_np_n.dat, .... not found!"
   exit
fi

ls md*bp_n.dat | sed 's/_bp_n.dat//' > filelist.txt 

# Optionally onty a fraction of data is used
if [ -z "$PERCEN" ]
then
   echo "Processing the whole data set"
   PERCEN="0"
else
   echo "Processing PERCEN=$PERCEN % of the data set"
fi
if [ "$PERCEN" -gt "100" ]
then
  echo 'PERCEN greater than 100' 
  echo 'Using all data !'
  PERCEN=0
fi
if [ "$PERCEN" -lt "-100" ]
then
  echo 'PERCEN greater than 100' 
  echo 'Using all data !'
  PERCEN=0
fi

if [ $PERCEN -ne 0 ] 
then 

echo "Using only $PERCEN % of the available data"
echo "  > 0  --> From the begining"
echo "  < 0  --> From the end"
if [ "$PERCEN" -gt 0 ]
then
    let  "NFILES_use=( $NFILES *  $PERCEN ) / 100 "
    head -${NFILES_use}  filelist.txt > tmp; mv tmp filelist.txt
    echo "Using only the first $NFILES_use trajectory files"
else
    let  "NFILES_use=($NFILES * ( - $PERCEN ) ) / 100  "
    tail -${NFILES_use}  filelist.txt > tmp; mv tmp filelist.txt
    echo "Using only the last $NFILES_use trajectory files"
fi

fi

REFPDB=$( grep REFPDB options.txt | sed 's/REFPDB=//' )
if [ ! -e  $REFPDB ]
then
   echo "$REFPDB not found"
   exit
fi 
grep 'ATOM  ' $REFPDB | awk '{printf("%s_%i\n", $4,$5)}' | uniq > reslabel.txt 

NRES_check=$(cat reslabel.txt | wc -l)
if [ $NRES -ne $NRES_check ]
then
   echo "NRES=$NRES in options.txt, but NRES=$NRES_check in $REFPDB"
   exit
fi

$OCTAVE -q $APTAMD/STRUCT/stat_bp_np.m 


