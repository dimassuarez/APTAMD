#!/bin/bash
#
# This script extracts DATA from nmode output files
# to evaluate RRHO entropy terms for MMPB(GB)SA calculations
#
# Usage:
#
# bash nmode_data_parser.sh LIST
#
# where LIST is a file containing the list of snapshot IDs. 
# For example,
#
#    snap_0001.pdb
#    snap_0002.pdb
#    snap_0003.pdb
#    snap_0003.pdb
#    ....
#
#  On output the script prints out several .dat files
#  that should be further processed using the stat_plot.sh script
# 


#
TEMP=0.300

function GET_ENTRO()
{

#GRMS threshold value check. IF GRMS > THRESHOLD, then all entropy values are zeroed.
THRESHOLD=0.00001

# Remember that by default variables in a BASH function are NOT local variables, but global.

  file=$1

  declare -a ENTRO_DAT=""
  ENTRO_DAT=( $(zcat $file | grep -A 4 " Total                " | awk '{print $4}' ) )

  ENTRO=${ENTRO_DAT[0]}
  ENTRO_TRANS=${ENTRO_DAT[1]}
  ENTRO_ROT=${ENTRO_DAT[2]}
  ENTRO_VIB=${ENTRO_DAT[3]}
  GRMS=$(zcat $file | grep "GRMS" |tail -1   |  awk '{print $9}' )

  CHECK=$(echo "$THRESHOLD $GRMS " | awk '{ if ($2 > 0.0) {c= int ( $1 / $2);} else {c=1;} ;printf("%d",c)} ')

  if [ "$CHECK" == "0" ]
  then
     printf '\n%s ' " Gradient too large in $file GRMS= $GRMS"
     printf '\n '
     ENTRO="0.000"
     ENTRO_ROT="0.000"
     ENTRO_VIB="0.000"
     return
  fi

  if [ "$ENTRO" == "" ]
  then
        ENTRO="0.000"
  fi
  if [ "$ENTRO_ROT" == "" ]
  then
        ENTRO_ROT="0.000"
  fi
  if [ "$ENTRO_VIB" == "" ]
  then
        ENTRO_VIB="0.000"
  fi
}

# Checking if this is a 5-calc or 3-calc MM-PB system
# Only the first filename in the list is checked

file=$(head -1 $1) 
prefix=${file%%.*}
cmplx=${prefix}.nmode.gz

if [ ! -e "$cmplx" ]
then 
   echo "$cmplx , which is the first file in $1, does not exist !"
   echo "This file is needed to setup data collection"
   echo "Fix $1 and rerun the script"
   exit
fi   

NFRAG=0
for ((i=4;i<=20;i++))
do
   if [ -e ${prefix}.nmode${i}.gz ]
   then
      let "NFRAG=$NFRAG+1"
   fi
done
if [ $NFRAG -gt 1 ]
then 
  echo "===============================================================================" 
  echo "Processing $NFRAG+1 calcs: cmplx, frag1, frag2, ... frag${NFRAG}"
  echo "===============================================================================" 
  fragments="1"
else
  fragments="0"
fi

# Output files 
rm -f ENTRO.dat ENTRO_DIFF.dat

if [ "$fragments" -eq "1" ]
then 
  SFRAGTXT=""
  for ((ifrag=1;ifrag<=NFRAG;ifrag++))
  do
     SFRAGTXT="$SFRAGTXT ENTRO(FRAG${ifrag}) (ENTRO_TRANS${ifrag})  (ENTRO_ROT${ifrag}) ENTRO_VIB(FRAG${ifrag})"
  done
  echo "# PREFIX  ENTRO(CMPLX) ENTRO_TRANS(CMPLX)  ENTRO_ROT(CMPLX) ENTRO_VIB(CMPLX) $SFRAGTXT " > ENTRO_COMP.dat 
  echo '# PREFIX  diff_ENTRO  diff_ENTRO_TRANS diff_ENRO_ROT  diff_ENTRO_VIB '  > ENTRO_COMP_DIFF.dat 
else
  echo "# PREFIX  ENTRO(CMPLX) ENTRO_TRANS(CMPLX)  ENTRO_ROT(CMPLX) ENTRO_VIB(CMPLX)  " > ENTRO_COMP.dat 
fi

for file in $(cat $1)
do

# Filenames
  prefix=${file%%.*}
  cmplx=${prefix}.nmode.gz 
  printf '\r%s:' "$prefix" 

# Get ENTRO data for CMPLX 
  GET_ENTRO  $cmplx
  ENTRO_cmplx=$ENTRO
  ENTRO_TRANS_cmplx=$ENTRO_TRANS
  ENTRO_ROT_cmplx=$ENTRO_ROT
  ENTRO_VIB_cmplx=$ENTRO_VIB

if [ "$fragments" -eq "1" ]
then 
   declare -a ENTRO_frag
   declare -a ENTRO_ROT_frag
   declare -a ENTRO_VIB_frag
# Get ENTRO data for FRAG
  for ((ifrag=1;ifrag<=NFRAG;ifrag++))
  do
    let " iout = $ifrag + 3 "
    frag=${prefix}.nmode${iout} 
    GET_ENTRO $frag
    ENTRO_frag["$ifrag"]=$ENTRO
    ENTRO_TRANS_frag["$ifrag"]=$ENTRO_TRANS
    ENTRO_ROT_frag["$ifrag"]=$ENTRO_ROT
    ENTRO_VIB_frag["$ifrag"]=$ENTRO_VIB
  done 

  SUM_ENTRO="0.0"
  SUM_ENTRO_TRANS="0.0"
  SUM_ENTRO_ROT="0.0"
  SUM_ENTRO_VIB="0.0"
  for ((ifrag=1;ifrag<=NFRAG;ifrag++))
  do 
  TER_ENTRO=${ENTRO_frag["$ifrag"]}
  TER_ENTRO_TRANS=${ENTRO_TRANS_frag["$ifrag"]}
  TER_ENTRO_ROT=${ENTRO_ROT_frag["$ifrag"]}
  TER_ENTRO_VIB=${ENTRO_VIB_frag["$ifrag"]}
 
  SUM_ENTRO=$(echo "$SUM_ENTRO $TER_ENTRO " | awk '{s= $1 + $2  ;printf("%12.4f",s)}')
  SUM_ENTRO_TRANS=$(echo "$SUM_ENTRO_TRANS $TER_ENTRO_TRANS " | awk '{s= $1 + $2  ;printf("%12.4f",s)}')
  SUM_ENTRO_ROT=$(echo "$SUM_ENTRO_ROT $TER_ENTRO_ROT " | awk '{s= $1 + $2  ;printf("%12.4f",s)}')
  SUM_ENTRO_VIB=$(echo "$SUM_ENTRO_VIB $TER_ENTRO_VIB " | awk '{s= $1 + $2  ;printf("%12.4f",s)}')
  done 

# Computing Differences

  ENTRO_diff=$(echo "$ENTRO_cmplx  $SUM_ENTRO" | awk '{s= $1 - $2 ;printf("%12.4f",s)}') 
  ENTRO_TRANS_diff=$(echo "$ENTRO_TRANS_cmplx  $SUM_ENTRO_TRANS" | awk '{s= $1 - $2 ;printf("%12.4f",s)}') 
  ENTRO_ROT_diff=$(echo "$ENTRO_ROT_cmplx  $SUM_ENTRO_ROT  " | awk '{s= $1 - $2 ;printf("%12.4f",s)}') 
  ENTRO_VIB_diff=$(echo "$ENTRO_VIB_cmplx  $SUM_ENTRO_VIB  " | awk '{s= $1 - $2 ;printf("%12.4f",s)}') 

# Printing out data
  SFRAGTXT=""
  for ((ifrag=1;ifrag<=NFRAG;ifrag++))
  do
       TER_ENTRO=${ENTRO_frag["$ifrag"]}
       TER_ENTRO_TRANS=${ENTRO_TRANS_frag["$ifrag"]}
       TER_ENTRO_ROT=${ENTRO_ROT_frag["$ifrag"]}
       TER_ENTRO_VIB=${ENTRO_VIB_frag["$ifrag"]}
       SFRAGTXT="$SFRAGTXT $TER_ENTRO $TER_ENTRO_TRANS $TER_ENTRO_ROT $TER_ENTRO_VIB "
  done
  echo "$prefix $ENTRO_cmplx $ENTRO_TRANS_cmplx $ENTRO_ROT_cmplx $ENTRO_VIB_cmplx $SFRAGTXT" >> ENTRO_COMP.dat
  echo "$prefix $ENTRO_diff $ENTRO_TRANS_diff $ENTRO_ROT_diff $ENTRO_VIB_diff " >> ENTRO_COMP_DIFF.dat

else
  echo "$prefix $ENTRO_cmplx $ENTRO_TRANS_cmplx $ENTRO_ROT_cmplx $ENTRO_VIB_cmplx " >> ENTRO_COMP.dat
fi

printf '\n'

done

if [ ! "$DO_ENTRO_NMODE_TAR" ]
then
  echo " Shall we tar all output files in OUTPUT.tar?  (Y/N) ?"
  read YES
  
  if [ "$YES" == "y" ] || [ "$YES" == "Y" ]
  then 
    DO_ENTRO_NMODE_TAR="1"
  else
    DO_ENTRO_NMODE_TAR="0"
  fi

fi

if [ "$DO_ENTRO_NMODE_TAR" -eq "1" ]
then
   rm -f OUTPUT.tar
   tar cvf OUTPUT.tar *.nmode.gz
   tar rvf OUTPUT.tar *.out_min.gz
   tar rvf OUTPUT.tar *.pdb_min.gz
   if [ $fragments -eq 1 ] 
   then 
       for ((i=1;i<=NFRAG;i++)); do let "j=$i+3"; tar rvf OUTPUT.tar *.nmode${j}.gz; done
   fi
   rm -f *.nmode.gz
   rm -f *.out_min.gz
   rm -f *.pdb_min.gz
   if [ "$fragments" -eq "1" ]
   then
     rm -f *.nmode?.gz rm -f *.nmode??.gz
   fi
   echo "All output files are archived into OUTPUT.tar" 
fi 


  
