#!/bin/bash

# Only 3 fragments: 

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


# Checking if this is a 3-calc or 1-calc ENTRO  system
# Only the first filename in the list is checked

file=$(head -1 $1) 
prefix=${file%%.*}
host_A=${prefix}.nmode_host.gz
lig_A=${prefix}.nmode_lig.gz
host_B=${prefix}.nmode4.gz
lig_B=${prefix}.nmode5.gz

if [ -e "$host_A" ] && [ -e "$lig_A"  ]
then  
    echo "==================================================" 
    echo "Processing 3 calcs: prot, host, lig." 
    echo "==================================================" 
    fragments="1"
    mode_frag="A"
elif [ -e "$host_B" ] && [ -e "$lig_B"  ]
then  
    echo "==================================================" 
    echo "Processing 3 calcs: prot, host, lig." 
    echo "==================================================" 
    fragments="1"
    mode_frag="B"
else 
    echo "==================================================" 
    echo "Processing 1 calcs: prot." 
    echo "==================================================" 
    fragments="0"
fi

# Output files 
rm -f ENTRO.dat ENTRO_DIFF.dat 

if [ "$fragments" -eq "1" ]
then 
  echo '# PREFIX  ENTRO(CMPLX) ENTRO_ROT(CMPLX) ENTRO_VIB(CMPLX) ENTRO(HOST) ENTRO_ROT(HOST) ENTRO_VIB(HOST) ENTRO(LIG)  ENTRO_ROT(LIG) ENTRO_VIB(LIG)' > ENTRO.dat 
  echo '# PREFIX  diff_ENTRO  -T*diff_ENTRO diff_ENTRO_ROT  -T*diff_ENTRO_ROT  diff_ENTRO_VIB -T*diff_ENTRO_VIB'  > ENTRO_DIFF.dat 
else
  echo '# PREFIX  ENTRO(CMPLX) ENTRO_ROT(CMPLX) ENTRO_VIB(CMPLX) ' > ENTRO.dat 
fi

for file in $(cat $1)
do

# Filenames
  prefix=${file%%.*}
  prot=${prefix}.nmode
 
  printf '\r%s:' "$prefix" 

# Get ENTRO
  GET_ENTRO  $prot
  ENTRO_prot=$ENTRO
  ENTRO_TRANS_prot=$ENTRO_TRANS
  ENTRO_ROT_prot=$ENTRO_ROT
  ENTRO_VIB_prot=$ENTRO_VIB

if [ "$fragments" -eq "1" ]
then 
  if [ $mode_frag == "A" ]
  then
      host=${prefix}.nmode_host
      lig=${prefix}.nmode_lig
  else
      host=${prefix}.nmode4
      lig=${prefix}.nmode5
  fi
  GET_ENTRO $host
  ENTRO_host=$ENTRO
  ENTRO_TRANS_host=$ENTRO_TRANS
  ENTRO_ROT_host=$ENTRO_ROT
  ENTRO_VIB_host=$ENTRO_VIB
  GET_ENTRO  $lig
  ENTRO_lig=$ENTRO
  ENTRO_TRANS_lig=$ENTRO_TRANS
  ENTRO_ROT_lig=$ENTRO_ROT
  ENTRO_VIB_lig=$ENTRO_VIB

# Computing Differences
 if [ $ENTRO_prot != "0.000" ] && [ $ENTRO_host != "0.000" ] && [ $ENTRO_lig != "0.000" ]
 then 
   ENTRO_diff=$(echo "$ENTRO_prot $ENTRO_host $ENTRO_lig " | awk '{s= $1 - $2 - $3 ;printf("%12.4f",s)}')
   ENTRO_TRANS_diff=$(echo "$ENTRO_TRANS_prot $ENTRO_TRANS_host $ENTRO_TRANS_lig " | awk '{s= $1 - $2 - $3 ;printf("%12.4f",s)}')
   ENTRO_ROT_diff=$(echo "$ENTRO_ROT_prot $ENTRO_ROT_host $ENTRO_ROT_lig " | awk '{s= $1 - $2 - $3 ;printf("%12.4f",s)}')
   ENTRO_VIB_diff=$(echo "$ENTRO_VIB_prot $ENTRO_VIB_host $ENTRO_VIB_lig " | awk '{s= $1 - $2 - $3 ;printf("%12.4f",s)}')
   TEMP_ENTRO_diff=$(echo "$ENTRO_diff $TEMP " | awk '{g= - $1 * $2  ;printf("%12.4f",g)}')
   TEMP_ENTRO_TRANS_diff=$(echo "$ENTRO_TRANS_diff $TEMP " | awk '{g= - $1 * $2  ;printf("%12.4f",g)}')
   TEMP_ENTRO_ROT_diff=$(echo "$ENTRO_ROT_diff $TEMP " | awk '{g= - $1 * $2  ;printf("%12.4f",g)}')
   TEMP_ENTRO_VIB_diff=$(echo "$ENTRO_VIB_diff $TEMP " | awk '{g= - $1 * $2  ;printf("%12.4f",g)}')
  else
    ENTRO_diff="0.000"
    ENTRO_ROT_diff="0.000"
    ENTRO_VIB_diff="0.000"
    TEMP_ENTRO_diff="0.000"
    TEMP_ENTRO_ROT_diff="0.000"
    TEMP_ENTRO_VIB_diff="0.000"
  fi

fi
# Printing out data
if [ "$fragments" -eq "1" ]
then 
    echo  "$prefix  $ENTRO_prot  $ENTRO_ROT_prot  $ENTRO_VIB_prot  $ENTRO_host  $ENTRO_ROT_host  $ENTRO_VIB_host    $ENTRO_lig  $ENTRO_ROT_lig  $ENTRO_VIB_lig  "  >>  ENTRO.dat
    echo  "$prefix  $ENTRO_diff  $TEMP_ENTRO_diff    $ENTRO_ROT_diff  $TEMP_ENTRO_ROT_diff  $ENTRO_VIB_diff  $TEMP_ENTRO_VIB_diff  "  >>  ENTRO_DIFF.dat

else 
    echo  "$prefix  $ENTRO_prot  $ENTRO_ROT_prot  $ENTRO_VIB_prot  "  >>  ENTRO.dat
fi

unset ENTRO ENTRO_ROT ENTRO_VIB

done

if [ "$fragments" -eq "1" ]
then 
  echo "# TRANSLATIONAL ENTROPIES (CMPLX_HOST_LIG) $ENTRO_TRANS_prot $ENTRO_TRANS_host $ENTRO_TRANS_lig  " >> ENTRO.dat
  echo "# TRANSLATIONAL ENTROPY DIFF  -T*DIFF  $ENTRO_TRANS_diff $TEMP_ENTRO_TRANS_diff  " >> ENTRO_DIFF.dat
else 
  echo "# TRANSLATIONAL ENTROPY  $ENTRO_TRANS_prot " >> ENTRO.dat
fi

printf '\n' 

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
   if [ "$fragments" -eq "1" ]
   then
     tar rvf OUTPUT.tar *.nmode_host.gz
     tar rvf OUTPUT.tar *.nmode_lig.gz
     for ((i=1;i<=20;i++)); do tar rvf OUTPUT.tar *.nmode${i}.gz; done
   fi
   rm -f *.nmode.gz
   rm -f *.out_min.gz
   rm -f *.pdb_min.gz
   if [ "$fragments" -eq "1" ]
   then
     rm -f *.nmode_host.gz
     rm -f *.nmode_lig.gz
     for ((i=1;i<=20;i++)); do rm -f  *.nmode${i}.gz; done
   fi
   echo "All output files are archived into OUTPUT.tar" 
fi 


  
