#!/bin/bash
#
# This script extracts DATA from sander output files
# to evaluate the MM-PB(GB) or  QMMM-PB(GB) calculations. 
# Two versions of MM-PB(GB) are considered:
#
#    1) (QM)MM-PB + Esolute-solvent(VdW) + gamma*A_molsurf 
#    2) (QM)MM-PB(GB)SA 
#
# For option 2, we also estimate the host-ligand interaction energies
# if the fragment data (i.e., host & ligand files) are available
#
# Both linear and non-linear PB calculations can be handled.
# For (QM)MMPBSA, the non-polar term is expected to be ECAVITY+EDISP
# as estimated by option inp=2 in PBSA/AMBER14
#
# Usage:
#
# bash do_energy_analyses.sh LIST
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
#  that should be further processed using the do_stat.sh script
# 

function GET_MMPBSA()
{
 
# Remember that by default variables in a BASH function are NOT local variables, but global.

  TEMPERATURE=0.300  
  RCTE=1.9872
  
  file=$1 

  declare -a DAT=""

  NDAT="0"
  BAD="0.000"
  VDWELEC="0.000"   
  GB="0.000"   
  VDWELEC_14="0.000"  
  PB="0.000"   
  PB1="0.000"   
  PB2="0.000"   
  NP="0.000"  
  ECAV="0.000"  
  EDIS="0.000"  
  ESCF="0.000"
  CMAP="0.000"
  AREA_MOLSURF="0.000"
  ASES_NUM="0.000"
  VSES_NUM="0.000"
  GCAV_CSPT="0.000"

  if [ ! -e ${file}.gz ] 
  then 
      echo "${file}.gz does not exist !"
      return
  fi

  DAT=( $(zcat $file  | grep -v scale |  sed -n '/FINAL/,+15p' |   \
  awk -v "pdie=$PDIE"  '{if ($1 == "BOND")    { bad=$3+$6+$9;  printf(" %12.4f",bad)}  }   \
       {if ($1 == "VDWAALS") { vdwelec=$3+$6/pdie; gb=$9;  printf(" %12.4f %12.4f",vdwelec,gb)}  }   \
       {if ($1 == "1-4" && $2 == "VDW") { vdwelec_14=$4+$8/pdie;  printf(" %12.4f  \n",vdwelec_14)}  }   \
       {if ($1 == "EELEC") { pb=$6;  printf(" %12.4f \n",pb)}  }   \
       {if ($1 == "ECAVITY=") { ecav=$2; edis=$5; np=$2+$5;  printf(" %12.4f %12.4f %12.4f \n",np,ecav,edis)}  }  \
       {if ($1 == "CMAP") { cmap=$3;  printf(" %12.4f \n",cmap) } } '|  \
       uniq   ) )
  NDAT=${#DAT[*]}
  if [ "$NDAT" -ne "9" ] && [ "$NDAT" -ne "10" ] 
  then 
     # NDAT=9 ( no CMAP) NDAT=10 (CMAP present)
     echo " Not enough data in $file NDAT= ${NDAT}" 
     echo ${DAT[*]} 
     echo 
     return
  fi 
  BAD=${DAT[0]}
  VDWELEC=${DAT[1]}
  GB=${DAT[2]}
  VDWELEC_14=${DAT[3]}
  if [ $NDAT -eq 9 ]
  then
    PB=${DAT[5]}
    NP=${DAT[6]}
    ECAV=${DAT[7]}
    EDIS=${DAT[8]}
  else
    CMAP=${DAT[4]}
    PB=${DAT[6]}
    NP=${DAT[7]}
    ECAV=${DAT[8]}
    EDIS=${DAT[9]}
    BAD=$(echo $BAD $CMAP | awk '{print $1+$2}') 
  fi

#
# Getting Surface Data
#
  unset DAT
  declare -a DAT=""
  DAT=($(zcat $file | grep -A6 AREA_MOLSURF | awk '{print $2}'))
  NDAT=${#DAT[@]} 
  if [ $NDAT -eq 7 ]
  then
      AREA_MOLSURF=${DAT[0]}
      ASES_NUM=${DAT[1]}
      VSES_NUM=${DAT[3]}
      GCAV_CSPT=${DAT[6]}
  else
      AREA_MOLSURF=${DAT[0]}
  fi
#
#   Getting SCF energy / Check is D3H4 is to be included
#  
  if [ "$NQMMM" -gt "0" ]
  then
      ESCF=$(zcat -f $file |  sed 's/ESCF =/ESCF=/'  | grep 'ESCF=' | head -1 |awk '{scf=$2;printf(" %12.4f \n",scf)}')
      D3H4=$(zcat -f $file |   grep 'D3H4 ' | tail -1 |awk '{printf(" %12.4f \n",$NF)}')
      if [ -z $D3H4 ]; then D3H4="0.000" ;fi
      ESCF_D3H4=$(echo $ESCF  $D3H4 | awk '{printf(" %12.4f \n",$1+$2)}')
      ESCF=$ESCF_D3H4
  else
      ESCF="0.000"
  fi

# EGAS= " $ESCF +  $BAD +  $VDWELEC +  $VDWELEC_14 + 3 * $RCTE * $TEMPERATURE "
  EGAS=$(echo "$ESCF $BAD $VDWELEC $VDWELEC_14 $RCTE $TEMPERATURE" |  \
   awk '{egas= $1 + $2 + $3 + $4 + ( 3.0 * $5 * $6 ) ;printf("%12.4f",egas)}')

  if [ "$EGAS" == "" ] 
  then 
        EGAS="0.000"
  fi
  if [ "$GB" == "" ] 
  then 
        GB="0.000"
  fi
  if [ "$NP" == "" ] 
  then 
        NP="0.000"
  fi
  if [ "$PB" == "" ] 
  then  
        PB="0.000"
  fi
  if [ "$ESCF" == "" ] 
  then  
        ESCF="0.000"
  fi
  if [ "$ECAV" == "" ] 
  then  
        ECAV="0.000"
  fi
  if [ "$EDIS" == "" ] 
  then  
        ECAV="0.000"
  fi
  if [ "$AREA_MOLSURF" == "" ] 
  then  
        AREA_MOLSURF="0.000"
  fi
  if [ "$ASES_NUM" == "" ] 
  then  
        ASES_NUM="0.000"
  fi
  if [ "$VSES_NUM" == "" ] 
  then  
        VSES_NUM="0.000"
  fi
  if [ "$GCAV_CSPT" == "" ] 
  then  
        GCAV_CSPT="0.000"
  fi
}


function GET_MMRISM()
{
 
# Remember that by default variables in a BASH function are NOT local variables, but global.

  TEMPERATURE=0.300  
  RCTE=1.9872
  
  file=$1 

  declare -a DAT=""

  NDAT="0"
  BAD="0.000"
  VDWELEC="0.000"   
  GB="0.000"   
  NP="0.000"   
  VDWELEC_14="0.000"  
  RISM="0.000"   

  if [ ! -e ${file}.gz ] 
  then 
      echo "${file}.gz does not exist !"
      return
  fi

  DAT=( $(zcat $file  | sed '1,/FINAL RESULTS/d' | sed '/TIMINGS/,$d' |  \
  awk '{if ($1 == "BOND")    { bad=$3+$6+$9;  printf(" %12.4f",bad)}  }   \
       {if ($1 == "VDWAALS") { vdwelec=$3+$6; gb=$9;  printf(" %12.4f %12.4f",vdwelec,gb)}  }   \
       {if ($1 == "1-4" && $2 == "VDW") { vdwelec_14=$4+$8;  printf(" %12.4f  \n",vdwelec_14)}  }   \
       {if ($1 == "ESURF") { np=$3;  printf(" %12.4f \n",np)}  } ' | \
       uniq   ) )
# echo ${DAT[*]} 
  NDAT=${#DAT[*]}
  if [ "$NDAT" -ne "5" ]
  then 
     echo " Not enough data in $file NDAT= ${NDAT}" 
     echo ${DAT[*]} 
     return
  fi 
  BAD=${DAT[0]}
  VDWELEC=${DAT[1]}
  GB=${DAT[2]}
  VDWELEC_14=${DAT[3]}
  NP=${DAT[4]}

  # RISM=$(zcat $file | grep 'ERISM' | head -1 | awk '{print $9}')
  RISM=$(zcat $file | grep 'rism_excessChemicalPotentialPCPLUS' | head -1 | awk '{printf("%12.4f  \n",$2)}')

  EGAS=$(echo "$ESCF $BAD $VDWELEC $VDWELEC_14 $RCTE $TEMPERATURE" |  \
   awk '{egas= $1 + $2 + $3 + $4 + ( 3.0 * $5 * $6 ) ;printf("%12.4f",egas)}')

  if [ "$EGAS" == "" ] 
  then 
        EGAS="0.000"
  fi
  if [ "$GB" == "" ] 
  then 
        GB="0.000"
  fi
  if [ "$NP" == "" ] 
  then 
        NP="0.000"
  fi
  if [ "$RISM" == "" ] 
  then  
        RISM="0.000"
  fi
}


# Checking if this is a 5-calc or 3-calc MM-PB system
# Only the first filename in the list is checked

file=$(head -1 $1) 
prefix=${file%%.*}
cmplx=${prefix}.out3.gz


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
   if [ -e ${prefix}.out${i}.gz ]
   then
      let "NFRAG=$NFRAG+1"
   fi
done
if [ $NFRAG -gt 1 ]
then 
  echo "===============================================================================" 
  echo "Processing $NFRAG+3 calcs: cmplxwat, wat, cmplx, frag1, frag2, ... frag${NFRAG}"
  echo "===============================================================================" 
  fragments="1"
else
  fragments="0"
fi

# Checking if this is RISM
NRISM=$(zcat ${cmplx} | grep -v 'memory allocation summary' |  grep -c '3D-RISM')
if [ "$NRISM" -gt "0" ]
then 

#####################################3 RISM begins here ############################
# RISM is treated in a completely different way

echo "==================================================" 
echo "Processing RISM calculations." 
echo "==================================================" 

# Output files 
rm -f ECOMP.dat ECOMP_DIFF.dat G_MMRISM.dat G_MMGBSA.dat NWAT.dat NONPOLAR.dat  NONPOLAR_DIFF.dat

if [ "$fragments" -eq "1" ]
then 
  EFRAGTXT=""
  GFRAGTXT=""
  for ((ifrag=1;ifrag<=NFRAG;ifrag++))
  do
     EFRAGTXT="$EFRAGTXT EGAS(FRAG${ifrag}) RISM(FRAG${ifrag}) GB(FRAG${ifrag}) NP(FRAG${ifrag})"
     GFRAGTXT="$GFRAGTXT G(FRAG${ifrag})" 
  done
  echo "# PREFIX  EGAS(CMPLX)  RISM(CMPLX) GB(CMPLX)  NP(CMPLX) $EFRAGTXT " > ECOMP.dat 
  echo '# PREFIX  diff_EGAS  diff_RISM diff_GB  diff_SURF  diff_G_RISM  diff_G_GB'  > ECOMP_DIFF.dat 
  echo "# PREFIX  G(CMPLX) $GFRAGTXT diff_G "  > G_MMGBSA.dat 
  echo "# PREFIX  G(CMPLX) $GFRAGTXT diff_G "  > G_MMRISM.dat 
  echo '# PREFIX  NWAT(CMPLX) '  > NWAT.dat 
else
  echo '# PREFIX  EGAS(CMPLX)  RISM(CMPLX) GB(CMPLX)  NP(CMPLX)   ' > ECOMP.dat 
  echo '# PREFIX  G_RISM(CMPLX)  G_GB(CMPLX) '  > G_MMRISMGBSA.dat 
  echo '# PREFIX  NWAT(CMPLX) '  > NWAT.dat 
fi

for file in $(cat $1)
do

# Filenames
  prefix=${file%%.*}
  cmplxwat=${prefix}.out1
  wat=${prefix}.out2
  cmplx=${prefix}.out3
 
  printf '\r%s:' "$prefix" 

# Get Egas, PB, GB and  Surf for CMPLX 
  GET_MMRISM $cmplx
  EGAS_cmplx=$EGAS
  ESCF_cmplx=$ESCF
  GB_cmplx=$GB
  NP_cmplx=$NP
  RISM_cmplx=$RISM

if [ "$fragments" -eq "1" ]
then 
   declare -a EGAS_frag
   declare -a ESCF_frag
   declare -a GB_frag
   declare -a RISM_frag
   declare -a NP_frag

# Get Egas, GB and  Surf for HOST

  for ((ifrag=1;ifrag<=NFRAG;ifrag++))
  do
    let " iout = $ifrag + 3 "
    frag=${prefix}.out${iout} 
    GET_MMRISM $frag
    EGAS_frag["$ifrag"]=$EGAS
    ESCF_frag["$ifrag"]=$ESCF 
    GB_frag["$ifrag"]=$GB
    NP_frag["$ifrag"]=$NP
    RISM_frag["$ifrag"]=$RISM
  done 

fi

# Computing G energies
 if [ "$EGAS_cmplx" != "0.000" ] &&  [ "$RISM_cmplx" != "0.000" ] &&  [ "$NP_cmplx" != "0.000" ] 
 then 
# GCMPLX_RISM= "$EGAS_cmplx + $RISM_cmplx "
  GCMPLX_RISM=$(echo "$EGAS_cmplx $RISM_cmplx " | awk '{g= $1 + $2 ;printf("%12.4f",g)}')
 else
    GCMPLX_RISM="0.000"
 fi
 if [ "$EGAS_cmplx" != "0.000" ] &&  [ "$GB_cmplx" != "0.000" ] &&  [ "$NP_cmplx" != "0.000" ] 
 then 
#   GCMPLX_GB="$EGAS_cmplx + $GB_cmplx + $NP_cmplx"
    GCMPLX_GB=$(echo "$EGAS_cmplx $GB_cmplx $NP_cmplx " | awk '{g= $1 + $2 + $3 ;printf("%12.4f",g)}')
 else
    GCMPLX_GB="0.000"
 fi

if [ "$fragments" -eq "1" ]
then 

 declare -a GFRAG_RISM=""
 declare -a GFRAG_GB=""

 SUM_G_RISM="0.0"
 SUM_G_GB="0.0"
 SUM_RISM="0.0"
 SUM_GB="0.0"
 SUM_EGAS="0.0"
 SUM_ESCF="0.0"
 SUM_NP="0.0"

 for ((ifrag=1;ifrag<=NFRAG;ifrag++))
 do 

 TER_EGAS=${EGAS_frag["$ifrag"]}
 TER_ESCF=${ESCF_frag["$ifrag"]}
 TER_RISM=${RISM_frag["$ifrag"]}
 TER_GB=${GB_frag["$ifrag"]}
 TER_NP=${NP_frag["$ifrag"]}

 if [ ${EGAS_frag["$ifrag"]} != "0.000" ] &&  [ ${RISM_frag["$ifrag"]} != "0.000" ] 
 then 
    GFRAG_RISM["$ifrag"]=$(echo "$TER_EGAS  $TER_RISM " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 else
    GFRAG_RISM["$ifrag"]="0.000"
 fi
 if [ ${EGAS_frag["$ifrag"]} != "0.000" ] &&  [ ${GB_frag["$ifrag"]} != "0.000" ] &&  [ ${NP_frag["$ifrag"]} != "0.000" ] 
 then 
    GFRAG_GB["$ifrag"]=$(echo "$TER_EGAS  $TER_GB $TER_NP " | awk '{g= $1 + $2 + $3 ;printf("%12.4f",g)}')
 else
    GFRAG_GB["$ifrag"]="0.000"
 fi

 TER_G_RISM=${GFRAG_RISM["$ifrag"]}
 TER_G_GB=${GFRAG_GB["$ifrag"]}
 SUM_G_RISM=$(echo "$SUM_G_RISM $TER_G_RISM " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_G_GB=$(echo "$SUM_G_GB $TER_G_GB " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
  

 SUM_RISM=$(echo "$SUM_RISM $TER_RISM " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_GB=$(echo "$SUM_GB $TER_GB " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_EGAS=$(echo "$SUM_EGAS $TER_EGAS " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_ESCF=$(echo "$SUM_ESCF $TER_ESCF " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_NP=$(echo "$SUM_NP $TER_NP " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')

done 

# Computing Differences

 G_RISM_diff=$(echo "$GCMPLX_RISM $SUM_G_RISM" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 G_GB_diff=$(echo "$GCMPLX_GB $SUM_G_GB" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 RISM_diff=$(echo "$RISM_cmplx  $SUM_RISM" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 GB_diff=$(echo "$GB_cmplx  $SUM_GB" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 EGAS_diff=$(echo "$EGAS_cmplx  $SUM_EGAS" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 ESCF_diff=$(echo "$ESCF_cmplx  $SUM_ESCF" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 NP_diff=$(echo "$NP_cmplx  $NP_ESCF" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 


# Printing out data
EFRAGTXT=""
RISMFRAGTXT=""
GBFRAGTXT=""
for ((ifrag=1;ifrag<=NFRAG;ifrag++))
do
       TER_EGAS=${EGAS_frag["$ifrag"]}
       TER_RISM=${RISM_frag["$ifrag"]}
       TER_GB=${GB_frag["$ifrag"]}
       TER_NP=${NP_frag["$ifrag"]}
       EFRAGTXT="$EFRAGTXT $TER_EGAS $TER_RISM $TER_GB $TER_NP "
       TER_G_RISM=${GFRAG_RISM["$ifrag"]}
       TER_G_GB=${GFRAG_GB["$ifrag"]}
       RISMFRAGTXT="$RISMFRAGTXT $TER_G_RISM "
       GBFRAGTXT="$GBFRAGTXT $TER_G_GB "
done
echo "$prefix $EGAS_cmplx $RISM_cmplx $GB_cmplx $NP_cmplx $EFRAGTXT" >> ECOMP.dat
echo "$prefix $EGAS_diff $RISM_diff $GB_diff $NP_diff $G_RISM_diff $G_GB_diff" >> ECOMP_DIFF.dat
echo "$prefix $GCMPLX_RISM  $RISMFRAGTXT $G_RISM_diff" >> G_MMRISM.dat
echo "$prefix $GCMPLX_GB $GBFRAGTXT $G_GB_diff" >> G_MMGBSA.dat
echo "$prefix $NWAT" >> NWAT.dat

else 

  echo "$prefix $EGAS_cmplx $RISM_cmplx $GB_cmplx $NP_cmplx " >> ECOMP.dat
  echo "$prefix $GCMPLX_RISM $GCMPLX_GB " >> G_MMRISMGBSA.dat
  echo "$prefix $NWAT" >> NWAT.dat

fi

done

printf '\n'

############################  RISM ends here ############################################
else 

# WARNING:
# The rest of the SCRIPT deaals with PB(GB) calcs.

# For MMPBSA we can check now different eps_in values
if [ ! "$PDIE" ]
then
    echo "PDIE not defined. Using PDIE=1.0" 
    PDIE="1.0"
else
    echo "PDIE defined. Using PDIE=$PDIE" 
fi

# Checking if PDIE matches intdiel
INTDIEL=$(zcat ${cmplx} | grep  'intdiel =' | head -1 | awk '{printf("%i",$NF)}')
if [ $INTDIEL -eq ${PDIE/.*} ]
then
    echo "PDIE matches INTDIEL = $PDIE, hence changing to PDIE=1" 
    PDIE="1.0"
fi   

# Checking if this is non-linear PB or linear PB
NPB=$(zcat ${cmplx} | grep -c -i 'pbsa')

if [ "$NPB" -eq  "0" ]
then 
    echo "==================================================" 
    echo "Processing only GB calculations." 
    echo "==================================================" 

else

# Checking if this is non-linear PB, non-linear PB with BCOPT=6 or linear PB
NPBOPT=$(zcat ${cmplx} | grep -c -i 'npbopt=1\|npbopt= 1\|npbopt=  1')
BCOPT=$(zcat ${cmplx}  | grep -c -i 'bcopt=6\|bcopt= 6\|bcopt=  6')

if [ "$NPBOPT" -gt "0" ] && [ "$BCOPT" -eq "0" ]
then
    echo "===================================================" 
    echo "Processing Non-linear PB calculations with BCOPT!=6" 
    echo "===================================================" 
fi
if [ "$NPBOPT" -gt "0" ] && [ "$BCOPT" -gt "0" ]
then
    echo "==================================================" 
    echo "Processing Non-linear PB calculations with BCOPT=6" 
    echo "==================================================" 
fi
if [ "$NPBOPT" -eq "0" ]
then 
    echo "=================================="
    echo "Processing linear PB calculations." 
    echo "=================================="
fi

fi

# Checking if this is QMMM 
NQMMM=$(zcat ${cmplx} | grep -c -i 'ifqnt=1\|ifqnt= 1\|ifqnt=  1')
if [ "$NQMMM" -gt "0" ]
then
    echo "===================================" 
    echo "Processing CMPLX QMMM calculations" 
    echo "===================================" 
# Checking if D3H4 correction is inlcuded
    ND3H4=$(zcat ${cmplx} |  grep -c 'D3H4')
    if [ "$ND3H4" -gt "0" ]
    then
      echo "=====================================================" 
      echo "ESCF in SCC-DFTB3 to be combined with D3H4 correction" 
      echo "=====================================================" 
      D3H4_INCL="YES"
    else
      D3H4_INCL="NO"
    fi
fi
declare -a NQMMM_FRAG=""
for ((ifrag=1;ifrag<=NFRAG;ifrag++))
do
  let " iout = $ifrag + 3 "
  frag=${prefix}.out${iout}.gz
  NQMMM_FRAG["$ifrag"]=$(zcat $frag | grep -c -i 'ifqnt=1\|ifqnt= 1\|ifqnt=  1')
  if [ ${NQMMM_FRAG["$ifrag"]} -gt "0" ]
  then 
    echo "===========================================" 
    echo "Processing FRAG $ifrag  QMMM calculations." 
    echo "===========================================" 
  fi
done 

# Output files 
rm -f ECOMP.dat ECOMP_DIFF.dat G_MMPB.dat G_MMPBSA.dat G_MMGBSA.dat NWAT.dat

if [ "$fragments" -eq "1" ]
then 
  EFRAGTXT=""
  GFRAGTXT=""
  NPFRAGTXT=""
  for ((ifrag=1;ifrag<=NFRAG;ifrag++))
  do
     if [ ${NQMMM_FRAG["$ifrag"]} -gt "0" ]
     then 
        EFRAGTXT="$EFRAGTXT EGAS(FRAG${ifrag}) ESCF(FRAG${ifrag}) PB(FRAG${ifrag}) GB(FRAG${ifrag}) NP(FRAG${ifrag})"
     else 
        EFRAGTXT="$EFRAGTXT EGAS(FRAG${ifrag}) PB(FRAG${ifrag}) GB(FRAG${ifrag}) NP(FRAG${ifrag})"
     fi
     NPFRAGTXT="$NPFRAGTXT AREA_MOLSURF(FRAG${ifrag}) GCAV_MOLSURF(FRAG${ifrag}) ECAV(FRAG${ifrag}) EDIS(FRAG${ifrag}) ASES_NUM(FRAG${ifrag}) VSES_NUM(FRAG${ifrag}) GCAV_CSPT(FRAG${ifrag})"
     GFRAGTXT="$GFRAGTXT G(FRAG${ifrag})" 
  done
  if [ "$NQMMM" -gt "0" ] 
  then 
    echo "# PREFIX  EGAS(CMPLX)  ESCF(CMPLX) PB(CMPLX) GB(CMPLX)  NP(CMPLX) $EFRAGTXT " > ECOMP.dat 
    echo '# Warning:  EGAS includes ESCF ' >> ECOMP.dat
    if [ $D3H4_INCL == "YES" ]; then echo '# Warning:  ESCF includes D3H4 ' >> ECOMP.dat ; fi
    echo '# PREFIX  EGAS(CMPLX) ESCF(CMPLX)  PB(CMPLX)  VDWint(CMPLX) GCAV_CSPT(CMPLX) GCAV_MOLSURF(CMPLX)  G(CMPLX)_SPT   G(CMPLX)_MOLSURF'  > G_MMPB.dat 
    echo '# Warning:  EGAS includes ESCF ' >> G_MMPB.dat
    echo '# PREFIX  diff_EGAS diff_ESCF  diff_PB diff_GB  diff_SURF  diff_G_PB  diff_G_GB'  > ECOMP_DIFF.dat 
    echo '# Warning:  diff_EGAS includes diff_ESCF ' >> ECOMP_DIFF.dat
  else
    echo "# PREFIX  EGAS(CMPLX)  PB(CMPLX) GB(CMPLX)  NP(CMPLX) $EFRAGTXT " > ECOMP.dat 
    echo '# PREFIX  EGAS(CMPLX)  PB(CMPLX)  VDWint(CMPLX) GCAV_CSPT(CMPLX) GCAV_MOLSURF(CMPLX)  G(CMPLX)_SPT  G(CMPLX)_MOLSURF '  > G_MMPB.dat 
    echo '# PREFIX  diff_EGAS  diff_PB diff_GB  diff_SURF  diff_G_PB  diff_G_GB'  > ECOMP_DIFF.dat 
  fi
  echo "# PREFIX VDWint(CMPLX)  AREA_MOLSURF(CMPLX) GCAV_MOLSURF(CMPLX)   ECAV(CMPLX) EDIS(CMPLX) ASES_NUM(CMPLX) VSES_NUM(CMPLX) GCAV_CSPT(CMPLX)  $NPFRAGTXT " >  NONPOLAR.dat 
  echo "# PREFIX  diff_AREA_MOLSURF diff_GCAV_MOLSURF  diff_ECAV  diff_EDIS diff_ASES_NUM diff_VSES_NUM diff_GCAV_CSPT " >  NONPOLAR_DIFF.dat 
  echo "# PREFIX  G(CMPLX) $GFRAGTXT diff_G "  > G_MMPBSA.dat 
  echo "# PREFIX  G(CMPLX) $GFRAGTXT diff_G "  > G_MMGBSA.dat 
  echo '# PREFIX  NWAT(CMPLX) '  > NWAT.dat 
else
  if [ "$NQMMM" -gt "0" ] 
  then 
    echo '# PREFIX  EGAS(CMPLX)  ESCF(CMPLX) PB(CMPLX) GB(CMPLX)  NP(CMPLX)   ' > ECOMP.dat 
    echo '# Warning:  EGAS includes ESCF ' >> ECOMP.dat
    if [ $D3H4_INCL == "YES" ]; then echo '# Warning:  ESCF includes D3H4 ' >> ECOMP.dat ; fi
    echo '# PREFIX  EGAS(CMPLX) ESCF(CMPLX)  PB(CMPLX)  VDWint(CMPLX) GCAV_CSPT(CMPLX) GCAV_MOLSURF(CMPLX)  G(CMPLX)_SPT   G(CMPLX)_MOLSURF '  > G_MMPB.dat 
    echo '# Warning:  EGAS includes ESCF ' >> G_MMPB.dat 
    if [ $D3H4_INCL == "YES" ]; then echo '# Warning:  ESCF includes D3H4 ' >> G_MMPB.dat ; fi
  else
    echo '# PREFIX  EGAS(CMPLX)  PB(CMPLX) GB(CMPLX)  NP(CMPLX)   ' > ECOMP.dat 
    echo '# PREFIX  EGAS(CMPLX)  PB(CMPLX)  VDWint(CMPLX) GCAV_CSPT(CMPLX) GCAV_MOLSURF(CMPLX)  G(CMPLX)_SPT   G(CMPLX)_MOLSURF '  > G_MMPB.dat 
  fi
  echo "# PREFIX VDWint(CMPLX) AREA_MOLSURF(CMPLX) GCAV_MOLSURF(CMPLX)  ECAV(CMPLX) EDIS(CMPLX) ASES_NUM(CMPLX) VSES_NUM(CMPLX) GCAV_CSPT(CMPLX) " >  NONPOLAR.dat 
  echo '# PREFIX  G_MMPBSA(CMPLX)   '  > G_MMPBSA.dat 
  echo '# PREFIX  G_MMGBSA(CMPLX)  '  > G_MMGBSA.dat 
  echo '# PREFIX  NWAT(CMPLX) '  > NWAT.dat 
fi

for file in $(cat $1)
do

# Filenames
  prefix=${file%%.*}
  cmplxwat=${prefix}.out1
  wat=${prefix}.out2
  cmplx=${prefix}.out3
 
  printf '\r%s:' "$prefix" 

# Get Egas, PB, GB and  Surf for CMPLX 
  GET_MMPBSA $cmplx
  EGAS_cmplx=$EGAS
  ESCF_cmplx=$ESCF
  GB_cmplx=$GB
  NP_cmplx=$NP
  ECAV_cmplx=$ECAV
  EDIS_cmplx=$EDIS
  PB_cmplx=$PB
  ASES_NUM_cmplx=$ASES_NUM
  VSES_NUM_cmplx=$VSES_NUM
  GCAV_CSPT_cmplx=$GCAV_CSPT
  AREA_MOLSURF_cmplx=$AREA_MOLSURF
  GCAV_MOLSURF_cmplx=$(echo $AREA_MOLSURF_cmplx | awk '{g=$1*0.069;printf(" %10.3f \n",g)}')
  if [ "$GCAV_MOLSURF_cmplx" == "" ] 
  then  
     GCAV_MOLSURF_cmplx="0.000"
  fi


if [ "$fragments" -eq "1" ]
then 
   declare -a EGAS_frag
   declare -a ESCF_frag
   declare -a GB_frag
   declare -a PB_frag
   declare -a NP_frag
   declare -a ECAV_frag
   declare -a EDIS_frag
   declare -a AREA_MOLSURF_frag
   declare -a GCAV_MOLSURF_frag
   declare -a ASES_NUM_frag
   declare -a VSES_NUM_frag
   declare -a GCAV_CSPT_frag

# Get Egas, GB and  Surf for HOST

  for ((ifrag=1;ifrag<=NFRAG;ifrag++))
  do
    let " iout = $ifrag + 3 "
    frag=${prefix}.out${iout} 
    GET_MMPBSA $frag
    EGAS_frag["$ifrag"]=$EGAS
    ESCF_frag["$ifrag"]=$ESCF 
    GB_frag["$ifrag"]=$GB
    NP_frag["$ifrag"]=$NP
    ECAV_frag["$ifrag"]=$ECAV
    EDIS_frag["$ifrag"]=$EDIS
    PB_frag["$ifrag"]=$PB
    AREA_MOLSURF_frag["$ifrag"]=$AREA_MOLSURF
    ASES_NUM_frag["$ifrag"]=$ASES_NUM
    VSES_NUM_frag["$ifrag"]=$VSES_NUM
    GCAV_CSPT_frag["$ifrag"]=$GCAV_CSPT
    GCAV_MOLSURF_frag["$ifrag"]=$(echo $AREA_MOLSURF | awk '{g=$1*0.069;printf(" %10.3f \n",g)}')
    if [ ${GCAV_MOLSURF_frag["$ifrag"]} == "" ] 
    then  
       GCAV_MOLSURF_frag["$ifrag"]="0.000"
    fi
  done 

fi

# Get solute-solvent VDW interaction energy
  NWAT=$(zcat $cmplxwat  | grep ' Number of triangulated 3-point waters found:'  | tail -1 | awk '{print $7}')
  EVDW_cmplxwat=$(zcat $cmplxwat  | grep 'VDWAALS ='  | tail -1 | awk '{print $3}')
  EVDW_wat=$(zcat $wat  | grep 'VDWAALS ='  | tail -1 | awk '{print $3}')
  EVDW_cmplx=$(zcat $cmplx  | grep 'VDWAALS ='  | tail -1 | awk '{print $3}')

  if [ "$EVDW_cmplxwat" == "" ] ||  [ "$EVDW_cmplx" == "" ]  ||  [ "$EVDW_wat" == "" ]
  then 
    EVDW_int="0.000"
  else 
#   EVDW_int="$EVDW_cmplxwat - $EVDW_wat - $EVDW_cmplx" 
    EVDW_int=$(echo "$EVDW_cmplxwat $EVDW_wat $EVDW_cmplx" | awk '{e= $1 - $2 - $3 ;printf("%12.4f",e)}')
  fi

# Computing G energies
 if [ "$EGAS_cmplx" != "0.000" ] &&  [ "$PB_cmplx" != "0.000" ] &&  [ "$NP_cmplx" != "0.000" ] 
 then 
# GCMPLX_PB= "$EGAS_cmplx + $PB_cmplx + $NP_cmplx"
  GCMPLX_PB=$(echo "$EGAS_cmplx $PB_cmplx $NP_cmplx " | awk '{g= $1 + $2 + $3 ;printf("%12.4f",g)}')
 else
    GCMPLX_PB="0.000"
 fi
 if [ "$EGAS_cmplx" != "0.000" ] &&  [ "$GB_cmplx" != "0.000" ] &&  [ "$NP_cmplx" != "0.000" ] 
 then 
#   GCMPLX_GB="$EGAS_cmplx + $GB_cmplx + $NP_cmplx"
    GCMPLX_GB=$(echo "$EGAS_cmplx $GB_cmplx $NP_cmplx " | awk '{g= $1 + $2 + $3 ;printf("%12.4f",g)}')
 else
    GCMPLX_GB="0.000"
 fi
# GCMPLX_MMPB_SPT="$EGAS_cmplx + $PB_cmplx+ $EVDW_int + $GCAV_CSPT"
  GCMPLX_MMPB_SPT=$(echo "$EGAS_cmplx $PB_cmplx $EVDW_int $GCAV_CSPT_cmplx " | awk '{g= $1 + $2 + $3 +$4 ;printf("%12.4f",g)}')
# GCMPLX_MMPB_MOLSRUF="$EGAS_cmplx + $PB_cmplx+ $EVDW_int + $GCAV_MOLSURF"
  GCMPLX_MMPB_MOLSURF=$(echo "$EGAS_cmplx $PB_cmplx $EVDW_int $GCAV_MOLSURF_cmplx " | awk '{g= $1 + $2 + $3 +$4 ;printf("%12.4f",g)}')

if [ "$fragments" -eq "1" ]
then 

 declare -a GFRAG_PB=""
 declare -a GFRAG_GB=""

 SUM_G_PB="0.0"
 SUM_G_GB="0.0"
 SUM_PB="0.0"
 SUM_GB="0.0"
 SUM_EGAS="0.0"
 SUM_ESCF="0.0"
 SUM_NP="0.0"
 SUM_ECAV="0.0"
 SUM_EDIS="0.0"
 SUM_AREA_MOLSURF="0.0"
 SUM_GCAV_MOLSURF="0.0"
 SUM_ASES_NUM="0.0"
 SUM_VSES_NUM="0.0"
 SUM_GCAV_CSPT="0.0"

 for ((ifrag=1;ifrag<=NFRAG;ifrag++))
 do 

 TER_EGAS=${EGAS_frag["$ifrag"]}
 TER_ESCF=${ESCF_frag["$ifrag"]}
 TER_PB=${PB_frag["$ifrag"]}
 TER_GB=${GB_frag["$ifrag"]}
 TER_NP=${NP_frag["$ifrag"]}
 TER_ECAV=${ECAV_frag["$ifrag"]}
 TER_EDIS=${EDIS_frag["$ifrag"]}
 TER_AREA_MOLSURF=${AREA_MOLSURF_frag["$ifrag"]}
 TER_GCAV_MOLSURF=${GCAV_MOLSURF_frag["$ifrag"]}
 TER_ASES_NUM=${ASES_NUM_frag["$ifrag"]}
 TER_VSES_NUM=${VSES_NUM_frag["$ifrag"]}
 TER_GCAV_CSPT=${GCAV_CSPT_frag["$ifrag"]}

 if [ ${EGAS_frag["$ifrag"]} != "0.000" ] &&  [ ${PB_frag["$ifrag"]} != "0.000" ] &&  [ ${NP_frag["$ifrag"]} != "0.000" ] 
 then 
    GFRAG_PB["$ifrag"]=$(echo "$TER_EGAS  $TER_PB  $TER_NP " | awk '{g= $1 + $2 + $3 ;printf("%12.4f",g)}')
 else
    GFRAG_PB["$ifrag"]="0.000"
 fi
 if [ ${EGAS_frag["$ifrag"]} != "0.000" ] &&  [ ${GB_frag["$ifrag"]} != "0.000" ] &&  [ ${NP_frag["$ifrag"]} != "0.000" ] 
 then 
    GFRAG_GB["$ifrag"]=$(echo "$TER_EGAS  $TER_GB $TER_NP " | awk '{g= $1 + $2 + $3 ;printf("%12.4f",g)}')
 else
    GFRAG_GB["$ifrag"]="0.000"
 fi

 TER_G_PB=${GFRAG_PB["$ifrag"]}
 TER_G_GB=${GFRAG_GB["$ifrag"]}
 SUM_G_PB=$(echo "$SUM_G_PB $TER_G_PB " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_G_GB=$(echo "$SUM_G_GB $TER_G_GB " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
  

 SUM_PB=$(echo "$SUM_PB $TER_PB " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_GB=$(echo "$SUM_GB $TER_GB " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_EGAS=$(echo "$SUM_EGAS $TER_EGAS " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_ESCF=$(echo "$SUM_ESCF $TER_ESCF " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_NP=$(echo "$SUM_NP $TER_NP " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_ECAV=$(echo "$SUM_ECAV $TER_ECAV " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_EDIS=$(echo "$SUM_EDIS $TER_EDIS " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_AREA_MOLSURF=$(echo "$SUM_AREA_MOLSURF $TER_AREA_MOLSURF " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_GCAV_MOLSURF=$(echo "$SUM_GCAV_MOLSURF $TER_GCAV_MOLSURF " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_ASES_NUM=$(echo "$SUM_ASES_NUM $TER_ASES_NUM " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_VSES_NUM=$(echo "$SUM_VSES_NUM $TER_VSES_NUM " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')
 SUM_GCAV_CSPT=$(echo "$SUM_GCAV_CSPT $TER_GCAV_CSPT " | awk '{g= $1 + $2  ;printf("%12.4f",g)}')

done 

# Computing Differences

 G_PB_diff=$(echo "$GCMPLX_PB $SUM_G_PB" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 G_GB_diff=$(echo "$GCMPLX_GB $SUM_G_GB" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 PB_diff=$(echo "$PB_cmplx  $SUM_PB" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 GB_diff=$(echo "$GB_cmplx  $SUM_GB" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 EGAS_diff=$(echo "$EGAS_cmplx  $SUM_EGAS" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 ESCF_diff=$(echo "$ESCF_cmplx  $SUM_ESCF" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 NP_diff=$(echo "$NP_cmplx  $SUM_NP" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 ECAV_diff=$(echo "$ECAV_cmplx  $SUM_ECAV" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 EDIS_diff=$(echo "$EDIS_cmplx  $SUM_EDIS" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 AREA_MOLSUR_diff=$(echo "$AREA_MOLSUR_cmplx  $SUM_AREA_MOLSUR" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 GCAV_MOLSUR_diff=$(echo "$GCAV_MOLSUR_cmplx  $SUM_GCAV_MOLSUR" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 ASES_NUM_diff=$(echo "$ASES_NUM_cmplx  $SUM_ASES_NUM" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 VSES_NUM_diff=$(echo "$VSES_NUM_cmplx  $SUM_VSES_NUM" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 
 GCAV_CSPT_diff=$(echo "$GCAV_CSPT_cmplx  $SUM_GCAV_CSPT" | awk '{g= $1 - $2 ;printf("%12.4f",g)}') 


# Printing out data
  if [ "$NQMMM" -gt  "0" ] 
  then 
     EFRAGTXT=""
     PBFRAGTXT=""
     GBFRAGTXT=""
     NPFRAGTXT=""
     for ((ifrag=1;ifrag<=NFRAG;ifrag++))
     do
       TER_EGAS=${EGAS_frag["$ifrag"]}
       TER_ESCF=${ESCF_frag["$ifrag"]}
       TER_PB=${PB_frag["$ifrag"]}
       TER_GB=${GB_frag["$ifrag"]}
       TER_NP=${NP_frag["$ifrag"]}
       TER_ECAV=${ECAV_frag["$ifrag"]}
       TER_EDIS=${EDIS_frag["$ifrag"]}
       EFRAGTXT="$EFRAGTXT $TER_EGAS $TER_ESCF $TER_PB $TER_GB $TER_NP "
       NPFRAGTXT="$NPFRAGTXT $TER_ECAV $TER_EDIS $TER_ASES_NUM $TER_VSES_NUM $TER_GCAV_CSPT "
       TER_G_PB=${GFRAG_PB["$ifrag"]}
       TER_G_GB=${GFRAG_GB["$ifrag"]}
       PBFRAGTXT="$PBFRAGTXT $TER_G_PB "
       GBFRAGTXT="$GBFRAGTXT $TER_G_GB "
     done
     echo "$prefix $EGAS_cmplx $ESCF_cmplx $PB_cmplx $GB_cmplx $NP_cmplx $EFRAGTXT" >> ECOMP.dat
     echo "$prefix $EVDW_int  $AREA_MOLSURF_cmplx $GCAV_MOLSURF_cmplx  $ECAV_cmplx $EDIS_cmplx $ASES_NUM_cmplx $VSES_NUM_cmplx $GCAV_CSPT_cmplx  $NPFRAGTXT" >>  NONPOLAR.dat
     echo "$prefix $EGAS_diff  $ESCF_diff $PB_diff $GB_diff $NP_diff $G_PB_diff $G_GB_diff" >> ECOMP_DIFF.dat
     echo "$prefix $ECAV_diff  $EDIS_diff  $ASES_NUM_diff $VSES_NUM_diff $GCAV_CSPT_diff " >>  NONPOLAR_DIFF.dat
     echo "$prefix $EGAS_cmplx $ESCF_cmplx $PB_cmplx $EVDW_int $GCAV_CSPT_cmplx $GCAV_MOLSURF_cmplx $GCMPLX_MMPB_SPT  $GCMPLX_MMPB_MOLSURF" >> G_MMPB.dat
  else
     EFRAGTXT=""
     PBFRAGTXT=""
     NPFRAGTXT=""
     GBFRAGTXT=""
     for ((ifrag=1;ifrag<=NFRAG;ifrag++))
     do
       TER_EGAS=${EGAS_frag["$ifrag"]}
       TER_PB=${PB_frag["$ifrag"]}
       TER_GB=${GB_frag["$ifrag"]}
       TER_NP=${NP_frag["$ifrag"]}
       TER_ECAV=${ECAV_frag["$ifrag"]}
       TER_EDIS=${EDIS_frag["$ifrag"]}
       NPFRAGTXT="$NPFRAGTXT $TER_ECAV $TER_EDIS $TER_ASES_NUM $TER_VSES_NUM $TER_GCAV_CSPT "
       EFRAGTXT="$EFRAGTXT $TER_EGAS $TER_PB $TER_GB $TER_NP "
       TER_G_PB=${GFRAG_PB["$ifrag"]}
       TER_G_GB=${GFRAG_GB["$ifrag"]}
       PBFRAGTXT="$PBFRAGTXT $TER_G_PB "
       GBFRAGTXT="$GBFRAGTXT $TER_G_GB "
     done
     echo "$prefix $EGAS_cmplx $PB_cmplx $GB_cmplx $NP_cmplx $EFRAGTXT" >> ECOMP.dat
     echo "$prefix $EVDW_int  $AREA_MOLSURF_cmplx $GCAV_MOLSURF_cmplx $ECAV_cmplx $EDIS_cmplx $ASES_NUM_cmplx $VSES_NUM_cmplx $GCAV_CSPT_cmplx  $NPFRAGTXT" >>  NONPOLAR.dat
     echo "$prefix $EGAS_diff $PB_diff $GB_diff $NP_diff $G_PB_diff $G_GB_diff" >> ECOMP_DIFF.dat
     echo "$prefix $AREA_MOLSURF_diff $GCAV_MOLSURF_diff $ECAV_diff  $EDIS_diff  $ASES_NUM_diff $VSES_NUM_diff $GCAV_CSPT_diff " >>  NONPOLAR_DIFF.dat
     echo "$prefix $EGAS_cmplx $PB_cmplx $EVDW_int $GCAV_CSPT_cmplx $GCAV_MOLSURF_cmplx $GCMPLX_MMPB_SPT  $GCMPLX_MMPB_MOLSURF" >> G_MMPB.dat
  fi 
  echo "$prefix $GCMPLX_PB $PBFRAGTXT $G_PB_diff" >> G_MMPBSA.dat
  echo "$prefix $GCMPLX_GB $GBFRAGTXT $G_GB_diff" >> G_MMGBSA.dat
  echo "$prefix $NWAT" >> NWAT.dat

else 

  if [ "$NQMMM" -gt  "0" ]  
  then 
    echo "$prefix $EGAS_cmplx $ESCF_cmplx $PB_cmplx $GB_cmplx $NP_cmplx" >> ECOMP.dat
    echo "$prefix $EGAS_cmplx $ESCF_cmplx $PB_cmplx $EVDW_int $GCAV_CSPT_cmplx $GCAV_MOLSURF_cmplx $GCMPLX_MMPB_SPT  $GCMPLX_MMPB_MOLSURF" >> G_MMPB.dat
  else
    echo "$prefix $EGAS_cmplx $PB_cmplx $GB_cmplx $NP_cmplx " >> ECOMP.dat
    echo "$prefix $EGAS_cmplx $PB_cmplx $EVDW_int $GCAV_CSPT_cmplx $GCAV_MOLSURF_cmplx $GCMPLX_MMPB_SPT  $GCMPLX_MMPB_MOLSURF" >> G_MMPB.dat
  fi
  echo "$prefix $EVDW_int $AREA_MOLSURF_cmplx $GCAV_MOLSURF_cmplx $ECAV_cmplx $EDIS_cmplx $ASES_NUM_cmplx $VSES_NUM_cmplx $GCAV_CSPT_cmplx  " >>  NONPOLAR.dat
  echo "$prefix $GCMPLX_PB " >> G_MMPBSA.dat
  echo "$prefix $GCMPLX_GB " >> G_MMGBSA.dat
  echo "$prefix $NWAT" >> NWAT.dat

fi

done

fi

printf '\n'

if [ ! "$DO_ENERGY_ANALYSES_TAR" ]
then
  echo " Shall we tar all output files in OUTPUT.tar?  (Y/N) ?"
  read YES
  
  if [ "$YES" == "y" ] || [ "$YES" == "Y" ]
  then 
    DO_ENERGY_ANALYSES_TAR="1"
  else
    DO_ENERGY_ANALYSES_TAR="0"
  fi

fi

if [ "$DO_ENERGY_ANALYSES_TAR" -eq "1" ]
then
   rm -f OUTPUT.tar
   tar cvf OUTPUT.tar *.out1.gz
   tar rvf OUTPUT.tar *.out2.gz
   tar rvf OUTPUT.tar *.out3.gz
   if [ $NFRAG -gt 1 ] 
   then 
       for ((i=1;i<=NFRAG;i++)); do let "j=$i+3"; tar rvf OUTPUT.tar *.out${j}.gz; done
   fi
   if [ "$NQMMM" -gt "0" ]
   then 
      tar rvf OUTPUT.tar *.out_relax_qmmm.gz
   fi
   rm -f *.out?.gz *.out??.gz
   if [ "$NQMMM" -gt "0" ]
   then 
      rm -f *.out_relax_qmmm.gz
   fi
   tar rvf OUTPUT.tar *.pdb3.gz
   rm -f *.pdb3.gz
   echo "All output files are archived into OUTPUT.tar" 
fi 


  
