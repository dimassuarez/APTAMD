#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
source $APTAMD/ENV/aptamd_env.sh

declare -a IARG=""

IARG=($*)
NARG=${#IARG[@]} 

let "NFRAG=$NARG-1"

# Last argument is PDB file to be processed
# The other arguments specify the number of ions to
# be captured for each fragment 

NUM_INA=0
declare -a INA_FRAG=""

# Counting how many INA ions are to be retained
rm -f tmp_na_list
for ((IFRAG=0;IFRAG<=NFRAG-1;IFRAG++))
do
   INA_FRAG["$IFRAG"]=${IARG["$IFRAG"]}
   let "NUM_INA=${NUM_INA}+${IARG["$IFRAG"]}"
   echo "$IFRAG ${IARG["$IFRAG"]}" >> tmp_na_list
done
# PDB file to be processed
PDB=${IARG["$NFRAG"]} 
echo "PDB=$PDB"

# INA_CHECK annotates which INA ions are selected
# For simplicity, the array extends to all solute + ion atoms
NAT=$(grep -v WAT $PDB | grep -c 'ATOM ' )
declare -a INA_CHECK=""
for ((IAT=0;IAT<=NAT;IAT++))
do
   INA_CHECK["$IAT"]=0
done

# Put aside solut atoms
grep -v 'Na+\|Cl-\|WAT\|HOH' $PDB > temp_prep.pdb

# For each fragment, we pick up its INA ions as determined
# by the cpptraj analysis summarized in the NA_SORTED_FRAG_X.INFO files
for IFRAG in $(sort -n -k 2 tmp_na_list | awk '{print $1}')
do

let "JFRAG=$IFRAG+1"
rm -f tmp_frag_ina_${JFRAG}
touch tmp_frag_ina_${JFRAG}

if [ ${INA_FRAG["$IFRAG"]} -gt 0 ]
then

declare -a INA=""
let "JCOL=$NUM_INA+1"
INA=($(grep "$PDB" $WORKDIR/NA_SORTED_${JFRAG}.INFO | cut -d, -f2-${JRES} | sed 's/,/  /g'))
#echo "${INA[*]}"

k=0
for ((I=0;I<=NUM_INA-1;I++))
do
   na=${INA["$I"]} 
   if [ ${INA_CHECK["$na"]} -eq  0 ]
   then
      let "k=$k+1"
      iline=$(grep -n "${na} Na+" $PDB | sed 's/:/ /' | awk '{print $1}')
      sed -i "${iline},${iline}s/Na+/INA/g" $PDB 
      sed -n "${iline},${iline}p" $PDB  >> tmp_frag_ina_${JFRAG}
      INA_CHECK["$na"]=1
      if [ ${k} -eq ${INA_FRAG["$IFRAG"]} ]; then break; fi
   fi
done

unset INA

fi

done

for ((JFRAG=1;JFRAG<=NFRAG;JFRAG++))
do
    cat tmp_frag_ina_${JFRAG} >> temp_prep.pdb 
done

# All solvent molecules are put back into the final PDB
grep 'WAT\|HOH' $PDB >> temp_prep.pdb
$TOOLS/reorder < temp_prep.pdb >  $PDB

rm -f temp_prep.pdb tmp_na_list tmp_frag_ina_*
