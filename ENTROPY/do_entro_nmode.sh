#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
if [ "$#" -eq 0 ]; then more $APTAMD/DOC/do_entro_nmode.txt ; exit; fi

source $APTAMD/ENV/aptamd_env.sh

# Checking options
extension="${1##*.}"
if [ $extension == "src" ]
then
    echo "Sourcing $1"
    source $1
else
#  Put into double quotes the list of trajectories to be processed
   MD_TRAJ=$1
#  GAMD or MD trajectory type should be specified 
   MD_TYPE=$2
#  Number of NA ions to be included
   NSODIUM_LIMIT=$3
fi

if [ ! -z "$MOL" ] && [  -z "$MD_TRAJ" ]; then MD_TRAJ="${MOL}"; echo $MD_TRAJ; unset MOL; fi

if [ -z "$MD_TRAJ" ]; then more $APTAMD/DOC/do_entro_nmode.txt ; exit ; fi
if [ -z "$MD_TYPE" ]; then more $APTAMD/DOC/do_entro_nmode.txt ; exit ; fi

# For compatibility with former datasets
if [ ! -z "$NSODIUM_LIMIT" ] &&  [ -z "$CNTION_LIST" ]  &&  [ -z "$NCNTION_LIMIT" ]
then
        CNTION_LIST="Na+"
        ix=0
        NCNTION_LIMIT=""
        for x in $(echo $NSODIUM_LIMIT)
        do
                let "ix=$ix+1"
                if [ $ix -eq 1 ]
                then
                NCNTION_LIMIT="${x}"
                else
                NCNTION_LIMIT="${NCNTION_LIMIT}:${x}"
                fi
        done
        if [ ! -z "$SODIUM_FRAG" ]  &&  [ -z "$CNTION_FRAG" ]; then CNTION_FRAG="$SODIUM_FRAG"; fi
        echo "Interpreting NSODIUM_LIMIT=$NSODIUM_LIMIT as "
        echo "             CNTION_LIST=Na+"
        echo "             NCNTION_LIMIT=${NCNTION_LIMIT}"
        echo "Avoid declaring NSODIUM_LIMIT"
fi


#  List of counterions
if [ -z "$CNTION_LIST" ]
then
        declare -a CNTION=("Na+")
else
        declare -a CNTION=($CNTION_LIST)
fi
NCNTION=${#CNTION[@]}
if [ -z "$NCNTION_LIMIT" ]
then
  declare -a NCNTION_LIMIT=(0)
else
  declare -a NCNTION_LIMIT=(${NCNTION_LIMIT})
fi
if [ ${NCNTION} -eq ${#NCNTION_LIMIT[@]} ]
then
     for ((I=0;I<=NCNTION-1;I++))
     do
        echo "${NCNTION_LIMIT["$I"]} ${CNTION["$I"]} ions to be selected "
     done
else
      echo "CNTION_LIST=$CNTION_LIST  and NCNTION_LIMIT=${NCNTION_LIMIT[*]} =not compatible!"
      exit
fi

# CNTION ions can be separated fragments or just aggregated to cmplx/receptor/ligand fragments
if [ -z "$CNTION_FRAG" ]; then CNTION_FRAG="NO"; fi
if [ "$CNTION_FRAG" == "NO" ] && [ $NCNTION_LIMIT  -gt  0 ]
then
  echo "counterions and hydrating waters are assigned to CMPLX/FRAGs"
elif [ "$CNTION_FRAG" == "YES" ] && [ ${NCNTION_LIMIT[0]}  -gt  0 ]
then
  echo "counterions and hydrating waters are treated as additional fragments."
fi

if [ -z "$NPROCS" ]
then
   NPROCS=$(cat /proc/cpuinfo | grep -c processor)
   echo "Using $NPROCS available processors"
else
   echo "Using $NPROCS processors as predefined "
fi

if [ -z "$NUMPDB_INCR" ]
then
   INCRLIST="NO"
else
   INCRLIST="YES"
   echo "Increment in the List of PDB files"
   echo "Processing more PDB files NUMPDB_INCR=$NUMPDB_INCR "
fi
if [ -z "$SIEVE" ]
then
   SIEVE="10"
   echo "Processing 1 out of $SIEVE snapshots"
else
   echo "Processing 1 out of $SIEVE snapshots as predefined"
fi    
if [ -z "$FINERSIEVE" ]
then
   FINERSIEVE="NO"
fi    
if [ "$FINERSIEVE" == "YES" ]
then
   echo "Former OUTPUT files are kept (useful for finer sieve parameter executions)"
fi    
if [ "$FINERSIEVE" == "YES" ] && [ $INCRLIST == "YES" ]
then
   echo "FINERSIEVE and INCRLIST options cannot be both activated"
   exit
fi    
if [ -z "$TEMPERATURE" ]
then
   TEMPERATURE="300.0"
fi    
echo "TEMPERATURE=${TEMPERATURE}"

# Although we read PDB files...the origin of those PDB files may be different
if [ -z "$SUFFIX_MDCRD" ]; then echo 'Considering SUFFIX_MDCRD=_solutewat.mdcrd';  SUFFIX_MDCRD="_solutewat.mdcrd"; else echo "Using SUFFIX_MDCRD=$SUFFIX_MDCRD"; fi
if [ -z "$PREFIX_MDCRD" ]; then echo 'Considering PREFIX_MDCRD=md_';  PREFIX_MDCRD="md_"; else echo "Using PREFIX_MDCRD=$PREFIX_MDCRD"; fi

if [ -z $SNAPSHOTS_DIR ]
then
 SNAPSHOTS_DIR="SNAPSHOTS"
 echo "Using SNAPSHOTS_DIR=${SNAPSHOTS_DIR}"
else
 echo "Using SNAPSHOTS_DIR=${SNAPSHOTS_DIR} as predefined"
fi
if [ -z "$ENTRO_DIR" ]; then ENTRO_DIR="ENTRO_NMODE"; fi

# Solvent mask including counterions
if [ -z $SLVNTMASK ]
then
   SLVNTMASK=":WAT,Na+,Cl-,MG"
fi
GREPSLVNTMASK=$(echo $SLVNTMASK | sed 's/://' | sed 's/,/\\|/g')

# Counting  ions
declare -a NCNTION_TOT=""
declare -a NCNTION_FRAG=""
for ((I=0;I<=NCNTION-1;I++))
do
        XTEMP=${NCNTION_LIMIT["$I"]}
        J=0
        K=0
        for X in ${XTEMP//:/ }
        do
                let "J=$J+$X"
                let "K=$K+1"
        done
        NCNTION_TOT["$I"]=$J
        NCNTION_FRAG["$I"]=$K
        echo "Considering a total of $J ions of type ${CNTION["$I"]} distributed on $K fragments"
done

WORKDIR_TRJ=$PWD

for MOL in $MD_TRAJ 
do
  if [ -e ${MOL}_${MD_TYPE} ]
  then
    cd ${MOL}_${MD_TYPE} 
  else
    echo "${MOL}_${MD_TYPE} does not exist. Exiting!"
    exit
  fi

  if [ ${SUFFIX_MDCRD} == ".mdcrd" ]
  then
     TOPOLOGY=${MOL}.top
  elif [ ${SUFFIX_MDCRD} == "_solute.mdcrd" ]
  then
     TOPOLOGY=${MOL}_solute.top
  elif [ ${SUFFIX_MDCRD} == "_solutewat.mdcrd" ]
  then
     TOPOLOGY=${MOL}_solutewat.top
  else
     echo "Unexpected SUFFIX_MDCRD=${SUFFIX_MDCRD}"
     if [ -z $TOPOLOGY ]
     then
        TOPOLOGY=${MOL}.top
        echo "Trying TOPOLOGY=$TOPOLOGY, but maybe it does not work!"
     fi
  fi
  echo "Selecting TOPOLOGY=$TOPOLOGY for SUFFIX_MDCRD=${SUFFIX_MDCRD}"

  if [ ! -e $TOPOLOGY ]
  then 
      echo "$TOPOLOGY not found in ${MOL}_MD. Exiting"
      exit
  fi
  SOLUTE_CRD=${MOL}_solute.crd 
  SOLUTE_TOP=${MOL}_solute.top
  SOLUTE_PDB=${MOL}_solute.pdb
  if [ ! -e 1.EDITION/$SOLUTE_CRD ] || [ ! -e 1.EDITION/$SOLUTE_TOP ]
  then 
     echo "$SOLUTE_CRD or $SOLUTE_TOP not found in ${MOL}_MD/1.EDITION. Exiting"
     exit
  else
     cd 1.EDITION
     CMAP=$(grep -i -c cmap $SOLUTE_TOP) 
     if [ $CMAP -gt 0 ]
     then 
         echo "WARNING:"
         echo "$SOLUTE_TOP contains CMAP dihedral parameters."
         echo "These parameters cannot be processed by the NMODE program"
         echo "However, since the internal conformation of the solute molecule"
         echo "is preserved by the fixed solvent layer, the NMODE calculations"
         echo "are performed anyway."
     fi
     $AMBERHOME/bin/cpptraj $SOLUTE_TOP <<EOF > mlog 
trajin $SOLUTE_CRD
trajout $SOLUTE_PDB  pdb pdbatom
go
EOF


  NRES=$(grep 'ATOM  ' $SOLUTE_PDB   | tail -1 | awk '{print $5}')
  NFRAG=$(grep -c 'TER' $SOLUTE_PDB)
  declare -a IRES=""
  declare -a JRES=""
  JRES=($(grep TER $SOLUTE_PDB | awk '{print $4}'))
  IRES[0]="1"
  for ((i=1;i<=NFRAG-1;i++))
  do  
      let "j=$i-1"
      let "resid=${JRES["$j"]}+1"
      IRES["$i"]=$resid
  done
  echo "Detected $NFRAG fragments"
  for ((i=0; i<=NFRAG-1;i++))
  do
         echo  "Fragment :${IRES["$i"]}-${JRES["$i"]}"
  done
  cd ../

  fi

  #  Checking CNTION frag info
  for ((I=0;I<=NCNTION-1;I++))
  do
    if [ $NFRAG -ne ${NCNTION_FRAG["$I"]}  ]
    then
      echo "Detected $NFRAG fragments in $MOL, but NCNTION_LIMIT of ${CNTION["$I"]} ions is specified for $NCNTION_FRAG["$I"] fragments"
      if [ $NCNTION_FRAG["$I"] -eq 1 ]
      then
           echo "Assuming that NCNTION=${NCNTION_LIMIT["$I"]} applies to FRAG 1"
      else
           exit
      fi
    fi
  done

  NRES_NOION=${NRES}
  NFRAG_NOION=${NFRAG}
  NCNTION_ACCUM=0
  rm -f tmp_CNTION.info
  declare -a FRAG_NWAT=""
  declare -a FRAG_WATMASK=""
  for ((i=0;i<=20;i++))
  do
          FRAG_NWAT["$i"]=0
          FRAG_WATMASK["$i"]=""
  done


  for ((I=0;I<=NCNTION-1;I++))
  do
  LINE_CNTION_INFO=""
  if [ ${NCNTION_TOT["$I"]} -gt 0  ]
  then
       ion=${CNTION["$I"]}

	if [ ${ion} == 'Na+' ]
	then
          ionup="INA"
	elif [ ${ion} == 'K+' ]
	then
          ionup="IK"
	elif [ ${ion} == 'MG' ]
	then
          ion="MG "
          ionup="IXG"
	fi

       LINE_CNTION_INFO=" ${CNTION["$I"]} ${ionup}  "

       echo "Fixing $TOPOLOGY original label = $ion new label = $ionup "
#
#  We process the TOPOLOGY file to modify the atom/resname for Na+, MG, ...
#  ---> ION ( ION atoms will be kept in ENTRO analyses)
#
       echo  "NCNTION=${NCNTION_TOT["$I"]} ${CNTION["$I"]} ions  to be considered as solute atoms"
       if [ ! -e ${MOL}_solutewat_ION.top ]
       then
          cp $TOPOLOGY ${MOL}_solutewat_ION.top
          TOPOLOGY="${MOL}_solutewat_ION.top"
       fi

       NTOP_CNTION=$(sed  '/FLAG CHARGE/,$d' $TOPOLOGY | tr " \t" "\n" | grep -c "${CNTION["$I"]}" )

       if [ ${NTOP_CNTION} -lt ${NCNTION_TOT["$I"]} ] 
       then
          echo "Not enough ${CNTION["$I"]} in $TOPOLOGY  NCNTION=${NCNTION_TOT["$I"]}  NTOP_CNTION=${NTOP_CNTION}"
          exit
       fi

       sed -n '1,/FLAG CHARGE/p' $TOPOLOGY > temp_A.top
       sed -n '/FLAG CHARGE/,/FLAG AMBER_ATOM_TYPE/p' $TOPOLOGY  | sed '1,1d' > temp_B.top
       sed  '1,/FLAG AMBER_ATOM_TYPE/d' $TOPOLOGY > temp_C.top

       for file in temp_A.top temp_B.top 
       do
            declare -a iline=""
            iline=($(grep -n "${CNTION["$I"]}" $file | sed 's/:/  /' | awk '{print $1}'))
            nlines=${#iline[@]} 
            NUM_CNT=0
            for ((i=1;i<=nlines;i++))
            do
               let "j=$i-1"
               jline=${iline["$j"]}
               NCNTION_LINE=$(sed -n "${jline},${jline}p"  $file | tr " \t" "\n" | grep -c "${CNTION["$I"]}") 
               if [ ${NUM_CNT} -lt ${NCNTION_TOT["$I"]} ]
               then
                   let "NUM_CNT_NEED=${NCNTION_TOT["$I"]}-${NUM_CNT}"
                   if [ ${NCNTION_LINE} -le ${NUM_CNT_NEED} ]
                   then 
                       sed -i "${jline}s/${ion}/${ionup}/g" $file
                       let "NUM_CNT=${NUM_CNT} + ${NCNTION_LINE}"
                   else
                       for ((k=1;k<=NUM_CNT_NEED;k++)); do sed -i "${jline}s/${ion}/${ionup}/" $file; done
                       let "NUM_CNT=${NUM_CNT} + ${NUM_CNT_NEED}"
                   fi
               else
                   break
               fi
             done
       done
       cat temp_A.top temp_B.top temp_C.top > $TOPOLOGY
       rm -f temp_A.top temp_B.top temp_C.top

#      Fixing some variables according to the number of counterions
       ENTRO_DIR=${ENTRO_DIR}_${ionup/I/}_${NCNTION_LIMIT["$I"]}

       declare -a FRAG_NWAT=""
       declare -a FRAG_WATMASK=""

       if [ ${CNTION_FRAG} == "YES" ]
       then

          for ((i=1;i<=${NCNTION_TOT["$I"]};i++))
          do
             let " j=${NRES} + $i "
             let " k=${NFRAG} -1 + $i "
             IRES["$k"]=$j
             JRES["$k"]=$j
             FRAG_NWAT["$k"]=6
             FRAG_WATMASK["$k"]=":${j}"
          done
          let "NFRAG=${NFRAG}+${NCNTION_TOT["$I"]}"
          let "NRES=${NRES}+${NCNTION_TOT["$I"]}"
	  LINE_CNTION_INFO="${LINE_CNTION_INFO}  ${NCNTION_TOT["$I"]} "
     
       else

        for ((IFRAG=0;IFRAG<=NFRAG-1;IFRAG++))
        do

         let "JFRAG=$IFRAG+1"

	 if [ ${NCNTION_FRAG["$I"]} -eq 1 ] && [ $IFRAG -gt 1 ]
         then 
            NCNTION_FRAG_SPLIT=0
	 else 
            NCNTION_FRAG_SPLIT=$(echo ${NCNTION_LIMIT["$I"]} | sed 's/:/\n/g' | head -${JFRAG} | tail -1)
	 fi

         if [ ${NCNTION_FRAG_SPLIT} -gt 0 ] 
         then 
            let " j=${NRES} + 1 "
            if [ ${I} -eq 0 ]  
	    then 
		 watmask=" :${j} "  
		 txt=" 6 "   # first cntion type
            else
                 txt="${txt} ; 6 "
                 watmask="${watmask} ; :${j} "
	    fi
            for ((i=2;i<=NCNTION_FRAG_SPLIT;i++))
            do
               let " j=${NRES} + $i + ${NCNTION_ACCUM}"
               txt="${txt} ; 6 "
               watmask="${watmask} ; :${j} "
            done
            FRAG_NWAT["$JFRAG"]=${txt}
            FRAG_WATMASK["$JFRAG"]=${watmask}
            let "NCNTION_ACCUM=${NCNTION_ACCUM}+${NCNTION_FRAG_SPLIT}"
	    LINE_CNTION_INFO="${LINE_CNTION_INFO}  ${NCNTION_FRAG_SPLIT} "
         else
            FRAG_NWAT["$JFRAG"]=""
            FRAG_WATMASK["$JFRAG"]=""
	    LINE_CNTION_INFO="${LINE_CNTION_INFO}  0 "
         fi

        done

        let "NRES=${NRES}+${NCNTION_TOT["$I"]}"  

       fi

       echo "${LINE_CNTION_INFO}" >> tmp_CNTION.info
       
  fi

  done  # End of loop over NCNTION
       
       
  if [ ! -e 6.ANALYSIS ]
  then 
     echo "6.ANALYSIS directory not found in ${MOL}_MD."
     exit             
  fi
  cd 6.ANALYSIS
  if [ ! -e $SNAPSHOTS_DIR ]
  then
      echo "$SNAPSHOTS_DIR directory not found in ${MOL}_MD/6.ANALYSIS"
      if [ $NSODIUM -eq 0 ]
      then
          SNAPSHOTS_DIR="SNAPSHOTS_NA"
          if [ ! -e $SNAPSHOTS_DIR ]
          then
              echo "$SNAPSHOTS_DIR directory not found either in ${MOL}_MD/6.ANALYSIS"
              exit
          else
             echo "Using $SNAPSHOTS_DIR directory "
          fi
      else
          exit
      fi
  else 
     cd $SNAPSHOTS_DIR
     ls ${MOL}*.pdb.gz | sed 's/.gz//' > LISTA
     cd ../
  fi

  if [ -e  ${ENTRO_DIR}  ]
  then 
      echo " ${ENTRO_DIR} found in ${MOL}_MD/6.ANALYSIS"
  else
      mkdir ${ENTRO_DIR}
  fi 

  cd ${ENTRO_DIR}
  WORKDIR=$PWD 
  cp ../${SNAPSHOTS_DIR}/LISTA .

  NPDB=$(cat LISTA | wc -l)
  if [ $INCRLIST == "YES" ]
  then
      echo "LISTA contains $NPDB filenames (before SIEVE)"
      echo "but only the last ${NUMPDB_INCR} files will be processed as requested."
      tail -${NUMPDB_INCR} LISTA > tmp; mv tmp LISTA
  fi


  if [ ${SIEVE} -gt 1 ]
  then
      rm -f tmp_LISTA
      for iline in $(seq 1 $SIEVE $NPDB)
      do
         sed -n "${iline},${iline}p" LISTA >> tmp_LISTA
      done
      mv -f tmp_LISTA LISTA
  fi

  if [  ${NCNTION} -gt 0 ]; then mv ../../${TOPOLOGY} . ; mv ../../tmp_CNTION.info CNTION.info ; fi

  for ((I=0;I<=NCNTION-1;I++))
  do
  ion=${CNTION["$I"]}
    if [ ${ion} == 'Na+' ]
  then
          ionup="INA"
  elif [ ${ion} == 'K+' ]
  then
          ionup="IK"
  elif [ ${ion} == 'MG' ]
  then
          ion="MG "
          ionup="IXG"
  fi

  for ((JFRAG=1;JFRAG<=NFRAG_NOION;JFRAG++))
  do
     if [ ${NCNTION_TOT["$I"]} -gt 0 ]  &&  [ ! -e ../${SNAPSHOTS_DIR}/${ionup}_SORTED_${JFRAG}.INFO  ]  \
                                        &&  [ ! -e ../${SNAPSHOTS_DIR}/${ionup/I/}_SORTED_${JFRAG}.INFO  ]
     then
         echo " ../${SNAPSHOTS_DIR}/${ionup}_SORTED_${JFRAG}.INFO does not exist!"
         echo " ../${SNAPSHOTS_DIR}/${ionup/I/}_SORTED_${JFRAG}.INFO does not exist!"
         exit
     fi
     if [ ${NCNTION_TOT["$I"]} -gt 0 ]
     then
            if [ -e ../${SNAPSHOTS_DIR}/${ionup}_SORTED_${JFRAG}.INFO ]
            then
                ln -s ../${SNAPSHOTS_DIR}/${ionup}_SORTED_${JFRAG}.INFO ${ionup}_SORTED_${JFRAG}.INFO
            else
                ln -s ../${SNAPSHOTS_DIR}/${ionup/I/}_SORTED_${JFRAG}.INFO ${ionup}_SORTED_${JFRAG}.INFO
            fi
    fi
  done


  done

  cp $APTAMD/ENTROPY/run_nmode.sh run_nmode.sh

  sed -i "s/DUMMY_CMPLX_MASK/:1-${NRES}/" run_nmode.sh
  sed -i "s/DUMMY_SLVNTMASK/${SLVNTMASK}/" run_nmode.sh


  if [ "$NFRAG" -eq 1 ]  &&  [ ${NCNTION} -gt 0 ] &&  [ ${CNTION_FRAG} != "YES" ]
  then

     txt="01"
     sed -i "s/DUMMY_NFRAG/${NFRAG}/"  run_nmode.sh
     sed -i "s/DUMMY_FRAG_${txt}/:1-${NRES}/" run_nmode.sh
     j=${FRAG_NWAT[1]} 
     mask=${FRAG_WATMASK[1]} 
     sed -i "s/DUMMY_NWAT_FRAG_${txt}/${j}/" run_nmode.sh
     sed -i "s/DUMMY_WAT_FRAG_${txt}/${mask}/" run_nmode.sh


  elif [ "$NFRAG" -gt 1 ] 
  then
     sed -i "s/DUMMY_NFRAG/${NFRAG}/"  run_nmode.sh
     if [ $NFRAG -gt 20 ]; then echo "Too many frags. Adapt run_nmode.sh!"; exit; fi

     NCNTION_ACCUM=0
     for ((I=0;I<=NCNTION-1;I++))
     do

     for ((i=0;i<=NFRAG-1;i++))
     do
         let "ifrag=$i+1"
         MASK=":${IRES["$i"]}-${JRES["$i"]}"
         if [ ${ifrag} -lt 10 ]; then txt="0${ifrag}"; else  txt=$ifrag; fi

         if [ ${NCNTION_TOT["$I"]} -gt 0 ] &&  [ ${CNTION_FRAG} != "YES" ]
         then
            if [ ${NCNTION_FRAG["$I"]} -eq 1 ] && [ $ifrag -gt 1 ]
            then
               NCNTION_FRAG_SPLIT=0
            else
               NCNTION_FRAG_SPLIT=$(echo ${NCNTION_LIMIT["$I"]} | sed 's/:/\n/g' | head -${ifrag} | tail -1)
            fi
            if [ $NCNTION_FRAG_SPLIT -gt 0 ]
            then 
               for ((j=1;j<=NCNTION_FRAG_SPLIT;j++))
               do  
                    let "k=${NRES_NOION}+$j+${NCNTION_ACCUM}"
                    MASK="${MASK},${k}"
               done
               j=${FRAG_NWAT["$ifrag"]} 
               mask=${FRAG_WATMASK["$ifrag"]} 
               sed -i "s/DUMMY_NWAT_FRAG_${txt}/${j}/" run_nmode.sh
               sed -i "s/DUMMY_WAT_FRAG_${txt}/${mask}/" run_nmode.sh
               let "NCNTION_ACCUM=${NCNTION_ACCUM}+${NCNTION_FRAG_SPLIT}"
            fi
         fi
         sed -i "s/DUMMY_FRAG_${txt}/${MASK}/" run_nmode.sh
     done

     done
  elif [ "$NFRAG" -eq 1 ]
  then
     sed -i "s/DUMMY_NFRAG/0/"  run_nmode.sh
  fi

  if [ ${NCNTION} -gt 0 ]  &&  [ ${CNTION_FRAG} == "YES" ]
  then 
     if [ $NFRAG -gt 20 ]; then echo "Too many frags. Adapt run_mmpbsa.sh!"; exit; fi
     for ((i=0;i<=NFRAG-1;i++))
     do
         let "ifrag=$i+1"
         if [ ${ifrag} -lt 10 ]; then txt="0${ifrag}"; else  txt=$ifrag; fi
         j=${FRAG_NWAT["$i"]} 
         mask=${FRAG_WATMASK["$i"]} 
         sed -i "s/DUMMY_NWAT_FRAG_${txt}/${j}/" run_nmode.sh
         sed -i "s/DUMMY_WAT_FRAG_${txt}/${mask}/" run_nmode.sh
     done

  fi

  sed -i "s/DUMMY_FRAG_..//" run_nmode.sh
  sed -i "s/DUMMY_WAT_FRAG_..//" run_nmode.sh
  sed -i "s/DUMMY_NWAT_FRAG_..//" run_nmode.sh

  sed -i 's/DUMMY_DO_PEEL/YES/' run_nmode.sh
  sed -i 's/DUMMY_PEEL/8.0/'    run_nmode.sh
  sed -i "s/DUMMY_TEMPERATURE/${TEMPERATURE}/"    run_nmode.sh
  TMP_PREPARE="${APTAMD}/MMPBSA/prepare_snap.sh ${WORKDIR_TRJ}/${MOL}_${MD_TYPE}/6.ANALYSIS/${ENTRO_DIR}/CNTION.info "
  TMP_PREPARE=${TMP_PREPARE//\//\\\/}
  sed -i 's/# PREPARE_SNAP=/PREPARE_SNAP=/' run_nmode.sh
  sed -i "s/DUMMY_PREPARE/${TMP_PREPARE}/" run_nmode.sh
    
  if [ $NCNTION -gt 0 ]
  then
      TMP_TOPOLOGY="${WORKDIR_TRJ}/${MOL}_${MD_TYPE}/6.ANALYSIS/${ENTRO_DIR}/${TOPOLOGY}"
  else
      TMP_TOPOLOGY="${WORKDIR_TRJ}/${MOL}_${MD_TYPE}/${TOPOLOGY}"
  fi 

  TMP_TOPOLOGY=${TMP_TOPOLOGY//\//\\\/}
  sed -i "s/DUMMY_TOPOLOGY/${TMP_TOPOLOGY}/" run_nmode.sh
  TMP_SNAP="${WORKDIR_TRJ}/${MOL}_${MD_TYPE}/6.ANALYSIS/${SNAPSHOTS_DIR}/"
  TMP_SNAP=${TMP_SNAP//\//\\\/}

  sed -i "s/DUMMY_SNAPSHOTS/${TMP_SNAP}/" run_nmode.sh
  sed -i "s/DUMMY_ISTRNG/${ISTRNG}/"  run_nmode.sh

  INPUT_LIST="LISTA"
  NINPUT=$(cat $INPUT_LIST | wc -l)
  let " NSPLIT =  ( $NINPUT / $NPROCS ) + 1 "

  TASK="./run_nmode.sh"
  chmod 755 $TASK
  TT=$(date +%N)
    
  rm -f TASK.sh
  rm -f temp_task_list_*
    
  split -l $NSPLIT -d $INPUT_LIST temp_task_list_
    
  for file in $(ls temp_task_list_*)
  do
    echo " $TASK  $file >  $SCRATCH/${file}_${TT}.log " >> TASK.sh
  done
    
  if [ "${FINERSIEVE}" == "YES" ]
  then
        echo "Previous output files (if available) are unpacked" 
        echo "so that this new execution does not repeat calcs"
        if [ -e  OUTPUT.tar ]; then cp OUTPUT.tar OUTPUT_prev.tar ; tar xvf OUTPUT.tar; rm -f OUTPUT.tar  ; fi
  fi

  echo "Running parallel $TASK across  $NPROCS  procs ..."
  cat TASK.sh  | $PARHOME/bin/parallel --no-notice  -t -j$NPROCS
    
  rm -f temp_task_list_* $SCRATCH/temp_task_list*_${TT}.log 

# Merging results 
  if [ "${INCRLIST}" == "YES" ]
  then
        echo "Merging new calcs with previous output files (if available)" 
        if [ -e  OUTPUT.tar ]; then cp OUTPUT.tar OUTPUT_prev.tar ; tar xvf OUTPUT.tar; rm -f OUTPUT.tar  ; fi
        ls *.nmode.gz | sed 's/nmode.gz/pdb/'  > LISTA
  fi
  if  [ "${FINERSIEVE}" == "YES" ] ; then ls *.nmode.gz | sed 's/nmode.gz/pdb/'  > LISTA;  fi
    
# Incorporating qRRHO approx.
  echo  '# PREFIX ENTRO  qRRHO_corr  ENTRO_qRRHO' > ENTRO_qRRHO.dat
  for file in $(ls *.nmode.gz)
  do
      entro=$(zcat ${file/.gz/} | grep 'Total            ' | awk '{print $4}')
      zcat ${file/.gz/} | sed '1,/vibrational     /d' | grep -v cpu | sed '/^$/d'  > $SCRATCH/temp.dat
      nmodes=$(cat $SCRATCH/temp.dat| wc -l)
      echo $nmodes > $SCRATCH/temp_S.dat
      awk '{print $2}' $SCRATCH/temp.dat >> $SCRATCH/temp_S.dat
      $APTAMD/AUXTOOLS/qRRHO < $SCRATCH/temp_S.dat > $SCRATCH/temp_qRRHO.dat
      corr=$(grep CORR  $SCRATCH/temp_qRRHO.dat | awk '{print $2}')
      entro_corr=$(echo $entro  $corr | awk '{print $1+$2}')
      echo " ${file/.nmode.gz/}  $entro   $corr  $entro_corr " >> ENTRO_qRRHO.dat
  done
  rm -f  $SCRATCH/temp_qRRHO.dat  $SCRATCH/temp_S.dat  $SCRATCH/temp.dat

  # Getting data for fragments
  if [ $NFRAG -gt 1 ]
  then
  for ((i=1;i<=$NFRAG;i++))
  do
	  let "j=$i+3"
  	  echo  '# PREFIX ENTRO  qRRHO_corr  ENTRO_qRRHO' > ENTRO_qRRHO_${i}.dat
	  for file in $(ls *.nmode${j}.gz)
 	  do
  	      entro=$(zcat ${file/.gz/} | grep 'Total            ' | awk '{print $4}')
	      zcat ${file/.gz/} | sed '1,/vibrational     /d' | grep -v cpu | sed '/^$/d'  > $SCRATCH/temp.dat
              nmodes=$(cat $SCRATCH/temp.dat| wc -l)
              echo $nmodes > $SCRATCH/temp_S.dat
              awk '{print $2}' $SCRATCH/temp.dat >> $SCRATCH/temp_S.dat
              $APTAMD/AUXTOOLS/qRRHO < $SCRATCH/temp_S.dat > $SCRATCH/temp_qRRHO.dat
              corr=$(grep CORR  $SCRATCH/temp_qRRHO.dat | awk '{print $2}')
              entro_corr=$(echo $entro  $corr | awk '{print $1+$2}')
              echo " ${file/.nmode${i}.gz/}  $entro   $corr  $entro_corr " >> ENTRO_qRRHO_${i}.dat
          done
          rm -f  $SCRATCH/temp_qRRHO.dat  $SCRATCH/temp_S.dat  $SCRATCH/temp.dat
   done

   # Getting differences
   i=0
   for file in $(ls ENTRO_qRRHO_*.dat)
   do
	   let "i=$i+1"
 	   if [ $i -eq 1 ]
	   then
		   grep -v '#' $file > temp_sum.dat
	   else
		   grep -v '#' $file > temp_frag.dat
		   paste -d '  '  temp_sum.dat temp_frag.dat | awk '{printf("%s  %f  %f  %f \n",$1,$2+$6,$3+$7,$4+$8)}'   > temp_sum_frag; mv -f temp_sum_frag temp_sum
	   fi
   done
   echo '# PREFIX ENTRO_diff  qRRHO_corr_diff  ENTRO_qRRHO_diff ' > ENTRO_qRRHO_DIFF.dat
   grep -v '#' ENTRO_qRRHO.dat > temp_frag.dat
   paste -d  ' '  temp_frag.dat temp_sum.dat  | awk '{printf("%s  %f  %f  %f \n",$1,$2-$6,$3-$7,$4-$8)}' >> ENTRO_qRRHO_DIFF.dat
   rm -f temp_sum.dat temp_frag.dat

  fi

  # Averaging
  export PERCEN="0"
  export OUTLYER="0"
  export DO_ENTRO_NMODE_TAR="1"
  export DO_STAT_PLOT="1"
  export DO_STAT_PERCEN=0
  $APTAMD/ENTROPY/nmode_data_parser.sh LISTA

  # T_weighted differences
  if [ $NFRAG -gt 1 ]
  then
	  echo '# PREFIX -T*ENTRO_diff  -T*qRRHO_corr_diff  -T*ENTRO_qRRHO_diff ' > T_ENTRO_qRRHO_DIFF.dat
	  grep -v '#' ENTRO_qRRHO_diff.dat | awk  -v T=$TEMPERATURE '{printf("%s  %f  %f  %f \n",$1,-T*$2/1000.0, -T*$3/1000.0, -T*$4/1000.0)}'  >> T_ENTRO_qRRHO_DIFF.dat
	  echo '# PREFIX -T*ENTRO_diff  -T*ENTRO_TRANS_diff  -T*ENTRO_ROT_diff -T*ENTRO_VIB_diff ' > T_ENTRO_COMP_DIFF.dat
	  grep -v '#' ENTRO_COMP_DIFF.dat | awk  -v T=$TEMPERATURE '{printf("%s  %f  %f  %f  %f \n",$1,-T*$2/1000.0, -T*$3/1000.0, -T*$4/1000.0, -T*$5/1000.0)}'  >> T_ENTRO_COMP_DIFF.dat
  fi

  # Statistics 
  $APTAMD/MMPBSA/stat_plot.sh

  cd $WORKDIR_TRJ 

done


