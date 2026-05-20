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
if [ -z "$NSODIUM_LIMIT" ]
then
   NSODIUM_LIMIT=0
else
   echo "Including $NSODIUM_LIMIT Na+ ions (plus Oh hydration shells)"
fi

if [ -n "$PBS_ENVIRONMENT" ] ; then
  NPROCS=$(cat $PBS_NODEFILE | wc -l)
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
   echo "Processing 1 out of $SIEVE snapshots"
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

# SODIUM ions can be separated fragments or just aggregated to FRAG=1 
if [ -z $SODIUM_FRAG ]
then 
   SODIUM_FRAG="NO"
fi

if [ -z "$ENTRO_DIR" ]; then ENTRO_DIR="ENTRO_NMODE"; fi

WORKDIR_TRJ=$PWD
# NSODIUM data
declare -a NSODIUM_LIST=""
NSODIUM_LIST=($(echo $NSODIUM_LIMIT))
NSODIUM_TERMS=${#NSODIUM_LIST[@]}
NSODIUM=0
for ((i=0;i<=NSODIUM_TERMS-1;i++))
do
    let "NSODIUM=${NSODIUM_LIST["$i"]}+$NSODIUM"
done
echo "Considering a total of $NSODIUM ions"

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
    if [ $NFRAG -ne $NSODIUM_TERMS ]
    then
       echo "Detected $NFRAG fragments in $MOL, but NSODIUM is specified for $NSODIUM_TERMS fragments"
       if [ $NSODIUM_TERMS -eq 1 ]
       then
           echo "Assuming that NSODIUM=${NSODIUM_LIST["0"]} applies to FRAG 1"
           for ((IFRAG=1;IFRAG<=NFRAG-1;IFRAG++))
           do
               NSODIUM_LIST["$IFRAG"]=0
           done
           NSODIUM_TERMS=$NFRAG
      else
           exit
      fi
    fi
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

  if [ ${NSODIUM} -gt 0  ]
  then
#      We process the TOPOLOGY file to modify the atom/resname for Na+
#      Na+---> INA ( INA atoms will be kept in ENTRO_NMODE analyses)
       echo  "NSODIUM=${NSODIUM} to be considered as solute atoms"
       cp $TOPOLOGY ${MOL}_INA.top
       TOPOLOGY="${MOL}_INA.top"

       NTOP_SODIUM=$(sed  '/FLAG CHARGE/,$d' $TOPOLOGY | tr " \t" "\n" | grep -c 'Na+')

       if [ ${NTOP_SODIUM} -lt ${NSODIUM} ] 
       then
          echo "Not enough Na+ in $TOPOLOGY  NSODIUM=${NSODIUM}  NTOP_SODIUM=${NTOP_SODIUM}"
          exit
       fi

       sed -n '1,/FLAG CHARGE/p' $TOPOLOGY > temp_A.top
       sed -n '/FLAG CHARGE/,/FLAG AMBER_ATOM_TYPE/p' $TOPOLOGY  | sed '1,1d' > temp_B.top
       sed  '1,/FLAG AMBER_ATOM_TYPE/d' $TOPOLOGY > temp_C.top

       for file in temp_A.top temp_B.top 
       do
            declare -a iline=""
            iline=($(grep -n 'Na+' $file | sed 's/:/  /' | awk '{print $1}'))
            nlines=${#iline[@]} 
            NUM_INA=0
            for ((i=1;i<=nlines;i++))
            do
               let "j=$i-1"
               jline=${iline["$j"]}
               NSODIUM_LINE=$(sed -n "${jline},${jline}p"  $file | tr " \t" "\n" | grep -c 'Na+') 
               if [ ${NUM_INA} -lt ${NSODIUM} ]
               then
                   let "NUM_INA_NEED=${NSODIUM}-${NUM_INA}"
                   if [ ${NSODIUM_LINE} -le ${NUM_INA_NEED} ]
                   then 
                       sed -i "${jline}s/Na+/INA/g" $file
                       let "NUM_INA=${NUM_INA} + ${NSODIUM_LINE}"
                   else
                       for ((k=1;k<=NUM_INA_NEED;k++)); do sed -i "${jline}s/Na+/INA/" $file; done
                       let "NUM_INA=${NUM_INA} + ${NUM_INA_NEED}"
                   fi
               else
                   break
               fi
             done
       done
       cat temp_A.top temp_B.top temp_C.top > $TOPOLOGY
       rm -f temp_A.top temp_B.top temp_C.top

#      Fixing some variables according to the number of sodium ions
       ENTRO_DIR="ENTRO_NMODE"

       if [ $NSODIUM -eq ${NSODIUM_LIST[0]} ]
       then
           ENTRO_DIR=${ENTRO_DIR}_${NSODIUM}
       else
           for ((IFRAG=0;IFRAG<=NFRAG-1;IFRAG++))
           do
              ENTRO_DIR=${ENTRO_DIR}_${NSODIUM_LIST["$IFRAG"]}
           done
       fi

       declare -a FRAG_NWAT=""
       declare -a FRAG_WATMASK=""

       if [ ${SODIUM_FRAG} == "YES" ]
       then

          for ((i=0;i<=NFRAG-1;i++))
          do
             FRAG_NWAT["$i"]=0
             FRAG_WATMASK["$i"]=""
          done
          for ((i=1;i<=NSODIUM;i++))
          do
             let " j=${NRES} + $i "
             let " k=${NFRAG} -1 + $i "
             IRES["$k"]=$j
             JRES["$k"]=$j
             FRAG_NWAT["$k"]=6
             FRAG_WATMASK["$k"]=":${j}"
          done
          let "NFRAG=${NFRAG}+${NSODIUM}"
          NRES_NOION=${NRES} 
          let "NRES=${NRES}+${NSODIUM}"
     
       else

        NSODIUM_ACCUM=0
        for ((IFRAG=0;IFRAG<=NFRAG-1;IFRAG++))
        do

         let "JFRAG=$IFRAG+1"
         NSODIUM_FRAG=${NSODIUM_LIST["$IFRAG"]}

         if [ $NSODIUM_FRAG -gt 0 ]
         then
            txt=" 6 "
            let " j=${NRES} + 1 + ${NSODIUM_ACCUM}  "
            watmask=" :${j} "
            for ((i=2;i<=NSODIUM_FRAG;i++))
            do
               let " j=${NRES} + $i + ${NSODIUM_ACCUM}"
               txt="${txt} ; 6 "
               watmask="${watmask} ; :${j} "
            done
            FRAG_NWAT["$JFRAG"]=${txt}
            FRAG_WATMASK["$JFRAG"]=${watmask}
            let "NSODIUM_ACCUM=${NSODIUM_ACCUM}+${NSODIUM_FRAG}"
         else
            FRAG_NWAT["$JFRAG"]=""
            FRAG_WATMASK["$JFRAG"]=""
         fi

        done

         NRES_NOION=${NRES} 
         let "NRES=${NRES}+${NSODIUM}"

       fi
       
  fi
       
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
  if [  ${NSODIUM} -gt 0 ]; then mv ../../${TOPOLOGY} . ; fi

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

  for ((JFRAG=1;JFRAG<=NFRAG;JFRAG++))
  do
     if [ ${NSODIUM} -gt 0 ]  &&  [ ! -e ../${SNAPSHOTS_DIR}/NA_SORTED_${JFRAG}.INFO  ]
     then
         echo " ../${SNAPSHOTS_DIR}/NA_SORTED_${JFRAG}.INFO does not exist!"
         exit
     fi
     if [ ${NSODIUM} -gt 0 ]; then ln -s ../${SNAPSHOTS_DIR}/NA_SORTED_${JFRAG}.INFO NA_SORTED_${JFRAG}.INFO; fi
  done

  cp $APTAMD/ENTROPY/run_nmode.sh run_nmode.sh

  sed -i "s/DUMMY_CMPLX_MASK/:1-${NRES}/" run_nmode.sh
  if [ "$NFRAG" -eq 1 ]  &&  [ ${NSODIUM} -gt 0 ] &&  [ ${SODIUM_FRAG} != "YES" ]
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
     NSODIUM_ACCUM=0
     for ((i=0;i<=NFRAG-1;i++))
     do
         let "ifrag=$i+1"
         MASK=":${IRES["$i"]}-${JRES["$i"]}"
         if [ ${ifrag} -lt 10 ]; then txt="0${ifrag}"; else  txt=$ifrag; fi
         if [ ${NSODIUM} -gt 0 ] &&  [ ${SODIUM_FRAG} != "YES" ]
         then
            NSODIUM_FRAG=${NSODIUM_LIST["$i"]}
            if [ $NSODIUM_FRAG -gt 0 ]
            then
               for ((j=1;j<=NSODIUM_FRAG;j++))
               do
                    let "k=${NRES_NOION}+$j+${NSODIUM_ACCUM}"
                    MASK="${MASK},${k}"
               done
               j=${FRAG_NWAT["$ifrag"]}
               mask=${FRAG_WATMASK["$ifrag"]}
               sed -i "s/DUMMY_NWAT_FRAG_${txt}/${j}/" run_nmode.sh
               sed -i "s/DUMMY_WAT_FRAG_${txt}/${mask}/" run_nmode.sh
               let "NSODIUM_ACCUM=${NSODIUM_ACCUM}+${NSODIUM_FRAG}"
            fi
         fi
         sed -i "s/DUMMY_FRAG_${txt}/${MASK}/" run_nmode.sh
     done
  elif [ "$NFRAG" -eq 1 ]
  then 
     sed -i "s/DUMMY_NFRAG/${NFRAG}/"  run_nmode.sh
  fi
  sed -i "s/DUMMY_FRAG_..//" run_nmode.sh
  sed -i "s/DUMMY_WAT_FRAG_..//" run_nmode.sh
  sed -i "s/DUMMY_NWAT_FRAG_..//" run_nmode.sh

  if [ ${NSODIUM} -gt 0 ]  &&  [ ${SODIUM_FRAG} == "YES" ]
  then 
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
  sed -i 's/DUMMY_DO_PEEL/YES/' run_nmode.sh
  sed -i 's/DUMMY_PEEL/8.0/' run_nmode.sh

  TMP_PREPARE="${APTAMD}/MMPBSA/prepare_snap.sh ${NSODIUM_LIMIT} "
  TMP_PREPARE=${TMP_PREPARE//\//\\\/}
  echo "${TMP_PREPARE}"
  sed -i 's/# PREPARE_SNAP=/PREPARE_SNAP=/'  run_nmode.sh
  sed -i "s/DUMMY_PREPARE/${TMP_PREPARE}/"   run_nmode.sh

  if [ $NSODIUM -eq 0 ]
  then 
     TMP_TOPOLOGY="${WORKDIR_TRJ}/${MOL}_${MD_TYPE}/${TOPOLOGY}"
  else
     TMP_TOPOLOGY="${WORKDIR_TRJ}/${MOL}_${MD_TYPE}/6.ANALYSIS/${ENTRO_DIR}/${TOPOLOGY}"
  fi
  TMP_TOPOLOGY=${TMP_TOPOLOGY//\//\\\/}
  echo $TMP_TOPOLOGY
  sed -i "s/DUMMY_TOPOLOGY/${TMP_TOPOLOGY}/" run_nmode.sh
  TMP_SNAP="${WORKDIR_TRJ}/${MOL}_${MD_TYPE}/6.ANALYSIS/${SNAPSHOTS_DIR}/"
  TMP_SNAP=${TMP_SNAP//\//\\\/}
  echo $TMP_SNAP
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
  echo  '# PREFIX S_RRHO  qRRHO_corr  S_qRRHO' > S_qRRHO.dat
  for file in $(ls *.nmode.gz)
  do
      echo $file
      entro=$(zcat ${file/.gz/} | grep 'Total            ' | awk '{print $4}')
      zcat ${file/.gz/} | sed '1,/vibrational     /d' | grep -v cpu | sed '/^$/d'  > $SCRATCH/temp.dat
      nmodes=$(cat $SCRATCH/temp.dat| wc -l)
      echo $nmodes > $SCRATCH/temp_S.dat
      awk '{print $2}' $SCRATCH/temp.dat >> $SCRATCH/temp_S.dat
      $APTAMD/AUXTOOLS/qRRHO < $SCRATCH/temp_S.dat > $SCRATCH/temp_qRRHO.dat
      corr=$(grep CORR  $SCRATCH/temp_qRRHO.dat | awk '{print $2}')
      entro_corr=$(echo $entro  $corr | awk '{print $1+$2}')
      echo " ${file/.nmode.gz/}  $entro   $corr  $entro_corr " >> S_qRRHO.dat
  done
  rm -f  $SCRATCH/temp_qRRHO.dat  $SCRATCH/temp_S.dat  $SCRATCH/temp.dat

# Averaging  
  export PERCEN="0"
  export OUTLYER="0"
  export DO_ENTRO_NMODE_TAR="1"
  export PERCEN="0"
  export DO_STAT_PLOT="1"
  export DO_STAT_PERCEN=0
  $APTAMD/ENTROPY/nmode_data_parser.sh LISTA 
  $APTAMD/MMPBSA/stat_plot.sh

  cd $WORKDIR_TRJ 

done


