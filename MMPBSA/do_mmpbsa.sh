#!/bin/bash

# This is the main script in the MMPBSA (and QM/MMPBSA) mess performed
# by the APTAMD suite. 

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
if [ "$#" -eq 0 ]; then more $APTAMD/DOC/do_mmpbsa.txt ; exit; fi

source $APTAMD/ENV/aptamd_env.sh

# Checking options
extension="${1##*.}"
if [ $extension == "src" ]
then
    echo "Sourcing $1"
    source $1
else
 # Put into double quotes the list of trajectories to be processed
 MD_TRAJ=$1
 #  GAMD or MD trajectory type should be specified 
 MD_TYPE=$2
fi

if [ ! -z "$MOL" ] && [  -z "$MD_TRAJ" ]; then MD_TRAJ="${MOL}"; echo $MD_TRAJ; unset MOL; fi

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

EXIT="NO"
if [ -z "$MD_TRAJ" ]; then  EXIT="YES"; fi
if [ -z "$MD_TYPE" ]; then  EXIT="YES" ; fi

if [ "$EXIT" == "YES" ]
then
   echo 'Usage: do_mmpbsa.sh "MOL1   MOL2 .."   [ MD | GAMD ]  "NSODIUM_FRAG_1  NSODIUM_FRAG_2 .." '
   echo
   echo 'Many, many other options about activating RISM, QMMM, ADCK, are controlled by defining' 
   echo 'external BASH variables. '
   echo 'See $APTAMD/DOC/do_mmpbsa.txt and $APTAMD/DOC/energy_analysis.txt for details.'
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
if [ ! -z $TOPOLOGY ]
then
  echo  "Using TOPOLOGY=$TOPOLOGY for ALL systems in ${MD_TRAJ}"
fi 

if [ -z "$NSODIUM_LIMIT" ] 
then 
  NSODIUM_LIMIT=0 ; 
  echo "Using NSODIUM_LIMIT=${NSODIUM_LIMIT}"
else
  echo "Using NSODIUM_LIMIT=${NSODIUM_LIMIT} as predefined"
fi
# SODIUM ions can be separated fragments or just aggregated to cmplx/receptor/ligand fragments
if [ -z "$SODIUM_FRAG" ]; then SODIUM_FRAG="NO"; fi
if [ "$SODIUM_FRAG" == "NO" ] && [ $NSODIUM_LIMIT  -gt  0 ]
then 
  echo "Na+ and hydrating waters are assigned to CMPLX/FRAGs"
elif [ "$SODIUM_FRAG" == "YES" ] && [ $NSODIUM_LIMIT  -gt  0 ] 
then
  echo "Na+ and hydrating waters are treated as additional fragments."
fi

# SOME  MMPBSA parameters 
if [ -z "$PDIE_LIST" ]
then
   PDIE_LIST="1" 
   echo "Using PDIE=${PDIE_LIST}"
else
   echo "Using PDIE=${PDIE_LIST} as predefined "
fi
if [ -z "$ISTRNG" ]
then
   ISTRNG="150"
   echo "Using ISTRNG=${ISTRNG}"
else
   echo "Using ISTRNG=${ISTRNG} as predefined "
fi
if [ -z "$PEEL" ]
then
   PEEL="12.0"
   echo "Using PEEL=${PEEL}"
else
   echo "Using PEEL=${PEEL} as predefined "
fi
if [ -z "$XVVFILE" ] 
then
   RISM="NO"
elif [ -e "$XVVFILE" ]
then
   RISM="YES"
   echo "XVVFILE=$XVVFILE"
   echo "Activating RISM option"
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
   SIEVE="1"
else 
   echo "Processing only a fraction of PDB files"
   echo "SIEVE=$SIEVE"
fi
if [ -z "$SKIPFRAG" ] 
then
   SKIPFRAG="NO"
elif [ "$SKIPFRAG" == "YES" ]
then
   echo "SKIP fragment calcs in composed systems."
fi
if [ -z "$SKIPCMPLX" ] 
then
   SKIPCMPLX="NO"
elif [ "$SKIPCMPLX" == "YES" ]
then
   echo "SKIP complex calcs in composed systems."
fi

if [ -z "$USEOUTPUT" ] 
then
   USEOUTPUT="NO"
elif [ "$USEOUTPUT" == "YES" ]
then
   if [ -z "$USEOUTPUT_MASK" ]
   then
     echo "USEOUTPUT=YES but USEOUTPUT_MASK not defined!" 
     stop
   fi
   echo "Using available OUTPUT.tar with MASK=$USEOUTPUT_MASK"
fi

# QMMM stuff 
if [ -z $QM_CMPLX_MASK ] && [ -z $QM_MASK ]
then
   QMMM="NO"
   QM_ORCA="NO"
else
   QMMM="YES"
   if [ -z $QM_CMPLX_MASK ]; then QM_CMPLX_MASK=$QM_MASK ; fi
   echo "QM_CMPLX_MASK=$QM_CMPLX_MASK QM/MM PBSA calcs will be done"
fi
if [ $QM_CHARGE ]  &&  [  $QMMM == "YES" ]
then 
   QM_CMPLX_CHARGE=$QM_CHARGE
fi
if [ $QM_CMPLX_CHARGE ]  &&  [  $QMMM == "YES" ]
then 
   echo "Reading QM_CMPLX_CHARGE=$QM_CMPLX_CHARGE "
fi
#  ORCA settings for QM/MM 
if [ -z $QM_ORCA ] && [ $QMMM == "YES" ]
then 
    QM_ORCA="NO"
    echo "QMMM calcs requested, using DFTB3 Hamiltonian"
fi   
if [ $QM_ORCA == "YES" ] && [ -z $QM_ORCA_LEVEL ]
then 
    echo "QMMM calcs requested using ORCA interface"
    echo "but QM_ORCA_LEVEL is not defined !"
    exit
elif [ $QM_ORCA == "YES" ]
then
    echo "QMMM calcs requested using ORCA interface"
    echo "QM_ORCA_LEVEL= ${QM_ORCA_LEVEL}"
fi   
if [ -z $NPROCS_QM ] && [ $QMMM == "YES" ] 
then 
    echo "QMMM calcs requested, but NPROCS_QM is not defined !"
    echo "Using NPROCS_QM=1"
    NPROCS_QM=1
elif [ $QMMM == "YES" ]
then
    echo "NPROCS_QM = ${NPROCS_QM}"
else 
   NPROCS_QM=1
fi   
# ORCA settings for QM/MM relaxation
if [ -z $QM_ORCA_RELAX ] 
then 
    QM_ORCA_RELAX="NO"
fi   
if [ -z $QM_ORCA_RELAX_LEVEL ] && [ $QM_ORCA_RELAX == "YES" ]
then 
    echo "QMMM relaxation calcs requested using ORCA interface"
    echo "but QM_ORCA_RELAX_LEVEL is not defined !"
    exit
elif [ $QM_ORCA_RELAX == "YES" ]
then
    echo "QMMM relaxation calcs requested using ORCA interface"
    echo "QM_ORCA_RELAX_LEVEL= ${QM_ORCA_RELAX_LEVEL}"
fi   
# TC settings for QM/MMM
if [ -z $QM_TC_RELAX ]  
then 
    QM_TC_RELAX="NO"
fi   
if [ -z $QM_TC_LEVEL ] && [ $QM_TC_RELAX == "YES" ]
then 
    echo "QMMM relaxation calcs requested using TC interface"
    echo "but QM_TC_LEVEL is not defined !"
    exit
elif [ $QM_TC_RELAX == "YES" ]
then
    echo "QMMM relaxation calcs requested using TC interface"
    echo "QM_TC_LEVEL= ${QM_TC_LEVEL}"
fi   
if [ -z "$GPUID" ] && [ $QM_TC_RELAX == "YES" ]
then
   GPUID="0"
fi
if [ $QM_TC_RELAX == "YES" ]; then export GPUID; fi

# Definition of QM FRAGs
if [ -z ${QM_FRAG_MASK// /} ] && [ $QMMM == "YES" ]
then 
   echo "QMMM calcs requested, but QM_FRAG_MASK undefined"
   echo "It is assumed that this is a 1 fragment calculation"
   QM_FRAG_MASK="$QM_CMPLX_MASK"
   QM_FRAG_CHARGE="$QM_CMPLX_CHARGE"
fi
if [ -z ${QM_FRAG_CHARGE// /} ] && [ $QMMM == "YES" ]
then 
   echo "QMMM calcs requested, but QM_CHARGE_MASK undefined"
   exit
fi
if [ -z ${QM_ADD_WAT_FRAG// /} ] && [ $QMMM == "YES" ]
then 
   echo "QMMM calcs requested, but QM_ADD_WAT_FRAG undefined"
   echo "Assuming that no Na+/WAT is added to the QM region" 
   QMMM_ADD_SOLV="NO"
elif [ ! -z ${QM_ADD_WAT_FRAG// /} ] && [ $QMMM == "YES" ]
then
   echo "Assuming that Na+/WAT may be added to the QM region" 
   echo "QM_ADD_WAT_FRAG=$QM_ADD_WAT_FRAG"
   QMMM_ADD_SOLV="YES"
fi

# For DFTB3 calcs, D3H4 empirical corrections calculated with
# CUBY4 can be incorporated
if [ -z $D3H4 ]
then
  D3H4="NO"
fi 
if [ $QMMM == "YES" ] && [ $QM_ORCA == "YES" ] 
then
   D3H4="NO"
elif [ $QMMM == "YES" ] 
then
   if [ $D3H4 == "YES" ]
   then 
      echo "D3H4 corrections to SCC-DFTB3 energies requested using CUBY4 "
      if [ -z $CUBY4 ] 
      then
           D3H4="NO"
           echo "However, the CUBY4 environmental variable is not defined"
           echo "in aptamd_env.sh so that D3H4 calcs can not be done."
      fi
   fi
fi

# ADCK options 
if [ -z $ADCK ]
then
   ADCK="NO"
elif [ $ADCK == "YES" ]
then
   echo "ADCK=YES activating Autodock Scoring"
   echo 'A PDBQT_SRCFILE is required. Each line defines the corresponding PDBQT_FRAG variable:'
   echo 'PDBQT_FRAG[1]= .... '
   echo 'PDBQT_FRAG[2]= .... '
fi
#
if [ -z $ADCK_ONLY ]
then
   ADCK_ONLY="NO"
elif [ $ADCK_ONLY == "YES" ]
then
   echo "ADCK_ONLY=YES no other calcs will be performed"
   echo 'A PDBQT_SRCFILE is required. Each line defines the corresponding PDBQT_FRAG variable:'
   echo 'PDBQT_FRAG[1]= .... '
   echo 'PDBQT_FRAG[2]= .... '
   ADCK="YES"
fi

if [ "$ADCK" == "YES" ] && [ -z $PDBQT_SRCFILE ] 
then
   echo 'ADCK=YES, but PDBQT_SRCFILE is not defined!'
   exit
fi
if [ "$ADCK" == "YES" ] && [ ! -e  $PDBQT_SRCFILE ]
then
   echo 'ADCK=YES, but PDBQT_SRCFILE=$PDBQT_SRCFILE does not exist !'
   exit
fi

# Directory names
if [ -z $MMPBSA_DIR ]
then 
   if [ $RISM == "YES" ]; then MMPBSA_DIR="MMRISM"; else MMPBSA_DIR="MMPBSA"; fi
   if [ "$QMMM" == "YES" ]; then MMPBSA_DIR="QM${MMPBSA_DIR}";fi 
   echo "MMPBSA_DIR=$MMPBSA_DIR default"
else
   echo "MMPBSA_DIR=$MMPBSA_DIR in the input file"
fi
if [ -z $SNAPSHOTS_DIR ]; then  SNAPSHOTS_DIR="SNAPSHOTS"; fi 


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

  if [ -z $TOPOLOGY ] 
  then

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
     $AMBERHOME/bin/cpptraj $SOLUTE_TOP <<EOF > mlog 
trajin $SOLUTE_CRD
trajout $SOLUTE_PDB  pdb pdbatom
go
EOF
    NRES=$(grep 'ATOM  ' $SOLUTE_PDB   | tail -1 | awk '{print $5}')
    NFRAG=$(grep -c 'TER' $SOLUTE_PDB)
    if [ "$SKIPFRAG" == "YES" ] && [ "$NFRAG" -gt 1 ]
    then
        echo "Detected $NFRAG fragments in $MOL, but SKIPFRAG=YES."
    else
        echo "Detected $NFRAG fragments"
    fi
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
    for ((i=0; i<=NFRAG-1;i++))
    do
         echo  "Fragment :${IRES["$i"]}-${JRES["$i"]}"
    done
    cd ../
  fi

  if [ ${NSODIUM} -gt 0  ]
  then
#
#  We process the TOPOLOGY file to modify the atom/resname for Na+
#  Na+---> INA ( INA atoms will be kept in MMPBSA analyses)
#
       echo  "NSODIUM=${NSODIUM} to be considered as solute atoms"
       cp $TOPOLOGY ${MOL}_solutewat_INA.top
       TOPOLOGY="${MOL}_solutewat_INA.top"

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
       
       if [ $NSODIUM -eq ${NSODIUM_LIST[0]} ]
       then
           MMPBSA_DIR=${MMPBSA_DIR}_NA_${NSODIUM}
       else
           MMPBSA_DIR=${MMPBSA_DIR}_NA
           for ((IFRAG=0;IFRAG<=NFRAG-1;IFRAG++))
           do
              MMPBSA_DIR=${MMPBSA_DIR}_${NSODIUM_LIST["$IFRAG"]}
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

  if [ -e  ${MMPBSA_DIR}  ]
  then 
      echo " ${MMPBSA_DIR} found in ${MOL}_MD/6.ANALYSIS"
  else
      mkdir ${MMPBSA_DIR}
  fi 

  cd ${MMPBSA_DIR}
  WORKDIR=$PWD 
  cp ../${SNAPSHOTS_DIR}/LISTA .
  if [  ${NSODIUM} -gt 0 ]; then mv ../../${TOPOLOGY} . ; fi 

  for ((JFRAG=1;JFRAG<=NFRAG;JFRAG++))
  do

     if [ ${NSODIUM} -gt 0 ]  &&  [ ! -e ../${SNAPSHOTS_DIR}/NA_SORTED_${JFRAG}.INFO  ]
     then 
         echo " ../${SNAPSHOTS_DIR}/NA_SORTED_${JFRAG}.INFO does not exist!"
         exit
     fi
    
     if [ ${NSODIUM} -gt 0 ]; then ln -s ../${SNAPSHOTS_DIR}/NA_SORTED_${JFRAG}.INFO NA_SORTED_${JFRAG}.INFO; fi

  done 

  cp $APTAMD/MMPBSA/run_mmpbsa.sh run_mmpbsa.sh

  sed -i "s/DUMMY_CMPLX_MASK/:1-${NRES}/" run_mmpbsa.sh

  if [ "$NFRAG" -eq 1 ]  &&  [ ${NSODIUM} -gt 0 ] &&  [ ${SODIUM_FRAG} != "YES" ]
  then

     txt="01"
     sed -i "s/DUMMY_NFRAG/${NFRAG}/"  run_mmpbsa.sh
     sed -i "s/DUMMY_FRAG_${txt}/:1-${NRES}/" run_mmpbsa.sh
     j=${FRAG_NWAT[1]} 
     mask=${FRAG_WATMASK[1]} 
     sed -i "s/DUMMY_NWAT_FRAG_${txt}/${j}/" run_mmpbsa.sh
     sed -i "s/DUMMY_WAT_FRAG_${txt}/${mask}/" run_mmpbsa.sh


  elif [ "$NFRAG" -gt 1 ] 
  then
     sed -i "s/DUMMY_NFRAG/${NFRAG}/"  run_mmpbsa.sh
     if [ $NFRAG -gt 20 ]; then echo "Too many frags. Adapt run_mmpbsa.sh!"; exit; fi
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
               sed -i "s/DUMMY_NWAT_FRAG_${txt}/${j}/" run_mmpbsa.sh
               sed -i "s/DUMMY_WAT_FRAG_${txt}/${mask}/" run_mmpbsa.sh
               let "NSODIUM_ACCUM=${NSODIUM_ACCUM}+${NSODIUM_FRAG}"
            fi
         fi
         sed -i "s/DUMMY_FRAG_${txt}/${MASK}/" run_mmpbsa.sh
     done
  elif [ "$NFRAG" -eq 1 ]
  then
     sed -i "s/DUMMY_NFRAG/0/"  run_mmpbsa.sh
  fi
  sed -i "s/DUMMY_FRAG_..//" run_mmpbsa.sh
  sed -i "s/DUMMY_WAT_FRAG_..//" run_mmpbsa.sh
  sed -i "s/DUMMY_NWAT_FRAG_..//" run_mmpbsa.sh

  if [ ${NSODIUM} -gt 0 ]  &&  [ ${SODIUM_FRAG} == "YES" ]
  then 
     for ((i=0;i<=NFRAG-1;i++))
     do
         let "ifrag=$i+1"
         if [ ${ifrag} -lt 10 ]; then txt="0${ifrag}"; else  txt=$ifrag; fi
         j=${FRAG_NWAT["$i"]} 
         mask=${FRAG_WATMASK["$i"]} 
         sed -i "s/DUMMY_NWAT_FRAG_${txt}/${j}/" run_mmpbsa.sh
         sed -i "s/DUMMY_WAT_FRAG_${txt}/${mask}/" run_mmpbsa.sh
     done

  fi

  sed -i 's/DUMMY_DO_PEEL/YES/' run_mmpbsa.sh 
  sed -i "s/DUMMY_PEEL/${PEEL}/" run_mmpbsa.sh 
  TMP_PREPARE="${APTAMD}/MMPBSA/prepare_snap.sh ${NSODIUM_LIMIT} "
  TMP_PREPARE=${TMP_PREPARE//\//\\\/}
  sed -i 's/# PREPARE_SNAP=/PREPARE_SNAP=/' run_mmpbsa.sh
  sed -i "s/DUMMY_PREPARE/${TMP_PREPARE}/" run_mmpbsa.sh
    
  if [ $NSODIUM -gt 0 ]
  then
      TMP_TOPOLOGY="${WORKDIR_TRJ}/${MOL}_${MD_TYPE}/6.ANALYSIS/${MMPBSA_DIR}/${TOPOLOGY}"
  else
      TMP_TOPOLOGY="${WORKDIR_TRJ}/${MOL}_${MD_TYPE}/${TOPOLOGY}"
  fi 
  TMP_TOPOLOGY=${TMP_TOPOLOGY//\//\\\/}
  sed -i "s/DUMMY_TOPOLOGY/${TMP_TOPOLOGY}/" run_mmpbsa.sh
  TMP_SNAP="${WORKDIR_TRJ}/${MOL}_${MD_TYPE}/6.ANALYSIS/${SNAPSHOTS_DIR}/"
  TMP_SNAP=${TMP_SNAP//\//\\\/}
  sed -i "s/DUMMY_SNAPSHOTS/${TMP_SNAP}/" run_mmpbsa.sh
  sed -i "s/DUMMY_ISTRNG/${ISTRNG}/"  run_mmpbsa.sh

  if [ $RISM == "YES" ]
  then
       sed -i 's/DUMMY_RISM/YES/' run_mmpbsa.sh
       TMP_XVVFILE=${XVVFILE//\//\\\/}
       sed -i "s/DUMMY_XVVFILE/${TMP_XVVFILE}/"  run_mmpbsa.sh
       # Solvent is just defined in the XVVFILE
       PDIE_LIST="1" 
  else 
       sed -i 's/DUMMY_RISM/NO/' run_mmpbsa.sh
  fi
 
  if [ $INCRLIST == "YES" ]
  then 
      NPDB=$(cat LISTA | wc -l)
      echo "LISTA contains $NPDB filenames (before SIEVE)"
      echo "but only the last ${NUMPDB_INCR} files will be processed as requested."
      tail -${NUMPDB_INCR} LISTA > tmp; mv tmp LISTA
  fi

  if [ "$SIEVE" -gt 1 ]
  then
      npdb=$(cat LISTA | wc -l)
      rm -f tmp_LISTA
      for I in $(seq 1 $SIEVE $npdb )
      do
          sed -n "${I},${I}p" LISTA >> tmp_LISTA
      done
      mv tmp_LISTA LISTA
  fi

  if [ "$SKIPFRAG" == "YES" ]
  then 
     sed -i "s/DUMMY_SKIPFRAG/YES/" run_mmpbsa.sh
  else
     sed -i "s/DUMMY_SKIPFRAG/NO/" run_mmpbsa.sh
  fi
  if [ "$SKIPCMPLX" == "YES" ]
  then 
     sed -i "s/DUMMY_SKIPCMPLX/YES/" run_mmpbsa.sh
  else
     sed -i "s/DUMMY_SKIPCMPLX/NO/" run_mmpbsa.sh
  fi

# Transforming QM values into arrays
  if [  "$QMMM" == "YES" ]
  then
     sed -i "s/DUMMY_QM_CMPLX_MASK/${QM_CMPLX_MASK}/" run_mmpbsa.sh

     declare -a QM_LIST=""
     QM_LIST=($(echo $QM_FRAG_MASK ))
     QM_TERMS=${#QM_LIST[@]}
     if [ $QM_TERMS  != $NFRAG ]
     then
         echo "Inconsistent definition of QM_FRAG_MASK=$QM_FRAG_MASK and NFRAG=$NFRAG"
         exit
     fi
     for (( i=0; i<=NFRAG-1; i++))
     do 
         let "ifrag=$i+1"
         if [ ${ifrag} -lt 10 ]; then txt="0${ifrag}"; else  txt=$ifrag; fi
         j=${QM_LIST["$i"]}
         if [ ${j} != "NONE" ] && [ ${j} != "none" ]
         then 
             sed -i "s/DUMMY_QM_FRAG_MASK_${txt}/${j}/" run_mmpbsa.sh
         else
             sed -i "s/DUMMY_QM_FRAG_MASK_${txt}//" run_mmpbsa.sh
         fi
     done
     unset QM_LIST

     declare -a QM_LIST=""
     QM_LIST=($(echo $QM_FRAG_CHARGE ))
     QM_TERMS=${#QM_LIST[@]}
     if [ $QM_TERMS  != $NFRAG ]
     then
         echo "Inconsistent definition of QM_FRAG_CHARGE=$QM_FRAG_CHARGE and NFRAG=$NFRAG"
         exit
     fi
     QM_CMPLX_CHARGE_tmp=0
     for (( i=0; i<=NFRAG-1; i++))
     do 
         let "ifrag=$i+1"
         if [ ${ifrag} -lt 10 ]; then txt="0${ifrag}"; else  txt=$ifrag; fi
         j=${QM_LIST["$i"]}
         sed -i "s/DUMMY_QM_FRAG_CHARGE_${txt}/${j}/" run_mmpbsa.sh
         let "QM_CMPLX_CHARGE_tmp=$QM_CMPLX_CHARGE_tmp + $j"
     done
     if [ -z QM_CMPLX_CHARGE ]
     then
        QM_CMPLX_CHARGE=$QM_CMPLX_CHARGE_tmp
     elif [ $QM_CMPLX_CHARGE -ne $QM_CMPLX_CHARGE_tmp ]
     then
         echo "Inconsistent definition of QM_CMPLX_CHARGE=$QM_CMPLX_CHARGE and " 
         echo "the sum of the QM_FRAG_CHARGE values=$QM_CMPLX_CHARGE_tmp"
         exit
     fi
  
     sed -i "s/DUMMY_QM_CMPLX_CHARGE/${QM_CMPLX_CHARGE}/" run_mmpbsa.sh
     unset QM_LIST

     sed -i "s/DUMMY_D3H4/${D3H4}/" run_mmpbsa.sh

     if [ "$QMMM_ADD_SOLV" == "YES" ]
     then
         declare -a QM_LIST=""
         QM_LIST=($(echo $QM_ADD_WAT_FRAG ))
         QM_TERMS=${#QM_LIST[@]}
         if [ $QM_TERMS  != $NFRAG ]
         then
             echo "Inconsistent definition of QM_ADD_WAT_FRAG=$QM_ADD_WAT_FRAG and NFRAG=$NFRAG"
             exit
         fi
         for (( i=0; i<=NFRAG-1; i++))
         do 
             let "ifrag=$i+1"
             if [ ${ifrag} -lt 10 ]; then txt="0${ifrag}"; else  txt=$ifrag; fi
             j=${QM_LIST["$i"]}
             sed -i "s/DUMMY_QM_ADD_WAT_FRAG_${txt}/${j}/" run_mmpbsa.sh
         done
         unset QM_LIST
     else 
         for (( i=0; i<=NFRAG-1; i++))
         do 
             let "ifrag=$i+1"
             if [ ${ifrag} -lt 10 ]; then txt="0${ifrag}"; else  txt=$ifrag; fi
             sed -i "s/DUMMY_QM_ADD_WAT_FRAG_${txt}/NO/" run_mmpbsa.sh
         done
     fi
     if [ ${QM_ORCA} == "YES" ]
     then 
        TMP_LEVEL=${QM_ORCA_LEVEL//\//\\\/}
        sed -i "s/DUMMY_QM_ORCA_LEVEL/${TMP_LEVEL}/" run_mmpbsa.sh
        sed -i "s/DUMMY_NPROCS_QM/${NPROCS_QM}/" run_mmpbsa.sh
        sed -i "s/DUMMY_QM_ORCA/YES/" run_mmpbsa.sh
        if [ ${QM_ORCA_RELAX} == "YES" ]
        then 
           TMP_LEVEL=${QM_ORCA_RELAX_LEVEL//\//\\\/}
           sed -i "s/DUMMY_QM_ORCA_RELAX_LEVEL/${TMP_LEVEL}/" run_mmpbsa.sh
           sed -i "s/DUMMY_QM_ORCA_RELAX/YES/" run_mmpbsa.sh
        fi
     else
        sed -i "s/DUMMY_QM_ORCA/NO/" run_mmpbsa.sh
        sed -i "s/DUMMY_NPROCS_QM/${NPROCS_QM}/" run_mmpbsa.sh
     fi
     if [ ${QM_TC_RELAX} == "YES" ]
     then 
        sed -i "s/DUMMY_QM_TC_RELAX/YES/" run_mmpbsa.sh
        TMP_LEVEL=${QM_TC_LEVEL//\//\\\/}
        sed -i "s/DUMMY_QM_TC_LEVEL/${TMP_LEVEL}/" run_mmpbsa.sh
     else
        sed -i "s/DUMMY_QM_TC_RELAX/NO/" run_mmpbsa.sh
     fi

  fi
  sed -i "s/DUMMY_QM_CMPLX_MASK//" run_mmpbsa.sh
  sed -i "s/DUMMY_QM_CMPLX_CHARGE//" run_mmpbsa.sh
  sed -i "s/DUMMY_QM_FRAG_MASK_..//" run_mmpbsa.sh
  sed -i "s/DUMMY_QM_FRAG_CHARGE_..//" run_mmpbsa.sh
  sed -i "s/DUMMY_QM_ADD_WAT_FRAG_..//" run_mmpbsa.sh
  sed -i "s/DUMMY_QM_ORCA/NO/" run_mmpbsa.sh
  sed -i "s/DUMMY_QM_TC_RELAX/NO/" run_mmpbsa.sh
  sed -i "s/DUMMY_QM_ORCA_RELAX/NO/" run_mmpbsa.sh
  sed -i "s/DUMMY_NPROCS_QM/1/" run_mmpbsa.sh


  if [ "$ADCK" == "YES" ]
  then 
      if [ "$ADCK_ONLY" == "YES" ]  
      then
           sed -i 's/DUMMY_ADCK_ONLY/YES/' run_mmpbsa.sh
      else
           sed -i 's/DUMMY_ADCK_ONLY/NO/' run_mmpbsa.sh
      fi
      sed -i 's/DUMMY_ADCK/YES/' run_mmpbsa.sh

      declare -a PDBQT_FRAG=""
      source $PDBQT_SRCFILE
      NPDBQT=${#PDBQT_FRAG[@]} 
      let "NPDBQT=$NPDBQT-1"
      if [ "$NPDBQT" -ne "$NFRAG" ]
      then
          echo 'The PDBQT_FRAG array does not contain enough elements'
          echo "NPDBQT=${NPDBQT}"
          echo 'PDBQT_FRAG='
          echo ${PDBQT_FRAG[*]}
          exit
      fi
      for ((ifrag=1;ifrag<=NFRAG;ifrag++))
      do
         if [ ${ifrag} -lt 10 ]; then txt="0${ifrag}"; else  txt=$ifrag; fi
         PDBQT_FILE=${PDBQT_FRAG["$ifrag"]}
         TMP_PDBQT_FILE=${PDBQT_FILE//\//\\\/}
         sed -i "s/DUMMY_PDBQT_FRAG_${txt}/${TMP_PDBQT_FILE}/" run_mmpbsa.sh
      done
      if [ "$NSODIUM" -gt 0 ] && [ -z $PDBQT_WAT ] 
      then
            echo " A PDBQT template for WAT is also required "
            echo " Define PDBQT_WAT variable"      
            exit
      fi
      if [ "$NSODIUM" -gt 0 ] && [ -z $PDBQT_NA ] 
      then
            echo " A PDBQT template for Na+ is also required "
            echo " Define PDBQT_NA variable"      
            exit
      fi
      if [ $NSODIUM -gt 0 ]
      then
           TMP_PDBQT_FILE=${PDBQT_WAT//\//\\\/}
           sed -i "s/DUMMY_PDBQT_WAT/${TMP_PDBQT_FILE}/" run_mmpbsa.sh
           TMP_PDBQT_FILE=${PDBQT_NA//\//\\\/}
           sed -i "s/DUMMY_PDBQT_NA/${TMP_PDBQT_FILE}/" run_mmpbsa.sh
      fi
  else
      sed -i 's/DUMMY_ADCK_ONLY/NO/' run_mmpbsa.sh
      sed -i 's/DUMMY_ADCK/NO/' run_mmpbsa.sh
  fi
  sed -i 's/DUMMY_PDBQT_WAT//' run_mmpbsa.sh
  sed -i 's/DUMMY_PDBQT_NA//' run_mmpbsa.sh
  sed -i "s/DUMMY_PDBQT_FRAG_..//" run_mmpbsa.sh
  
  INPUT_LIST="LISTA"
  NINPUT=$(cat $INPUT_LIST | wc -l)
  let " NSPLIT =  ( $NINPUT / $NPROCS ) + 1 "

  for PDIE in $PDIE_LIST 
  do

     TASK="./run_mmpbsa_${PDIE}.sh"
     sed  "s/DUMMY_PDIE/${PDIE}/"  run_mmpbsa.sh > $TASK 
     chmod 755 $TASK
    
     rm -f TASK.sh
     rm -f temp_task_list_*
    
     split -l $NSPLIT -d $INPUT_LIST temp_task_list_
    
     for file in $(ls temp_task_list_*)
     do
       echo " $TASK  $file >  $SCRATCH/${file}.log " >> TASK.sh
     done
    
     echo "Running parallel $TASK across  $NPROCS  procs ..."
     if [ $QMMM == "YES" ]
     then
         let "NPROCS=$NPROCS/$NPROCS_QM"
     fi
     cat TASK.sh  | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS
    
     rm -f temp_task_list_* $SCRATCH/temp_task_list_*log $TASK
    
     if [ "${USEOUTPUT}" == "YES" ]  && [ -e PDIE_${PDIE}/OUTPUT.tar ]
     then 
        tar xvf PDIE_${PDIE}/OUTPUT.tar  ${USEOUTPUT_MASK} 
     elif [ "${USEOUTPUT}" == "YES" ]  && [ ! -e PDIE_${PDIE}/OUTPUT.tar ] 
     then 
        echo "USEOUTPUT=YES, but PDIE_${PDIE}/OUTPUT.tar does not exist!"
     fi

     export PDIE="$PDIE"
     export PERCEN="0"
     export OUTLYER="0"
     export DO_ENERGY_ANALYSES_TAR="1"
     export DO_STAT_PLOT="1"
     export DO_STAT_PERCEN=0
     if [ ! -e PDIE_${PDIE} ]; then mkdir PDIE_${PDIE} ; fi 
     if [ "${INCRLIST}" == "YES" ]
     then
        echo "Merging new calcs with previous output files (if available)" 
        if [ "$ADCK_ONLY" == "NO" ] && [ -e  PDIE_${PDIE}/OUTPUT.tar ]; then cp PDIE_${PDIE}/OUTPUT.tar . ; tar xvf OUTPUT.tar; rm -f OUTPUT.tar  ; fi
        if [ "$ADCK" == "YES" ] && [ -e  PDIE_${PDIE}/OUTPUT_DLG.tar ]; then cp PDIE_${PDIE}/OUTPUT_DLG.tar . ; tar xvf OUTPUT_DLG.tar; rm -f OUTPUT_DLG.tar  ; fi
        if [ "$ADCK_ONLY" == "YES" ]
        then
            ls *.dlg1 | sed 's/dlg1/pdb/'  > LISTA
        else 
            ls *.out3.gz | sed 's/out3.gz/pdb/'  > LISTA
        fi
     fi
     if [ "$ADCK_ONLY" == "NO" ]; then $APTAMD/MMPBSA/energy_data_parser.sh LISTA ; fi
     if [ "$ADCK" == "YES" ]
     then 
          for ((isite=1;isite<=NFRAG-1;isite++))
          do
              echo "# PREFIX E_INT_ADCK_SITE_${isite}" > E_INT_ADCK_SITE${isite}.dat
              grep 'Total   ' *.dlg${isite} | grep -v Elec | sed "s/.dlg${isite}://" | awk '{printf(" %s %10.3f \n",$1,$3)}' >> E_INT_ADCK_SITE${isite}.dat
              tar uvf OUTPUT_DLG.tar *.dlg${isite}
              rm -f *.dlg${isite}
          done
          mv OUTPUT_DLG.tar PDIE_${PDIE}/
     fi
     $APTAMD/MMPBSA/stat_plot.sh 
     mv *.dat*  *.med *.png OUTPUT.tar PDIE_${PDIE}/
     

  done

  cd $WORKDIR_TRJ 

  unset TOPOLOGY 


done
