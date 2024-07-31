#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
if [ "$#" -eq 0 ]; then more $APTAMD/DOC/do_cluster.txt ; exit; fi

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

if [ -z "$MD_TRAJ" ]; then echo 'Usage: do_cluster.sh "MOL1   MOL2 .."   [ MD | GAMD ] '; exit ; fi
if [ -z "$MD_TYPE" ]; then echo 'Usage: do_cluster.sh "MOL1   MOL2 .."   [ MD | GAMD ] '; exit ; fi
if [ -z "$MD_PROD" ]; then MD_PROD="5.PRODUCTION"; else echo "Assuming MD_PROD=$MD_PROD"; fi
if [ -z "$SOLUTE_TYPE" ]; then SOLUTE_TYPE="UNKNOWN"; else echo "Guessing SOLUTE_TYPE"; fi
if [ -z "$SIEVE" ]; then echo 'Considering SIEVE=50';  SIEVE=50; else echo "Using SIEVE=$SIEVE"; fi
if [ -z "$CLUST_ALL" ]; then echo 'Considering CLUST_ALL=NO';  CLUST_ALL="NO"; else echo "Using CLUST_ALL=$CLUST_ALL"; fi
if [ -z "$EPS" ]; then echo 'Considering EPS=2';  EPS=2; else echo "Using EPS=$EPS"; fi
if [ -z "$SUFFIX_MDCRD" ]; then echo 'Considering SUFFIX_MDCRD=_solute.mdcrd';  SUFFIX_MDCRD="_solute.mdcrd"; else echo "Using SUFFIX_MDCRD=$SUFFIX_MDCRD"; fi
if [ -z "$PREFIX_MDCRD" ]; then echo 'Considering PREFIX_MDCRD=md_';  PREFIX_MDCRD="md_"; else echo "Using PREFIX_MDCRD=$PREFIX_MDCRD"; fi
if [ -z "$CLUSTER_DIR" ]; then echo 'Considering CLUSTER_DIR=CLUSTER';  CLUSTER_DIR="CLUSTER"; else echo "Using CLUSTER_DIR=$CLUSTER_DIR"; fi
if [ -z "$MASK" ]; then MASK='NONE'; MASK_READ="NO"; echo "Guessing MASK for RMSD";  else echo "Using MASK=$MASK"; MASK_READ="YES";   fi
if [ -z "$TSNE" ]; then TSNE='NO';  echo "TSNE clustering is not performed";  fi
if [ -z "$DSSP" ]; then DSSP='NO';  echo "DSSP analsysis will not be performed";  fi
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


NEPS=$(echo $EPS | wc -w)
if [ $NEPS -gt 1 ]
then
   echo "$NEPS different values of EPS will be considered!"
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
export OMP_NUM_THREADS=$NPROCS

if [ $DSSP == "YES" ] 
then
  echo "After CPPTRAJ clustering, DSSP analysis will be executed. DSSP output will be placed"
  echo "in a separate DSSP directory"
  echo "This can be a cpu time consuming task so that a careful choice of the SIEVE variable"
  echo "is recommended!" 
fi

if [ $TSNE == "YES" ] 
then
  echo "After CPPTRAJ clustering, TSNE clustering will be executed. TSNE output will be placed"
  echo "in a separate TSNE directory"
  echo "This can be a cpu time consuming task so that a careful choice of the SIEVE variable"
  echo "is recommended!" 
fi
if [ $TSNE == "ONLY" ] 
then
  echo "CPPTRAJ clustering is omitted and only TSNE clustering will be executed. TSNE output will be placed"
  echo "in a separate TSNE directory"
  echo "This can be a cpu time consuming task so that a careful choice of the SIEVE variable"
  echo "is recommended!" 
fi
if [ $TSNE == "ONLY" ]  || [ $TSNE == "YES" ]
then
  if [ -z "$PERP0" ]; then PERP0=100; echo "Lower Perpelixity value= $PERP0"; else "Using Lower Perpelixity value= $PERP0";  fi
  if [ -z "$PERP1" ]; then PERP1=1000;echo "Upper Perpelixity value= $PERP1"; else "Using Upper Perpelixity value= $PERP1";  fi
  if [ -z "$PERPX" ]; then PERPX=100;  echo "Perpelixity increment= $PERPX"; else "Using Perpelixity increment = $PERPX";  fi
  if [ -z "$KCLUS0" ]; then KCLUS0=10; echo "Lower K-clus number value= $KCLUS0"; else "Using Lower K-clus number value= $KCLUS0";  fi
  if [ -z "$KCLUS1" ]; then KCLUS1=60; echo "Upper K-clus number value= $KCLUS1"; else "Using Upper K-clus number value= $KCLUS1";  fi
  if [ -z "$KCLUSX" ]; then KCLUSX=10; echo "K-clus number increment= $KCLUSX"; else "Using K-clus number increment = $KCLUSX";  fi
fi

if [ ${SOLUTE_TYPE} == "APT" ] && [ ${MASK} == "NONE" ]  
then
   MASK="P,OP1,OP2,O5',C5',C4',C3',C2',C1',O4',O3'"
   echo "SOLUTE_TYPE=''APT''   MASK=$MASK"
elif [ ${SOLUTE_TYPE} == "PEP" ] && [ ${MASK} == "NONE" ]
then
   MASK="C,O,CA,N"
   echo "SOLUTE_TYPE=''PEP''   MASK=$MASK"
fi

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

  if [ -z $REFPDB ]
  then 
     REFPDB=${WORKDIR_TRJ}/${MOL}_${MD_TYPE}/1.EDITION/${MOL}_solute.pdb 
  fi 

  if [ ! -e $REFPDB ]
  then 
     echo "$REFPDB not found. Exiting."
     exit
  fi
  NRES=$(grep 'ATOM  '  $REFPDB | tail -1 | awk '{print $5}')
  if [ $NRES -gt 10 ]
  then
     IRES=3
     let "JRES=$NRES-2"
  else 
     IRES=1
     JRES=$NRES
  fi
# Getting information from REFPDB to select the proper mask!
  if [ ${MASK_READ} == "NO" ] &&   [ ${SOLUTE_TYPE} == "PEP" ]
  then
     MOL_MASK=":${IRES}-${JRES}@${MASK}"
  elif [ ${MASK_READ} == "NO" ] &&   [ ${SOLUTE_TYPE} == "APT" ]
  then
     MOL_MASK=":${IRES}-${JRES}@${MASK}"
  elif [ ${MASK_READ} == "YES" ] &&   [ ${SOLUTE_TYPE} == "PEP" ]
  then
     MOL_MASK="${MASK}"
  elif [ ${MASK_READ} == "YES" ] &&   [ ${SOLUTE_TYPE} == "APT" ]
  then
     MOL_MASK="${MASK}"
  elif [  ${SOLUTE_TYPE} == "UNKNOWN" ] 
  then
     if [ $MASK_READ == "YES" ]
     then
        echo "Guessing SOLUTE_TYPE "
     else
        echo "Guessing both SOLUTE_TYPE and MASK "
     fi
     grep 'ATOM ' $REFPDB | awk '{printf(" %s %3i\n",$4,$5)}' | uniq > tmp.res
     nd=0
     for res in DA DT DC DG 
     do 
        n=$(grep -c $res tmp.res)
        let "nd=$nd+$n"
     done
     naa=0 
     for res in ACE ALA ARG ASH ASN ASP CYM CYS CYX GLH GLN GLU GLY \
     HID HIE HIP HYP ILE LEU LYN LYS MET PHE PRO SER THR TRP TYR VAL NME 
     do 
        n=$(grep -c $res tmp.res)
        let "naa=$naa+$n"
     done
     rm -f tmp.res
     if [ $NRES -eq $nd ]
     then 
        echo "Good. This seems to be a nice aptamer"
        SOLUTE_TYPE="APT"
        if [ "$MASK_READ" == "NO" ]
        then
           MASK="P,OP1,OP2,O5',C5',C4',C3',C2',C1',O4',O3'"
           MOL_MASK=":${IRES}-${JRES}@${MASK}"
        else
           MOL_MASK=${MASK} 
        fi
     elif [ $NRES -eq $naa ]
     then
        echo "Good. This seems to be a nice peptide"
        SOLUTE_TYPE="PEP"
        if [ "$MASK_READ" == "NO" ]
        then
          MASK="C,O,CA,N"
          MOL_MASK=":${IRES}-${JRES}@${MASK}"
        else
           MOL_MASK=${MASK} 
        fi
     elif  [ "$MASK_READ" == "YES" ]
     then
        echo "Unknown solute type"
        SOLUTE_TYPE="UNKNOWN"
        MOL_MASK=${MASK}
     else
        echo "Sorry, no luck in guessing solute type and no MASK available"
        echo "Better define your own MASK for clustering."
        exit
     fi
  else
     echo "No MASK for RSMD Clustering, no way"
     exit
  fi

  if [ ! -e 6.ANALYSIS ]
  then 
     echo "6.ANALYSIS directory not found in ${MOL}_${MD_TYPE}. Making it"
     mkdir 6.ANALYSIS 
  fi
  cd 6.ANALYSIS
  if [ ! -e $CLUSTER_DIR ]
  then 
      echo "Creating $CLUSTER_DIR directory in ${MOL}_${MD_TYPE}/6.ANALYSIS"
      mkdir $CLUSTER_DIR  
  else
      echo "Found CLUSTER  directory in ${MOL}_${MD_TYPE}/6.ANALYSIS"
  fi
  cd $CLUSTER_DIR
  WORKDIR=$PWD 

  # Printing out options
  cat <<EOF > options.txt
MOL=$MOL     
TOPOLOGY=$TOPOLOGY
NRES=$NRES
SOLUTE_TYPE=$SOLUTE_TYPE

CLUSTER_DIR=$CLUSTER_DIR
MASK=$MASK
MOL_MASK=$MOL_MASK
SIEVE=$SIEVE
CLUST_ALL=$CLUST_ALL
EPS="${EPS}"
REFPDB=$REFPDB

MD_TYPE=$MD_TYPE
MD_PROD=$MD_PROD
SUFFIX_MDCRD=$SUFFIX_MDCRD
PREFIX_MDCRD=$PREFIX_MDCRD
PERCEN=$PERCEN

TSNE=$TSNE
PERP0=$PERP0
PERP1=$PERP1
PERPX=$PERPX
KCLUS0=$KCLUS0
KCLUS1=$KCLUS1
KCLUSX=$KCLUSX

DSSP=$DSSP

NPROCS=$NPROCS
EOF


# Preparing CPPTRAJ input
  rm -f trajin.txt 
  if [ $TSNE != "ONLY" ]; then rm -f rep*pdb cnumvtime.dat CpptrajPairDist summary.dat info.dat RMSD.dat ; fi 
  ls  ../../${MD_PROD}/${PREFIX_MDCRD}?${SUFFIX_MDCRD} | sort  > tmp.list
  ls  ../../${MD_PROD}/${PREFIX_MDCRD}??${SUFFIX_MDCRD}| sort  >> tmp.list
  ls  ../../${MD_PROD}/${PREFIX_MDCRD}???${SUFFIX_MDCRD} | sort  >> tmp.list
  if [ $PERCEN -ne 0 ]
  then
     NFILES=$(cat tmp.list | wc -l)
     echo "Using only $PERCEN % of the available data"
     echo "  > 0  --> From the begining"
     echo "  < 0  --> From the end"
     if [ "$PERCEN" -gt 0 ]
     then
       let  "NFILES_use=( $NFILES *  $PERCEN ) / 100 "
       head -${NFILES_use}  tmp.list > tmp; mv tmp tmp.list
       echo "Using only the first $NFILES_use trajectory files"
       cat tmp.list
     else
       let  "NFILES_use=($NFILES * ( - $PERCEN ) ) / 100  "
       tail -${NFILES_use}  tmp.list > tmp; mv tmp tmp.list
       echo "Using only the last $NFILES_use trajectory files"
       cat tmp.list
     fi
  fi
  for file in $(cat tmp.list)
  do 
        if [ ${CLUST_ALL} == "YES" ]
        then
            echo "trajin  $file " >> trajin.txt 
        else 
            nsnap=$(ncdump -h  $file  |grep frame | grep UNLI | awk '{print $6}' | sed 's/(//')
            echo "trajin  $file 1 ${nsnap} ${SIEVE} " >>  trajin.txt 
        fi
  done
  rm -f tmp.list
  if [ ${CLUST_ALL} == "YES" ]
  then
    SIEVE_CLUST=$SIEVE
  else 
    SIEVE_CLUST=1
  fi
  echo "autoimage" >> trajin.txt 
  echo "rms first ${MOL_MASK} out RMSD.dat" >>  trajin.txt 
  echo "radgyr  :1-${NRES} out RGYR.dat " >> trajin.txt 

  if [ ${TSNE} == "ONLY" ] && [ ! -e CpptrajPairDist ]
  then

  cp trajin.txt cluster.in
  echo "rms2d ${MOL_MASK} out RMS2D.dat mass"  >> cluster.in
  echo "go " >> cluster.in

  TSNE_INPUT_MATRIX="RMS2D.dat"

  $AMBERHOME/bin/cpptraj.OMP ../../$TOPOLOGY < cluster.in > cluster.log

  elif [ ${TSNE} == "ONLY" ] && [ -e CpptrajPairDist ]
  then

     echo "CpptrajPairDist found and TSNE=ONLY"
     echo "Remove CpptrajPairDist if you do not want to use it!"
     TSNE_INPUT_MATRIX="CpptrajPairDist"

  else

  declare -a EPSARRAY=""
  EPSARRAY=($EPS)
  for ((IEPS=1;IEPS<=NEPS;IEPS++))
  do
  let "J=$IEPS-1"
  EPS_CLUST=${EPSARRAY["$J"]} 

  echo "Running clustering calcs with EPS=$EPS_CLUST"
   
  if [ $NEPS -gt 1 ]
  then
     mkdir EPS_${EPS_CLUST}
  fi

  cp trajin.txt cluster.in 
  cat <<EOF >> cluster.in

# Adapted from  Recipe from Chetham-III paper 
# doi:10.1016/j.bbagen.2014.09.007
# See cpptraj manual !
# Uncomment cluster option lines for further (a lot of)
# output  

cluster hieragglo averagelinkage epsilon ${EPS_CLUST} \\
        rms mass ${MOL_MASK} \\
        sieve ${SIEVE_CLUST} random out cnumvtime.dat \\
        repout rep repfmt pdb \\
        summary summary.dat info info.dat 
#       savepairdist pairdist CpptrajPairDist \\
#       cpopvtime cpop.agr normframe \\
#       singlerepout singlerep.nc singlerepfmt netcdf  \\
#       clusterout cluster clusterfmt netcdf lifetime

go
EOF

  if [ $TSNE == "YES" ] || [ ${NEPS} -gt 1 ]
  then 
    TSNE_INPUT_MATRIX="CpptrajPairDist" 
    sed  -i 's/info.dat/info.dat \\/' cluster.in 
    sed  -i 's/\#       savepairdist pairdist CpptrajPairDist \\/       savepairdist pairdist CpptrajPairDist \n /' cluster.in 
    if [ ${IEPS} -gt 1 ]
    then
       sed -i 's/ savepairdist/ loadpairdist/' cluster.in
    fi
  fi

  $AMBERHOME/bin/cpptraj.OMP ../../$TOPOLOGY < cluster.in > cluster.log

  if [ $NEPS -gt 1 ]
  then
     mv cluster.log rep.*.pdb summary.dat info.dat cnumvtime.dat RMSD.dat RGYR.dat EPS_${EPS_CLUST}
     cp cluster.in EPS_${EPS_CLUST}
  fi

done

fi

if  [ $DSSP == "YES" ]
then
   mkdir DSSP
   cd DSSP 
   if [ $SOLUTE_TYPE != "PEP" ]
   then
      echo "DSSP analysis of peptide/protein Secondary Structure will be executed, "
      echo "but SOLUTE_TYPE=$SOLUTE_TYPE may be not adequate!"
   fi 
   cp $APTAMD/STRUCT/run_dssp.sh .
   export SIEVE=$SIEVE
   export SUFFIX_MDCRD=$SUFFIX_MDCRD
   export PERCEN=$PERCEN
   chmod 755 run_dssp.sh
   ./run_dssp.sh ../../../$TOPOLOGY $NRES $NPROCS 0
   cd ../
fi

if [ $TSNE == "YES" ]  || [ $TSNE == "ONLY" ]
then
   if [ ! -e TSNE ]; then mkdir TSNE ; fi 
   cd TSNE 
   cp $APTAMD/STRUCT/TSNE_KMEANS_CLUST.py  .
   sed -i "s/DUMMY_FILE/${TSNE_INPUT_MATRIX}/" TSNE_KMEANS_CLUST.py
   sed -i "s/DUMMY_PERP0/${PERP0}/" TSNE_KMEANS_CLUST.py
   sed -i "s/DUMMY_PERP1/${PERP1}/" TSNE_KMEANS_CLUST.py
   sed -i "s/DUMMY_PERPX/${PERPX}/" TSNE_KMEANS_CLUST.py
   sed -i "s/DUMMY_KCLUS0/${KCLUS0}/" TSNE_KMEANS_CLUST.py
   sed -i "s/DUMMY_KCLUS1/${KCLUS1}/" TSNE_KMEANS_CLUST.py
   sed -i "s/DUMMY_KCLUSX/${KCLUSX}/" TSNE_KMEANS_CLUST.py
   $PYTHON/python3 TSNE_KMEANS_CLUST.py  > tsne.log 
#  rm -f ../${TSNE_INPUT_MATRIX} 
   cd ../
   grep 'trajin\|autoimage\|rms first' cluster.in > pdb_tsne.in
   REPFRAME=$(grep -v '#' TSNE/summary_optimal_TSNE.txt | grep BEST | awk '{printf("%i,",$1)}')
   REPCLUS=$(grep -v '#' TSNE/summary_optimal_TSNE.txt | grep BEST | awk '{printf("%s,",$2)}')
   REPFRAME=${REPFRAME:0:-1}
   echo "trajout rep_tsne.pdb multi include_ep onlyframes ${REPFRAME}" >> pdb_tsne.in
   $AMBERHOME/bin/cpptraj.OMP ../../${TOPOLOGY}  < pdb_tsne.in > pdb_tsne.log
   mv rep_tsne*  TSNE/
   cp cluster.in TSNE/
   cd TSNE
   REPFRAME_ARRAY=(${REPFRAME//,/ })
   REPCLUS_ARRAY=(${REPCLUS//,/ })
   NFRAME=${#REPFRAME_ARRAY[@]}
   for ((I=1;I<=NFRAME;I++))
   do
       let "J=$I-1"
       IFRAME=${REPFRAME_ARRAY["$J"]}
       ICLUS=${REPCLUS_ARRAY["$J"]}
       mv rep_tsne.pdb.${IFRAME} rep_tsne.${ICLUS}.pdb
   done
   tar cvfz  TSNE_data.tgz kmeans_* tsnep*
   if [ -e ../DSSP/info_for_gnuplot.dat ]
   then
      ln -s ../DSSP/info_for_gnuplot.dat info_for_gnuplot.dat
   fi
   rm -f  kmeans_* tsnep* status.txt
   cd ../

fi 

  cd $WORKDIR_TRJ 

done


   
  
 
