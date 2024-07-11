#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
if [ "$#" -eq 0 ]; then more $APTAMD/DOC/do_reweight_gamd.txt ; exit; fi

source $APTAMD/ENV/geuo_env.sh
PYTHON=/opt/apps/SL7/anaconda3/bin/python3
if [ -z "$PYTHON" ]; then echo 'PYTHON is not available, but needed'; exit; fi

# Checking options
extension="${1##*.}"
if [ $extension == "src" ]
then
    echo "Sourcing $1"
    source $1
else
# Put into double quotes the list of trajectories to be processed
 MD_TRAJ=$1
fi

# GAMD TYPE ONLY
MD_TYPE="GAMD"

if [ ! -z "$MOL" ] && [ -z "$MD_TRAJ" ]; then MD_TRAJ="${MOL}"; echo $MD_TRAJ; unset MOL; fi


if [ -z "$MD_TRAJ" ]; then more $APTAMD/DOC/do_reweight_gamd.txt ; exit ; fi

# Options for energy reweighting. These may have to be changed. For example, the EMAX
# parameter should be adapted to the depth of the free energy wells.

if [ -z "$EMAX" ]; then EMAX="8.0"; fi
if [ -z "$FINEGRID" ]; then FINEGRID="NO"; fi
if [ -z "$NGRID_FINE" ]; then NGRID_FINE="25"; fi
if [ -z "$NGRID_REG" ]; then NGRID_REG="15"; fi
if [ -z "$MD_PROD" ]; then MD_PROD="5.PRODUCTION"; fi
if [ -z "$RW_DIR" ]; then RW_DIR="RW_GAMD_ENE"; fi
if [ -z "$DO_INF_RMSD" ]; then DO_INF_RMSD="YES"; fi
if [ -z "$DO_RGYR_RMSD" ]; then DO_RGYR_RMSD="NO"; fi
if [ -z "$DATA_DIR" ]; then DATA_DIR="STRUCT"; fi
if [ -z "$PREFIX_MDCRD" ]; then echo 'Considering PREFIX_MDCRD=md_';  PREFIX_MDCRD="md_"; else echo "Using PREFIX_MDCRD=$PREFIX_MDCRD"; fi
if [ -z "$SUFFIX_MDCRD" ]; then echo 'Considering SUFFIX_MDCRD=_solute.mdcrd';  SUFFIX_MDCRD="_solute.mdcrd"; else echo "Using SUFFIX_MDCRD=$SUFFIX_MDCRD"; fi
if [ -z "$SIEVE" ]; then echo 'Considering SIEVE=1';  SIEVE=1; else echo "Using SIEVE=$SIEVE"; fi
if [ -z "$SUFFIX_GAMD" ]; then echo 'Considering SUFFIX_GAMD=.gamd_log';  SUFFIX_GAMD=".gamd_log"; else echo "Using SUFFIX_GAMD=$SUFFIX_GAMD"; fi


echo "GAMD ENERGY REWEIGHTING JOB"
echo "==========================="
echo "Considering EMAX=$EMAX"
echo "            FINEGRID=$FINEGRID"
echo "            NGRID_FINE=$NGRID_FINE"
echo "            NGRID_REG=$NGRID_REG"

WORKDIR_TRJ=$PWD

for MOL in $MD_TRAJ 
do

  echo $MOL
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
        echo ".... maybe this is not going not work!"
        TOPOLOGY=${MOL}_solute.top
     fi
  fi

  INDIR="${MOL}_${MD_TYPE}/6.ANALYSIS/${RW_DIR}"
  if [[ $WORKDIR_TRJ == *"${INDIR}"* ]]
  then

     echo "We are already located in $INDIR"
     if [ ! -e ../../$TOPOLOGY ]
     then 
      echo "$TOPOLOGY not found in ${MOL}_MD. Exiting"
      exit
     fi

  else

     if [ -e ${MOL}_${MD_TYPE} ]
     then
       cd ${MOL}_${MD_TYPE} 
     else
       echo "${MOL}_${MD_TYPE} does not exist. Exiting!"
       exit
     fi
     if [ ! -e $TOPOLOGY ]
     then 
         echo "$TOPOLOGY not found in ${MOL}_MD. Exiting"
         exit
     fi
     if [ ! -e 6.ANALYSIS ]
     then 
        echo "6.ANALYSIS directory not found in ${MOL}_MD. Exiting"
        exit
     fi
     cd 6.ANALYSIS
     if [ ! -e ${DATA_DIR}  ]
     then 
        echo "${DATA_DIR} directory not found in ${MOL}_MD/6.ANALYSIS. Exiting"
        exit
     fi
     if [ ! -e ${RW_DIR} ]
     then 
         echo "Creating ${RW_DIR} directory in ${MOL}_MD/6.ANALYSIS"
         mkdir ${RW_DIR} 
     else
         echo "Found ${RW_DIR} directory in ${MOL}_MD/6.ANALYSIS"
     fi
     cd ${RW_DIR} 

  fi

  # Printing out options
  cat <<EOF > options.txt
MOL=$MOL     
TOPOLOGY=$TOPOLOGY

RW_DIR=$RW_DIR
DATA_DIR=$DATA_DIR

EMAX=$EMAX
FINEGRID=$FINEGRID
NGRID_FINE=$NGRID_FINE
NGRID_REG=$NGRID_REG
SIEVE=$SIEVE
DO_INF_RMSD=$DO_INF_RMSD
DO_RGYR_RMSD=$DO_RGYR_RMSD

MD_PROD=$MD_PROD
SUFFIX_MDCRD=$SUFFIX_MDCRD
PREFIX_MDCRD=$PREFIX_MDCRD

EOF

  WORKDIR=$PWD 

# Pasting data

  if [ -e gamd.log ] && [ -e inf_rmsd.dat ] && [ -e ID.dat ] && [ -e rgyr_rmsd.dat ]
  then 

      echo "Found gamd.log inf_rmsd.dat ID.dat rgyr_rmsd.dat "
      echo "These will be reutilized for the 2D map"
      echo "If you wish to rebuild them, just erase them before"
      echo "running this script"
  
  else
   
  rm -f gamd.log inf_rmsd.dat ID.dat rgyr_rmsd.dat

  file=$(ls ../../${MD_PROD}/md_001.gamd_log | head -1)
  if [ -z $file ]; then echo " ../../${MD_PROD}/md_001.gamd_log missing" ; exit ;fi
  head -3 $file > gamd.log
 
  ls ../../${MD_PROD}/${PREFIX_MDCRD}???${SUFFIX_MDCRD} | sort > LIST_MDCRD
  ls ../../${MD_PROD}/${PREFIX_MDCRD}???${SUFFIX_GAMD} | sort > LIST_LOG
  ls   ../${DATA_DIR}/*.inf > LIST_INF
  ls   ../${DATA_DIR}/*.rmsd > LIST_RMSD
  ls   ../${DATA_DIR}/*.rgyr > LIST_RGYR
  nmdcrd=$(cat LIST_MDCRD | wc -l) 
  ninf=$(cat LIST_INF | wc -l) 
  nrmsd=$(cat LIST_RMSD | wc -l) 
  nrgyr=$(cat LIST_RGYR | wc -l) 

  if [ $nmdcrd -ne $ninf ]  
  then
     echo "There are $nmdcrd mdcrd files, but $ninf INF files"
     echo "Cannot do RMSD / INF reweight"
     DO_INF_RMSD="NO"
  fi 
  if [ $nmdcrd -ne $nrmsd ]  
  then
     echo "There are $nmdcrd mdcrd files, but $nrmsd RMSD files"
     echo "Can do nothing"
     exit 									
  fi 

  nsnap=0
  rm -f SNAP.dat
  echo '# icrd, n_(snap_in_INF)  n_(snap_in_RMSD)   nsnap_tot(snap_in_MDCRD)' > SNAP.dat
  for ((i=1;i<=nmdcrd;i++))
  do
    rm -f temp.rmsd temp.rgyr temp.inf  temp.log 

    file=$(sed -n "${i},${i}p" LIST_MDCRD)
    mdcrd=$(basename $file) 
    nsnap_tot=$(ncdump -h $file | grep frame | grep UNLI | awk '{print $6}' | sed 's/(//')
    file=$(sed -n "${i},${i}p" LIST_LOG)
    grep -v '#'  $file  > temp.log 
    if [ ${SIEVE} -gt 1 ]
    then
        $OCTAVE -q <<EOF
A=load('temp.log');
s=${SIEVE};
n=length(A(:,1));
B=A([1:s:n],:);
fid=fopen('temp_sieve.log','w');
for i=[1:fix(n/s)]
   fprintf(fid," %15i %15i %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f \n",B(i,:));
end
fclose(fid);
EOF
        mv -f temp_sieve.log temp.log
    fi

    file=$(sed -n "${i},${i}p" LIST_RMSD)
    id=$(basename $file) 
    id=${id%%.*} 
    grep -v '#' ../${DATA_DIR}/${id}.rmsd    | awk        '{print $2}'  > temp.rmsd
    grep -v '#' ../${DATA_DIR}/${id}.rgyr    | awk        '{print $2}'  > temp.rgyr
    nrmsd=$(cat temp.rmsd| wc -l)
    if [ ${DO_INF_RMSD} == "YES" ]
    then 
      grep -v '#' ../${DATA_DIR}/${id}.inf | awk -F ',' '{print $(NF)*10}'  > temp.inf 
      ninf=$(cat temp.inf  | wc -l)
      if [  "$nrmsd" -eq "$ninf" ]
      then 
         echo "md_${i} has same # of data in INF $ninf and RMS $nrmsd" 
         paste temp.inf  temp.rmsd | sed 's/\t/   /g' >> inf_rmsd.dat 
         paste temp.rgyr temp.rmsd | sed 's/\t/   /g' >> rgyr_rmsd.dat 
         cat temp.log >> gamd.log 
         declare -a ID=""
         for ((j=1;j<=nrmsd;j++))
         do 
           let "k=$nsnap + $j "
           ID["$j"]=$k
         done
         printf "%i \n" ${ID[*]} >>ID.dat 
         unset ID 
         let "nsnap=$nsnap+$nrmsd"
      elif [ "$ninf" -lt "$nrmsd" ] 
      then
         echo "md_${i} has different # of data in INF $ninf and RMS $nrmsd"
         echo "Getting only $nrmsd data from RMSD and GAMD log"
         head -$ninf temp.rmsd > temp ; mv temp temp.rmsd
         paste temp.inf temp.rmsd | sed 's/\t/   /g' >> inf_rmsd.dat
         head -$ninf temp.rgyr > temp ; mv temp temp.rgyr
         paste temp.rgyr temp.rmsd | sed 's/\t/   /g' >> rgyr_rmsd.dat 
         head -$ninf temp.log> temp ; mv temp temp.log
         cat temp.log >> gamd.log 
         declare -a ID=""
         for ((j=1;j<=ninf;j++))
         do 
           let "k=$nsnap+$j"
           ID["$j"]=$k
         done
         printf "%i \n" ${ID[*]} >>ID.dat 
         unset ID 
         let "nsnap=$nsnap+$nrmsd"
      else 
         echo "md_${i} has different # of data in INF $ninf > RMS $nrmsd" 
         echo "Not used!"
      fi
      echo "$i  $ninf   $nrmsd   $nsnap_tot " >> SNAP.dat 
    else
         paste temp.rgyr temp.rmsd | sed 's/\t/   /g' >> rgyr_rmsd.dat 
         cat temp.log >> gamd.log 
         declare -a ID=""
         for ((j=1;j<=nrmsd;j++))
         do 
           let "k=$nsnap+$j"
           ID["$j"]=$k
         done
         printf "%i \n" ${ID[*]} >>ID.dat 
         unset ID 
         let "nsnap=$nsnap+$nrmsd"
    fi 
  done

  rm -f temp.rmsd temp.log temp.rgyr temp.inf

  fi

# Using python scripts for energy reweighting 


  nlines=$(cat ID.dat | wc -l)   # number of data points used for reweighting
 
####################################################
# Prepare input file "weights.dat" in the following format: 
# Column 1: dV in units of kbT; column 2: timestep; column 3: dV in units of kcal/mol

  if [ ! -e weights.dat ]
  then 
    tail -n $nlines gamd.log | awk 'NR%1==0' | awk '{print ($8+$7)/(0.001987*300)"                " $2  "             " ($8+$7)}' > weights.dat
  else
  echo "Found weights.dat file. Reutilizing it. Delete it if you wish to rebuild it"
  
  fi

  LIST_COORD=""
  if [ ${DO_INF_RMSD} == "YES" ]  
  then 
     LIST_COORD="inf_rmsd" 
     # Pruning out failed x3dna analyses (-1) 
     ilines=$(grep -n '\-1' inf_rmsd.dat |  sed 's/:/   /' | awk '{printf("%sd;", $1)}' )
     if [ ${ilines} ]
     then
         sed -i "${ilines}"  inf_rmsd.dat
         sed -i "${ilines}"  rgyr_rmsd.dat
         sed -i "${ilines}"  weights.dat
         sed -i "${ilines}"  SNAP.dat
         sed -i "${ilines}"  ID.dat
         head -3 gamd.log > temp.heading
         sed -i '1,3d' gamd.log
         sed -i "${ilines}" gamd.log
         cat temp.heading gamd.log > temp; mv -f  temp gamd.log 
         rm -f temp.heading
     fi
  fi

  if [ ${DO_RGYR_RMSD} == "YES" ] ; then LIST_COORD="${LIST_COORD} rgyr_rmsd" ;fi

  for COORD in $(echo "${LIST_COORD}")   
  do

  if [ ${FINEGRID} == "YES" ]
  then
     COORD_DIR=${COORD}_${NGRID_FINE}
  else
     COORD_DIR=${COORD}_${NGRID_REG}
  fi


  if [ ! -e ${COORD_DIR} ]
  then
      mkdir ${COORD_DIR}
      cd ${COORD_DIR}
          ln -s ../${COORD}.dat ${COORD}.dat
          ln -s ../weights.dat  weights.dat 
          ln -s ../SNAP.dat     SNAP.dat 
          ln -s ../ID.dat       ID.dat 
          ln -s ../gamd.log     gamd.log 
  else
       cd ${COORD_DIR}
  fi

  $OCTAVE -q <<EOF > ${COORD}.max_min
A=load("${COORD}.dat");
printf(' Max  %f  %f  \n',max(A))
printf(' Min  %f  %f  \n',min(A))
printf(' FineGrid  %f  %f  \n',abs(max(A)-min(A))/${NGRID_FINE})
printf(' RegGrid  %f  %f  \n',abs(max(A)-min(A))/${NGRID_REG})
quit
EOF

  Xmax=$(grep Max ${COORD}.max_min | awk '{print $2}') 
  Ymax=$(grep Max ${COORD}.max_min | awk '{print $3}') 
  Xmin=$(grep Min ${COORD}.max_min | awk '{print $2}') 
  Ymin=$(grep Min ${COORD}.max_min | awk '{print $3}') 

  if [ "$FINEGRID" == "YES" ]
  then 
   Xinc=$(grep FineGrid ${COORD}.max_min | awk '{print $2}') 
   Yinc=$(grep FineGrid ${COORD}.max_min | awk '{print $3}') 
  else
   Xinc=$(grep RegGrid ${COORD}.max_min | awk '{print $2}') 
   Yinc=$(grep RegGrid ${COORD}.max_min | awk '{print $3}') 
  fi

# Reweighting using cumulant expansion 
  PY2=$($PYTHON -V | grep -c ' 2.')
  PY3=$($PYTHON -V | grep -c ' 3.')
  echo $PY2  $PY3
  PY3=1
  if [  $PY2 -ge 1 ]
  then
     sed "s/DUMMY_EMAX/${EMAX}/" $APTAMD/RWGAMD/PyReweighting-2D_python2.py > PyReweighting-2D.py
  elif [ $PY3 -ge 1 ]
  then
     sed "s/DUMMY_EMAX/${EMAX}/" $APTAMD/RWGAMD/PyReweighting-2D.py > PyReweighting-2D.py
  else
     echo "$PYTHON -V is not recognized!"
     $PYTHON -V
     exit
  fi
  echo "$PYTHON PyReweighting-2D.py -cutoff 10 -input ${COORD}.dat -Xdim $Xmin $Xmax  -discX $Xinc -Ydim $Ymin $Ymax -discY $Yinc -Emax $EMAX -job amdweight_CE -weight weights.dat | tee -a reweight_variable.log "  > PyReweighting-2D.log 
  $PYTHON PyReweighting-2D.py -cutoff 10 -input ${COORD}.dat -Xdim $Xmin $Xmax  -discX $Xinc -Ydim $Ymin $Ymax -discY $Yinc -Emax $EMAX -job amdweight_CE -weight weights.dat | tee -a reweight_variable.log  >> PyReweighting-2D.log 

  mv -v pmf-c1-${COORD}.dat.xvg pmf-2D-${COORD}-reweight-CE1.xvg
  mv -v pmf-c2-${COORD}.dat.xvg pmf-2D-${COORD}-reweight-CE2.xvg
  mv -v pmf-c3-${COORD}.dat.xvg pmf-2D-${COORD}-reweight-CE3.xvg
  if [ -e 2D_Free_energy_surface.png ]; then mv -v 2D_Free_energy_surface.png pmf-2D${COORD}-reweight-CE2.png; fi 

# Cumulant expansion 2nd order is our choice.... We treat the data file accordingly 
  head -1  pmf-2D-${COORD}-reweight-CE2.xvg > CE2.dat
  sed '1,5d'  pmf-2D-${COORD}-reweight-CE2.xvg >> CE2.dat

#  Using octave to draw the binned 2D MAPs 
  sed "s/DUMMY_EMAX/${EMAX}/" $APTAMD/RWGAMD/energy_reweight.m |sed "s/DUMMY_COORD/${COORD}/" > energy_reweight.m
  if  [ $COORD == "inf_rmsd" ]
  then 
     sed -i 's/DUMMY_RC_1/INF/' energy_reweight.m
  else
     sed -i 's/DUMMY_RC_1/Rgy/' energy_reweight.m
  fi
  sed -i 's/DUMMY_RC_2/RMSD/' energy_reweight.m
  if [ "$FINEGRID" == "YES" ]
  then 
     sed -i 's/DUMMY_GRID/1/' energy_reweight.m
  else
     sed -i 's/DUMMY_GRID/0/' energy_reweight.m
  fi

  $OCTAVE -q energy_reweight.m 

  # Selecting representativa snapshots of the lower energy bins 
   declare -a BIN=""
   declare -a ISNAP=""
   declare -a ICRD=""
   
   BIN=($(grep 'CRD_FILE' 2D_RW.dat | awk '{print $2}'))
   ISNAP=($(grep 'CRD_FILE' 2D_RW.dat | awk '{print $4}'))
   ICRD=($(grep 'CRD_FILE' 2D_RW.dat | awk '{print $6}'))

   NPDB="${#ICRD[@]}"

   for ((ipdb=0; ipdb<=NPDB-1; ipdb++))
   do
     
      i=${ICRD["$ipdb"]}
      file=$(sed -n "${i},${i}p" ../LIST_MDCRD)
      echo "trajin  ../${file}  ${ISNAP["$ipdb"]}  ${ISNAP["$ipdb"]}  1 " > pdb.in 
      echo "trajout  snap_${BIN["$ipdb"]}_${ISNAP["$ipdb"]}.pdb  pdb" >> pdb.in
      if [ $ipdb -gt  0 ] 
      then
          echo "reference snap_${BIN[0]}_${ISNAP[0]}.pdb" >> pdb.in 
	  echo "rms reference !@H=" >> pdb.in
       fi
       echo "go" >> pdb.in 
   
      $AMBERHOME/bin/cpptraj ../../../$TOPOLOGY < pdb.in

      echo "REMARK   SNAP= ${ISNAP["$ipdb"]} BIN= ${BIN["$ipdb"]} CRD= md_${ICRD["$ipdb"]}_prot.mdcrd" > tmp.pdb
      cat snap_${BIN["$ipdb"]}_${ISNAP["$ipdb"]}.pdb >> tmp.pdb
      mv tmp.pdb snap_${BIN["$ipdb"]}_${ISNAP["$ipdb"]}.pdb


    done

    cd ../

  done

  cd $WORKDIR_TRJ 

done


   
  
 
