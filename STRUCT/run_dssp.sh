#/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
source $APTAMD/ENV/aptamd_env.sh

if [ ! -e $DSSPHOME/bin/mkdssp ]; then echo "mkdsssp binary needed, but not available"; exit; fi

# Topology file 
TOPOLOGY=$1
# Number of residues
NRES=$2
#  Number of processes
NPROCS=$3
#  Option for dealing with coordinates
OPTION_CRD=$4
#  Option for dealing with coordinates
PDBLIST=$5

# Most of the data analysis and plotting is performed using octave scripts
GNUPLOT=$(which gnuplot | grep -v alias)

if [ -z "$TOPOLOGY" ] ||  [ -z "$NRES" ] ||  [ -z "$NPROCS" ] ||  [ -z "$OPTION_CRD" ]
then
   echo " "
   echo "Usage: run_dssp.sh  topology_file NRES NPROCS OPTION_CRD [PDBLIST] "
   echo "   OPTION_CRD= 0  --> Use all md_*_solute.crd files" 
   echo "   OPTION_CRD= N  --> Add data from last N  md_*_solute.crd files" 
   echo "   OPTION_CRD= -1 --> Redo statistical analysis for files in OUTPUT.tar"
   exit
fi

if [ -n "$OPTION_CRD" ] && [ "$OPTION_CRD" -eq "$OPTION_CRD" ] 2>/dev/null; then
  if [ ${OPTION_CRD} -eq 0 ]
  then 
      INCR_LIST="NO"
      REDO_ANALYSIS="NO"
      FRAME="frame"
  elif [ $OPTION_CRD -eq "-1" ]
  then
      INCR_LIST="NO"
      REDO_ANALYSIS="YES"
      echo "DSSP files in OUTPUT.tar will be reanalyzed"
  elif [ $OPTION_CRD -gt 0  ]
  then
      INCR_LIST="YES"
      REDO_ANALYSIS="NO"
      echo "Only the last $OPTION_CRD MDCRD will be processed"
      FRAME="frame"
  fi
else
  echo "OPTION_CRD = $OPTION_CRD is not a number" 
fi

if [ $INCR_LIST == "YES" ] && [ ! -e OUTPUT.tar ]
then
   echo "INCR_LIST=YES, but OUTPUT.tar does not exist!"
   INCR_LIST="NO"
fi
if [ $REDO_ANALYSIS == "YES" ] && [ ! -e OUTPUT.tar ]
then
   echo "REDO_ANALYSIS=YES, but OUTPUT.tar does not exist!"
   exit
fi

if [ ! -z "${PDBLIST}" ]
then
   nfile=0
   ngzip=0
   for file in $(cat $PDBLIST)
   do
      let "nfile=nfile+1"
      if [ ! -e $file ] && [ ! -e ${file}.gz ]
      then
         echo "PDBLIST supplied, but $file does not exist!"
         exit
      fi
      if [ -e ${file}.gz ]; then let "ngzip=ngzip+1" ; fi
   done
   DO_PDB="YES"
   if [ $nfile -eq $ngzip ]; then PDBGZIP="YES"; else PDBGZIP="NO"; fi 
else
   DO_PDB="NO"
fi

if [ -z "${SIEVE}" ]
then
  SIEVE=500
else
  echo "Using SIEVE=$SIEVE"
fi
if [ -z "${PREFIX_MDCRD}" ]
then
   PREFIX_MDCRD="md_"
else
   echo "Using  PREFIX_MDCRD=${PREFIX_MDCRD}"
fi
if [ -z "${SUFFIX_MDCRD}" ]
then
   SUFFIX_MDCRD="_solute.mdcrd"
else
   echo "Using  SUFFIX_MDCRD=${SUFFIX_MDCRD}"
fi
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

# Temp directory
SCRATCH=/dev/shm 
TT=$(date +%N)
TMPDIR=${SCRATCH}/TMPDIR_${TT}

WORKDIR=$PWD
mkdir $TMPDIR

if [ ${REDO_ANALYSIS} == "NO" ]
then 

cat <<EOF > top_pdb
HEADER    DUMMY                                   01-JAN-01   1DUM              
TITLE     DUMMY
COMPND    MOL_ID: 1;
SOURCE    MOL_ID: 1;
KEYWDS    DUMMY
EXPDTA    THEORETICAL MODEL
AUTHOR    DUMMY
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          
EOF

rm -f pdb.in 
if [ ${DO_PDB} == "YES" ]
then
   for file in $(cat $PDBLIST)
   do
       if [ $PDBGZIP == "YES" ] 
       then 
           zcat ${file}.gz | grep -v 'WAT\|Na+\|Cl-' > $TMPDIR/$file 
       else
           cat ${file} | grep -v 'WAT\|Na+\|Cl-' > $TMPDIR/$file 
       fi
       echo "trajin  $TMPDIR/$file  " >> pdb.in 
   done
else 

   if [ -e ../../5.PRODUCTION/ ]
   then
      PROD_DIR="../../5.PRODUCTION/"
   elif [ -e ../../../5.PRODUCTION/ ]
   then
      PROD_DIR="../../../5.PRODUCTION/"
   else
      echo "5.PRODUCTION directory not found !"
      exit
   fi
   for file in $(ls ${PROD_DIR}/${PREFIX_MDCRD}?${SUFFIX_MDCRD} ) 
   do
      echo "trajin  $file  1 last $SIEVE" >> pdb.in 
   done
   for file in $(ls ${PROD_DIR}/${PREFIX_MDCRD}??${SUFFIX_MDCRD} )
   do
      echo "trajin  $file  1 last $SIEVE" >> pdb.in 
   done
   for file in $(ls ${PROD_DIR}/${PREFIX_MDCRD}???${SUFFIX_MDCRD} )
   do
      echo "trajin  $file  1 last $SIEVE" >> pdb.in 
   done
fi

if [ $INCR_LIST == "YES" ]
then
   tail -$OPTION_CRD pdb.in > tmp; mv tmp pdb.in
fi

if [ $PERCEN -ne 0 ]
then
  NFILES=$(cat pdb.in | wc -l)
  echo "Using only $PERCEN % of the available data"
  echo "  > 0  --> From the begining"
  echo "  < 0  --> From the end"
  if [ "$PERCEN" -gt 0 ]
  then
    let  "NFILES_use=( $NFILES *  $PERCEN ) / 100 "
    head -${NFILES_use}  pdb.in > tmp; mv tmp pdb.in
    echo "Using only the first $NFILES_use trajectory files"
  else
    let  "NFILES_use=($NFILES * ( - $PERCEN ) ) / 100  "
    tail -${NFILES_use}  pdb.in > tmp; mv tmp pdb.in
    echo "Using only the last $NFILES_use trajectory files"
  fi
fi

echo "trajout  ${TMPDIR}/${FRAME}.pdb  pdb  multi pdbv3 chaind A keepext  pdb " >> pdb.in 
echo "go" >> pdb.in

mpirun -np $NPROCS cpptraj.MPI  $TOPOLOGY < pdb.in > pdb.log 

rm -f TASK.sh

for file in $(ls ${TMPDIR}/${FRAME}*pdb  | sort -n -t '.' -k 2 )
do
  echo " cat $WORKDIR/top_pdb $file > ${file}_ext; $APTAMD/AUXTOOLS/mkdssp --output-format dssp  ${file}_ext   ${file/pdb/out}  " >> TASK.sh
done
echo "Running parallel $TASK across  $NPROCS  procs ..."
cat TASK.sh | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS  
rm -f $TMPDIR/${FRAME}*pdb

fi   # ENDIF of REDO_ANALYSIS == "NO"

if [ $INCR_LIST == "YES" ] || [ $REDO_ANALYSIS == "YES" ]
then
    echo "Unpackaging previous OUTPUT.tar file"
    cd $TMPDIR
    NUMFRA_PREV=$(tar tvf $WORKDIR/OUTPUT.tar | wc -l) 
    NSLASH=$(tar tvf $WORKDIR/OUTPUT.tar  | awk '{print $NF}' | head -1 | awk -F '/' '{print NF-1 }')
    for file in $(ls ${FRAME}*out  | sort -n -t '.' -k 2 )
    do
       idfile=$(echo $file |   sed 's/[^0-9]//g' )
       let  idfile="$idfile + $NUMFRA_PREV"
       mv $file ${FRAME}.${idfile}.out
    done
    tar --strip-components=${NSLASH} -xvf $WORKDIR/OUTPUT.tar
    cd $WORKDIR
fi

echo "$NRES" > LISTA_OUT
ls ${TMPDIR}/${FRAME}*.out |  sort -n -t '.' -k 2   >> LISTA_OUT
echo "1 $NRES" > info_for_gnuplot.dat 
$APTAMD/AUXTOOLS//stat_dssp < LISTA_OUT > summary.dat

if [ $REDO_ANALYSIS == "NO" ]
then
  rm -f OUTPUT.tar
  tar cvf OUTPUT.tar ${TMPDIR}/${FRAME}*.out 
fi

# Preparing GNUPLOT script 
YTIC=""
let "nlines=$NRES+1" 
declare -a AACODE=""
AACODE=($(grep -A${nlines} 'AVERAGES' summary.dat | awk '{print $4}'))
for ((ires=1;ires<=NRES;ires++))
do
    let "j=$ires-1"
    TMP=' " '${AACODE["$j"]}' " '${ires}','
    YTIC="${YTIC} ${TMP}"
done
YTIC=${YTIC:0:-1}

NFRAMES=$(cat LISTA_OUT | wc -l)

if [ $NRES -lt 10 ]
then
   PS="10"
else
  PS="5"
fi
PT="4"

cat <<EOF > input.gp
#
set terminal png nocrop enhanced font "/usr/share/fonts/dejavu/DejaVuSerifCondensed.ttf"  22 size 1280,960

set output 'plot_notic.png'
set border 4095 lt -1 lw 3.000

unset xtics
unset ytics
unset border
unset colorbox
unset key
set palette model RGB defined (0 'grey30', 1 'grey70', 2 'red', 3 'blue', 4 'magenta')
#  grey30 --> COIL (T+C+B)
#  grey70 --> BEND (S)
#  red    --> HELIX (H+G+I)
#  blue   --> beta-STRAND (I)
#  magenta   --> PPII (I)

set xtics 0,500,${NFRAMES}
set ytics 1,1,${NRES} 
unset ytics
set ytics  ( $YTIC )

plot[0:${NFRAMES}][0.5:${NRES}+0.5]  'info_for_gnuplot.dat'  u 1:2:(\$3) w p pt ${PT} ps ${PS} palette, \\
EOF

for ((i=2;i<=NRES-1;i++))
do
    let "j=2*$i"
    let "k=$j+1"
    cat <<EOF >> input.gp
'' u 1:DUMMY_I:(\$DUMMY_J) w p pt ${PT} ps ${PS} palette, \\
EOF
    sed -i "s/DUMMY_I/${j}/" input.gp
    sed -i "s/DUMMY_J/${k}/" input.gp
done

let "j=2*$NRES"
let "k=$j+1"
cat <<EOF >> input.gp
'' u 1:DUMMY_I:(\$DUMMY_J) w p pt ${PT} ps ${PS} palette
EOF
    sed -i "s/DUMMY_I/${j}/" input.gp
    sed -i "s/DUMMY_J/${k}/" input.gp

if [ ! -z "${GNUPLOT}" ]
then 
    $GNUPLOT < input.gp
fi
rm -f -r ${TMPDIR}

