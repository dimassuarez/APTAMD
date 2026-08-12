#!/bin/bash
#
# This script helps estimate the 
# conformational entropy from AMBER topology file and trajectory files
#
if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
source $APTAMD/ENV/aptamd_env.sh
##########################################################################3
# Other programs being used

#CENCALC programs
CENCALC="$CENCALC_PATH/cencalc_ccmla"
PREP="$CENCALC_PATH/cencalc_prep"
GETTOR="$CENCALC_PATH/get_tor.py" 
export LD_LIBRARY_PATH=$DISLIN:$LD_LIBRARY_PATH

#  NCDUMP converts netCDF file to text format. NETCDF library is required (usually bundled with AMBER).
NCDUMP=$(which  ncdump| grep -v alias)
if [ -z "$NCDUMP" ]; then echo 'ncudmp seems to be not available, but needed!';  DO_EXIT=1; fi

# AMBER trajectory analysis software
CPPTRAJ="$AMBERHOME/bin/cpptraj"
CPPTRAJ_OMP="$AMBERHOME/bin/cpptraj.OMP"

#  GNU Parallel builds and runs commands in parallel. 
PARALLEL="$PARHOME/bin/parallel --no-notice "


##########################################################################3
# The following environmental variables specify the solute mask, 
# topology file, trajectory data, details of CC-MLA, ...
echo "calc_sconform.sh is running...."
DO_EXIT=0

# CCMLA options
if [ -z "$DO_CC_MLA" ]
then
   DO_CC_MLA="0"
fi
if [ -z "$DO_S2" ]
then
   DO_S2="0"
fi
if [ $DO_CC_MLA -eq 1 ]; then echo "Processing trajectory and doing CC-MLA calcs "; fi
if [ $DO_CC_MLA -eq 2 ]; then echo "Doing CC-MLA calcs (previously processsed trajectory) "; fi

if [ $DO_S2 -eq 1 ]; then echo "Processing trajectory and doing S2 calcs "; fi
if [ $DO_S2 -eq 2 ]; then echo "Doing S2 calcs (previously processsed trajectory) "; fi

# Composite CC-MLA ?
if [ -z "$DO_COMP_CC_MLA" ]
then
   DO_COMP_CC_MLA="0"
fi
if [ $DO_CC_MLA -gt 0 ] && [ $DO_COMP_CC_MLA -eq  1 ]; then echo "Trying composite CC-MLA calcs "; fi

# Cutoff values for CCMLA   (CUTOFF=-1 ---> No cutoff)
if [ -z "$CUTOFF" ]
then
   export CUTOFF="6 8 10"
   echo "Using CUTOFF=$CUTOFF values for CCMLA calcs (if requested)"
else
   echo "Using CUTOFF=$CUTOFF values for CCMLA calcs as predefined"
fi

if [ $DO_CC_MLA -le 1 ] && [ $DO_S2 -le 1 ] 
then

# MASK DEFINITION  
if [ -z "$SLVNTMASK" ]
then
   SLVNTMASK=":WAT,Na+,Cl-,MG"
fi
echo "SLVNTMASK=${SLVNTMASK}  including counterions"

if [ -z "$MOL_MASK" ]
then
   echo "MOL_MASK not defined... then keeping only solute atoms"
   MOL_MASK="!${SLVNTMASK}"
fi
echo "Using MOL_MASK=$MOL_MASK"
# GETTOR OPTIONS 
if [ -z "$GETTOR_OPT" ]
then
   echo "GETTOR_OPT not defined... then selecting -puck -noMet "
   GETTOR_OPT="-puck -noMet"
fi
echo "Using GETTOR_OPT=$GETTOR_OPT"

# TOPOLOGY FILE (Full path name)
if [ -z "$REFTOP" ]
then
   echo "Topology file not defined, but needed." 
   DO_EXIT="1"
else 
   if [ ! -e  "$REFTOP" ]
   then 
       echo "Topology file $REFTOP does not exist"
       DO_EXIT="1"
   fi
   echo "Using REFTOP=$REFTOP"
fi

# Directory containing trajectory files
if [ -z "$TRAJDIR" ]
then
   echo "Trajectory directory not defined, but needed." 
   DO_EXIT=1
else 
   if  ls $TRAJDIR 1> /dev/null 2>&1 
   then
      echo "Using TRAJDIR=$TRAJDIR"
   else 
       echo "Traj directory $TRAJDIR does not exist"
       if [ $DO_CC_MLA -le 1 ];then  DO_EXIT=1; fi 
   fi
fi
if [ -z "$PREFIX_MDCRD" ]
then
   echo "Prefix of mdcrd files not defined, but needed." 
   DO_EXIT=1
else 
   echo "Assuming PREFIX_MDCRD=$PREFIX_MDCRD"
fi
if [ -z "$SUFFIX_MDCRD" ]
then
   echo "Suffix of mdcrd files not defined, but recommended." 
   SUFFIX_MDCRD=".mdcrd"
   echo "Attempting SUFFIX_MDCRD=$SUFFIX_MDCRD"
else 
   echo "Assuming SUFFIX_MDCRD=$SUFFIX_MDCRD"
fi

fi     # endif of DO_CC_MLA -lt 1

# TEMP Directories
if [  -z "$SCRATCH" ] 
then
   echo "SCRATCH directory not defined, but needed." 
   DO_EXIT=1
else 
   if [ ! -e  "$SCRATCH" ]
   then 
       echo "SCRATCH=$SCRATCH does not exist"
       DO_EXIT=1
   else 
       TT=$(date +%N)
       TMPDIR=${SCRATCH}/TMPDIR_${TT}
       echo "Created $TMPDIR temporal directory"
       echo "Warning: /dev/shm (RAM) will be also used!"
   fi
fi

# Number of processors
if [ -z "$NPROCS" ]
then
   export NPROCS=$(cat /proc/cpuinfo | grep -c processor)
   echo "Using all $NPROCS available processors"
else
   echo "Using $NPROCS processors as predefined "
fi
export OMP_NUM_THREADS=$NPROCS 
export OMP_STACKSIZE="2G"    
 
# The following variables remain normally unchanged
# PERCEN 
if [ -z "$PERCEN" ]
then
   echo "Processing the whole data set"
   PERCEN="0"
else
   echo "Processing $PERCEN % of the data set"
fi

# OFFSET   
if [ -z "$OFFSET" ]
then
   echo "Processing all snapshots"
   OFFSET="1"
else
   echo "Processing snapshots with OFFSET=$OFFSET"
fi

# NINTERVAL
if [ -z "$NINTERVAL" ]
then
   NINTERVAL="20"
   echo "Convergence plot using ${NINTERVAL} points"
else
   echo "Convergence plot using ${NINTERVAL} points"
fi

# TEMPERATURE
# Concentration Parameter
if [ -z "$TEMPERATURE" ]
then
   TEMPERATURE="300"
fi
echo "TEMPERATURE (K)=$TEMPERATURE"

# Concentration Parameter
if [ -z "$KPARAM" ]
then
   KPARAM="0.50"
fi
echo "Concentration parameter K=$KPARAM"

# Maximum number of conformations (3 recommended)
if [ -z "$MAXNUMCONF" ]
then
   MAXNUMCONF="3"
fi
echo "Maximum number of conformations per torsion=$MAXNUMCONF"

# VERBOSE option for CCMLA
if [ -z "$VERBOSE" ]
then
   VERBOSE=""
else
   VERBOSE="-verbose"
fi

# Cutoffs for Composite Calcs
if [ -z "$CUTS2F1" ]
then
    CUTS2F1="0.01"
fi
if [ -z "$CUTS2F2" ]
then
    CUTS2F2="0.10"
fi

# DISTANCE MATRIX CALCULATION (ON/OFF)
DISTMAT="ON" 

if [ $DO_EXIT -eq 1 ]; then  echo 'Sorry, cannot continue'; exit; fi

##########################################################################3
# Run chicken run 

#WORK directory
WORKDIR=$PWD

if [ "${DO_CC_MLA}" -le 1 ] && [ "${DO_S2}" -le 1 ]  
then

# TRAJECTORY 
i=0
declare -a TRAJ=""
for trajfile in $(ls ${TRAJDIR}/${PREFIX_MDCRD}?${SUFFIX_MDCRD}  )
do
   let "i=$i+1"
   TRAJ["$i"]=$trajfile 
done
for trajfile in $(ls ${TRAJDIR}/${PREFIX_MDCRD}??${SUFFIX_MDCRD} )
do
   let "i=$i+1"
   TRAJ["$i"]=$trajfile 
done
for trajfile in $(ls ${TRAJDIR}/${PREFIX_MDCRD}???${SUFFIX_MDCRD}  )
do
   let "i=$i+1"
   TRAJ["$i"]=$trajfile 
done
NTRAJ="$i"

echo "After doing ls $TRAJDIR/${PREFIX_MDCRD}*${SUFFIX_MDCRD} .... we get"
echo "number of trajectory files = $NTRAJ "  

if [ $NTRAJ -eq 0 ]
then 
    echo "No trajectory files"; exit;
fi

# Working in the local scratch may improve efficiency 
mkdir $TMPDIR
cd    $TMPDIR 

# Number of frames in the trajectory file
for ((index=1; index <= NTRAJ; index++))
do

NFRAMES[$index]=$($NCDUMP -h  ${TRAJ[$index]} | grep 'frame = ' | sed -e 's/.*[^0-9]\([0-9]\+\)[^0-9]*$/\1/')

INIT[$index]="1"

echo "Trajectory: ${TRAJ[$index]} Init= ${INIT[$index]} nframes= ${NFRAMES[$index]}  offset= ${OFFSET}" 

done

# Preparing cpptraj input 

$PYTHON/python3 $GETTOR ${GETTOR_OPT} -sel ${MOL_MASK/:/} $REFTOP  > temp.in 

for ((index=1; index <= NTRAJ; index++))
do
   echo  "trajin ${TRAJ[$index]}  ${INIT[$index]}  ${NFRAMES[$index]}  10"   >>  matrix.in 
done

sed  -n '/matrix/,$ p' temp.in   >> matrix.in 

for ((index=1; index <= NTRAJ; index++))
do
  echo  "trajin ${TRAJ[$index]}  ${INIT[$index]}  ${NFRAMES[$index]}  ${OFFSET} "    >>  torsion.in 
done

grep -v trajin temp.in | sed '/matrix/,$ d' >> torsion.in
rm -f  temp.in

# Running CPPTRAJ
if [ $NPROCS -gt "1" ]
then 

# Splitting the torsion.in file 
grep -v 'trajin' torsion.in | grep -v '#' > temp.in
nlines=$(cat temp.in | wc -l)
if [ $nlines -gt $OMP_NUM_THREADS ]
then
   let " nlines_per_task = ( $nlines / $OMP_NUM_THREADS ) "
else
   nlines_per_task=1
fi
split -l $nlines_per_task  -d  temp.in   temp_X
rm -f input.torsion 
itor=0
for file in $(ls temp_X*) 
do
  let "itor=$itor+1"
  for ((index=1; index <= NTRAJ; index++))
  do
    echo  "trajin ${TRAJ[$index]}  ${INIT[$index]}  ${NFRAMES[$index]}  ${OFFSET} "    >>  torsion_${itor}.in
  done
  cat $file >> torsion_${itor}.in
  echo " ${CPPTRAJ}  $REFTOP  <  torsion_${itor}.in > torsion_${itor}.out " >>input.torsion 
  rm -f $file 
done 

cat input.torsion | $PARALLEL -t -j $NPROCS
    
else

echo "Running cpptraj torsion/pucker ..."
$CPPTRAJ  $REFTOP < torsion.in > torsion_cpptraj.out 

fi

# CPPTRAJ CPU TIME 
CPU_TOR=$(grep 'Total execution time:' torsion*out | awk '{print SUM+=$5}' | tail -1)

# Processing dihedral data 
ls d????.dat > LISTA

rm -f TASK.sh
echo "Running cencalc_prep individual tasks (time consuming)"
for file in $(cat LISTA)
do
  dihed=${file/.dat/}
  echo  "$PREP -k $KPARAM  -maxconf $MAXNUMCONF -nocut -ag yes -plot yes -wrbigmat ${dihed}.bm ${file} > ${dihed}.out" >> TASK.sh
done

chmod 755 TASK.sh
export OMP_NUM_THREADS=1
cat TASK.sh | $PARALLEL -j $NPROCS

CPU_PREP=$(grep ' CPU-TIME ' d????.out | awk '{print SUM+=$4}' | tail -1)
paste d????.bm | sed 's/\t//g' > BIGMAT.dat
grep -h 'Minima' d????.out  > MINIMA.dat

export OMP_NUM_THREADS=$NPROCS

if [ $DISTMAT == "ON" ]
then

  echo "Running cpptraj matrix job ..."
  $CPPTRAJ_OMP   $REFTOP < matrix.in > matrix_cpptraj.out 
  CPU_DISTMAT=$(grep 'Total execution time:' matrix_cpptraj.out | awk '{print $5}' )

  #Preparing CENCALC data
  echo "Running cencalc_prep global ..."
  $PREP -k $KPARAM  -maxconf $MAXNUMCONF -ag yes -plot yes -rdbigmat BIGMAT.dat  d????.dat > cencalc_prep.out 

else

#Preparing CENCALC data
  echo "Running cencalc_prep global ..."
  $PREP -k $KPARAM  -maxconf $MAXNUMCONF -ag yes -plot yes -nocut  -rdbigmat BIGMAT.dat  d????.dat > cencalc_prep.out 

fi

# Copy important stuff in WORKDIR
cp cencalc_prep.out MATRIX.dat MINIMA.dat reduced_dist_matrix.dat distance_matrix.dat \
   matrix_cpptraj.out atoms_in_tor.info torsion.in d*.png    $WORKDIR/

echo "CPU TIME INFO : (s)" >  $WORKDIR/CPU_TIME.dat
echo "===================" >> $WORKDIR/CPU_TIME.dat
echo "Cpptraj-dihedral : $CPU_TOR " >> $WORKDIR/CPU_TIME.dat
echo "Cpptraj-distmat  : $CPU_DISTMAT " >> $WORKDIR/CPU_TIME.dat
echo "Cencalc-Prep     : $CPU_PREP " >> $WORKDIR/CPU_TIME.dat

cd $WORKDIR/

rm -r -f $TMPDIR 

fi 

# Preparing cencalc executions
TT=$(date +%N)
TMPSHM=/dev/shm/TMPDIR_${TT}
mkdir $TMPSHM 
NFRAMES=$(wc -l MATRIX.dat | awk '{print $1}')

#  Optionaly, only a fraction of the data in MATRIX.dat is used
if [ "$PERCEN" -gt "100" ]
then
   echo 'PERCENTAGE greater than 100' 
   echo 'Using all data !'
   PERCEN=0
fi
if [ "$PERCEN" -lt "-100" ]
then
   echo 'PERCENTAGE greater than 100' 
   echo 'Using all data !'
   PERCEN=0
fi
if [ $PERCEN -eq 0 ] 
then 
   MATRIX_FILE="MATRIX.dat"
   SUFFIX_TAB=""
else
   echo "Using only $PERCEN % of the available data"
   echo "  > 0  --> From the begining"
   echo "  < 0  --> From the end"
   if [ "$PERCEN" -gt 0 ]
   then
       let  " NFRAMES_use  =  ( $NFRAMES *  $PERCEN ) / 100  "
       head -${NFRAMES_use}  MATRIX.dat > $TMPSHM/MATRIX_${PERCEN}.dat
       echo "Using only the first $NFRAMES_use frames"
       SUFFIX_TAB="_first_${PERCEN}"
   else
       let  " NFRAMES_use  =  ( $NFRAMES * ( - $PERCEN ) ) / 100  "
       tail -${NFRAMES_use}  MATRIX.dat > $TMPSHM/MATRIX_${PERCEN}.dat
       echo "Using only the last $NFRAMES_use frames"
       SUFFIX_TAB="_last_${PERCEN/-/}"
   fi 
   MATRIX_FILE="$TMPSHM/MATRIX_${PERCEN}.dat"
fi 

#Running cencalc 
echo "Running cencalc S1 ..."
NFRAMES=$(wc -l $MATRIX_FILE | awk '{print $1}')
let "NSTEP = ${NFRAMES} / ${NINTERVAL}"
$CENCALC -s1  -c -1 -data $MATRIX_FILE  -ns $NSTEP $NFRAMES $NSTEP -t s1_plot.tab${SUFFIX_TAB} > s1.out${SUFFIX_TAB}

CPU_S1=$(grep ' Total ' s1.out  | grep CPU | awk '{print $4}') 
echo "Cencalc-S1       : $CPU_S1 " >> $WORKDIR/CPU_TIME.dat

if [ "$DO_CC_MLA" -eq 1 ]  || [ "$DO_CC_MLA" -eq 2 ]
then

   if [ ! -e $MATRIX_FILE ];  then  echo "DO_CC_MLA=1,2 but $MATRIX_FILE does not exist"; exit; fi
   NFRAMES=$(wc -l $MATRIX_FILE | awk '{print $1}')
   let "NSTEP = ${NFRAMES} / ${NINTERVAL}"

   if [ "$PERCEN" -eq 0 ]
   then 
       cp MATRIX.dat $TMPSHM/
       MATRIX_FILE="$TMPSHM/MATRIX.dat"
   fi

   rm -f ccmla_plot.tab${SUFFIX_TAB}
   touch ccmla_plot.tab${SUFFIX_TAB} 
   cuttoff_plot=""

   if [ ! -e s1_plot.tab${SUFFIX_TAB} ] 
   then 
      $CENCALC -s1 -c -1 -data $MATRIX_FILE  -ns $NSTEP $NFRAMES $NSTEP -t s1_plot.tab${SUFFIX_TAB} > s1.out${SUFFIX_TAB}
   fi

   if [ $DO_COMP_CC_MLA -eq  0 ]
   then 

        for cutoff in $(echo $CUTOFF)
        do
            echo "Running cencalc CCMLA with cutoff=${cutoff} "
            if [ ! -e  s_ccmla_${cutoff}.out${SUFFIX_TAB} ]
            then
            $CENCALC $VERBOSE -ccmla  -c ${cutoff} -data $MATRIX_FILE   -dist  reduced_dist_matrix.dat \
             -ns $NSTEP  $NFRAMES  $NSTEP  -t s_ccmla_${cutoff}.tab${SUFFIX_TAB} >  s_ccmla_${cutoff}.out${SUFFIX_TAB}
            else
               echo "File: s_ccmla_${cutoff}.out${SUFFIX_TAB} exists. If you wish to recompute it, please do rm -f s_ccmla_${cutoff}.* and execute run_sconform.sh again"
            fi
            paste ccmla_plot.tab${SUFFIX_TAB} s_ccmla_${cutoff}.tab${SUFFIX_TAB}  > tmp; mv tmp ccmla_plot.tab${SUFFIX_TAB}
            cutoff_plot="${cutoff_plot} ${cutoff} ,"
            CPU_CC=$(grep ' Total ' s_ccmla_${cutoff}.out${SUFFIX_TAB}  | grep CPU | awk '{print $4}') 
            echo "Cencalc-CCMLA ${cutoff}  : $CPU_CC " >> CPU_TIME.dat
        done
        cutoff_plot=${cutoff_plot::-1}
        sed -i 's/\t/ /g' ccmla_plot.tab${SUFFIX_TAB}

   else

       touch ccmla_plot.tab${SUFFIX_TAB} 
       cuttoff_plot=""
       for cutoff in $(echo $CUTOFF)
       do
           echo "Running cencalc composite CCMLA cutoff=${cutoff} "
           if [ ! -e  s_ccmla_comp_${cutoff}.out${SUFFIX_TAB} ]
           then
              $CENCALC $VERBOSE  -ccmla  -c ${cutoff} -data $MATRIX_FILE   -dist  reduced_dist_matrix.dat \
             -s2filt $CUTS2F1  $CUTS2F2  -prs2mat s_ccmla_comp.s2mat \
            -ns $NSTEP  $NFRAMES  $NSTEP  -t s_ccmla_comp_${cutoff}.tab${SUFFIX_TAB} >  s_ccmla_comp_${cutoff}.out${SUFFIX_TAB}
           else
               echo "File: s_ccmla_comp_${cutoff}.out${SUFFIX_TAB} exists. If you wish to recompute it, please do rm -f s_ccmla_comp_${cutoff}.* and execute run_sconform.sh again"
           fi
           paste ccmla_plot.tab${SUFFIX_TAB} s_ccmla_comp_${cutoff}.tab${SUFFIX_TAB}  > tmp; mv tmp ccmla_plot.tab${SUFFIX_TAB}
           cutoff_plot="${cutoff_plot} ${cutoff} ,"
           CPU_CC=$(grep ' Total ' s_ccmla_comp_${cutoff}.out${SUFFIX_TAB}  | grep CPU | awk '{print $4}') 
           echo "Cencalc-CCMLA comp $CUTS2F1  $CUTS2F2  ${cutoff} : $CPU_CC " >> CPU_TIME.dat
       done
       cutoff_plot=${cutoff_plot::-1}
       sed -i 's/\t/ /g' ccmla_plot.tab${SUFFIX_TAB}

fi

cp $APTAMD/ENTROPY/sconform_plot.py .

sed -i  "s/DUMMY_CUTOFF/${cutoff_plot}/g" sconform_plot.py
sed -i  "s/DUMMY_METHOD/ccmla/g" sconform_plot.py
sed -i  "s/DUMMY_TEMP/${TEMPERATURE}/g" sconform_plot.py
sed -i  "s/DUMMY_SUFFIX_TAB/${SUFFIX_TAB}/g" sconform_plot.py
sed -i  "s/DUMMY_COMPOSITE/${DO_COMP_CC_MLA}/g" sconform_plot.py

$PYTHON/python3 sconform_plot.py

fi


if [ "$DO_S2" -eq 1 ]  || [ "$DO_S2" -eq 2 ]
then

   if [ ! -e $MATRIX_FILE ];  then  echo "DO_S2=1,2 but $MATRIX_FILE does not exist"; exit; fi
   NFRAMES=$(wc -l $MATRIX_FILE | awk '{print $1}')
   let "NSTEP = ${NFRAMES} / ${NINTERVAL}"

   if [ "$PERCEN" -eq 0 ]
   then 
       cp MATRIX.dat $TMPSHM/
       MATRIX_FILE="$TMPSHM/MATRIX.dat"
   fi

   rm -f s2_plot.tab${SUFFIX_TAB}
   touch s2_plot.tab${SUFFIX_TAB} 
   cuttoff_plot=""

   if [ ! -e s1_plot.tab${SUFFIX_TAB} ]; then $CENCALC -s1 -c -1 -data $MATRIX_FILE  -ns $NSTEP $NFRAMES $NSTEP -t s1_plot.tab${SUFFIX_TAB} > s1.out${SUFFIX_TAB}; fi

        for cutoff in $(echo $CUTOFF)
        do
            echo "Running cencalc S2 with cutoff=${cutoff} "
            if [ ! -e  s2_${cutoff}.out${SUFFIX_TAB} ]
            then
                $CENCALC $VERBOSE -s2 -shuffle  -c ${cutoff} -data $MATRIX_FILE   -dist  reduced_dist_matrix.dat \
                -ns $NSTEP  $NFRAMES  $NSTEP  -t s2_${cutoff}.tab${SUFFIX_TAB} >  s2_${cutoff}.out${SUFFIX_TAB}
            else
               echo "File: s2_${cutoff}.out${SUFFIX_TAB} exists. If you wish to recompute it, please do rm -f s2_${cutoff}.* and execute run_sconform.sh again"
            fi
            paste s2_plot.tab${SUFFIX_TAB} s2_${cutoff}.tab${SUFFIX_TAB}  > tmp; mv tmp s2_plot.tab${SUFFIX_TAB}
            cutoff_plot="${cutoff_plot} ${cutoff} ,"
            CPU_S2=$(grep ' Total ' s2_${cutoff}.out${SUFFIX_TAB}  | grep CPU | awk '{print $4}') 
            echo "Cencalc-S2 ${cutoff}  : $CPU_S2 " >> CPU_TIME.dat
        done
        sed -i 's/\t/ /g' s2_plot.tab${SUFFIX_TAB}
     
    cp $APTAMD/ENTROPY/sconform_plot.py .

    sed -i  "s/DUMMY_CUTOFF/${cutoff_plot}/g" sconform_plot.py
    sed -i  "s/DUMMY_METHOD/s2/g" sconform_plot.py
    sed -i  "s/DUMMY_TEMP/${TEMPERATURE}/g" sconform_plot.py
    sed -i  "s/DUMMY_SUFFIX_TAB/${SUFFIX_TAB}/g" sconform_plot.py
    sed -i  "s/DUMMY_COMPOSITE/0/g" sconform_plot.py

    $PYTHON/python3 ./sconform_plot.py



fi

rm -r -f $TMPSHM 

#=======================================================================================
 

