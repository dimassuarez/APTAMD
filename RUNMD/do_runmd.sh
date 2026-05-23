#!/bin/bash 

# do_runmd.sh : drives the execution of GAMD and MD simulations

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
if [ "$#" -eq 0 ]; then more $APTAMD/DOC/do_runmd.txt ; exit; fi

source $APTAMD/ENV/aptamd_env.sh

# Checking options
extension="${1##*.}"
if [ $extension == "src" ]
then
    echo "Sourcing $1"
    source $1
else
    MOL=$1
    MD_TYPE=$2
fi

if [ -z $MOL ]; then more $APTAMD/DOC/do_runmd.txt ; exit; fi
if [ -z $MD_TYPE ]; then more $APTAMD/DOC/do_runmd.txt ; exit; fi
if [ $MD_TYPE != "MD" ] && [ $MD_TYPE != "GAMD" ]
then 
     echo "MD_TYPE=$MD_TYPE, but must be MD or GAMD."
     exit
fi

if [ -n "$PBS_ENVIRONMENT" ] ; then
  NPROCS=$(cat $PBS_NODEFILE | wc -l)
fi

if [ -z "$NPROCS" ]
then
   NPROCS=$(cat /proc/cpuinfo | grep -c processor)
   echo "Using $NPROCS available processors"
else
   echo "Using $NPROCS processors "
fi

if [ -z $DO_PROD ]; then DO_PROD="YES"; fi 
if [ -z $USE_GPU ]; then USE_GPU="YES"; fi 
if [ -z $INIT ]; then INIT=0; fi 
if [ -z $NJOBS ]; then NJOBS=50; fi 
if [ -z $PEEL ]; then PEEL=14; fi 
if [ -z $PEELONLY ]; then PEELONLY="NO"; fi 
if [ -z $SIEVE ]; then SIEVE=20 ; fi 
if [ -z $MULTI_MD ]; then MULTI_MD=1 ; fi 

# Filenames 
REFTAR="$APTAMD/RUNMD/inputfiles_${MD_TYPE}.tar" 
TOPOLOGY="${MOL}.top"
MD_DIR="${MOL}_${MD_TYPE}"
if [ ! -e $MD_DIR ]; then echo "$MD_DIR does not exist!"; exit; fi

# Entering  MD directory (1.EDITION must be prepared in advance)
cd $MD_DIR 
if [ ! -e 1.EDITION ]; then echo "1.EDITION does not exist"; exit; fi 
if [ -e 2.RELAX_SOLVENT ]; then echo "2.RELAX_SOLVENT is present. Better to remove it manually."; exit; fi 
if [ -e 3.THERMALIZATION ]; then echo "3.THERMALIZATION  is present. Better to remove it manually."; exit; fi 
if [ -e 4.PRESSURIZATION ] ; then echo "4.PRESSURIZATION is present. Better to remove it manually."; exit; fi 
if [ -e 5.PRODUCTION ]; then echo "5.PRODUCTION is present. Better to remove it manually."; exit; fi 

# rm -r -f 2.RELAX_SOLVENT 3.THERMALIZATION 4.PRESSURIZATION 5.PRODUCTION

tar xvf $REFTAR

cd 1.EDITION/
if [ ! -e $TOPOLOGY ]; then echo "$TOPOLOGY does not exist"; exit; fi 
if [ ! -e ${MOL}.pdb ]; then echo "${MOL}.pdb does not exist"; exit; fi 
if [ ! -e ${MOL}.crd ]; then echo "${MOL}.crd does not exist"; exit; fi 

cd ../

if [ ! -e $TOPOLOGY ]; then ln -s  1.EDITION/$TOPOLOGY $TOPOLOGY; fi
if [ ! -e initial.pdb ]; then ln -s  1.EDITION/${MOL}.pdb initial.pdb ; fi
if [ ! -e initial.crd ]; then ln -s  1.EDITION/${MOL}.crd initial.crd ; fi

NAT=$(grep -c 'ATOM  ' initial.pdb) 

export NPROCS
export USE_GPU

# THERMALIZATION AND PRESSURIZATION 
cd 2.RELAX_SOLVENT
touch ../running_2.RELAX_SOLVENT
bash run_relax.sh ../$TOPOLOGY 
cp min_all.rst ../3.THERMALIZATION
rm -f ../running_2.RELAX_SOLVENT
   
cd ../3.THERMALIZATION
touch ../running_3.THERMALIZATION
bash run_therma.sh ../$TOPOLOGY 
cp therma_300.rst ../4.PRESSURIZATION
rm -f ../running_3.THERMALIZATION

cd ../4.PRESSURIZATION 
touch ../running_4.PRESSURIZATION
bash run_press.sh ../$TOPOLOGY 
rm -f ../running_4.PRESSURIZATION

# Once that pressurization is run, obtain  peeling info and auxiliary topologies
if [ $PEEL -gt 0 ] 
then 
# The mdcrd file is processed to generate PEEL.info and  ${MOL}_solutewat.top
  $APTAMD/RUNMD/peel_mdcrd_cpptraj.sh ../$TOPOLOGY press.mdcrd $PEEL 1
  mv ${MOL}_solutewat.top   ../
  cp ${MOL}_solutewat.info  PEEL.info
  mv ${MOL}_solutewat.info  ../
  rm -f press_solutewat.mdcrd 

else

echo ' 0  0  0  0 ' >  PEEL.info

fi

if [ $SIEVE -gt 1 ]
then 

# Only solute atoms topology
  $AMBERHOME/bin/parmed  ../$TOPOLOGY  <<EOF
strip ":WAT,Cl-,Na+"
parmout ../${MOL}_solute.top
go
EOF

fi

#Getting initial coordinates for production jobs
if [ "$MULTI_MD" -gt 50 ] 
then
   echo "Only half of the pressurization phase is considered"
   echo "for extracting the initial snaps of Multi-MD jobs"
   MULTI_MD=50
fi

ilast=100
let "ifirst=$ilast-  $MULTI_MD + 1  " 
$AMBERHOME/bin/cpptraj ../$TOPOLOGY  <<EOF
trajin    press.mdcrd   $ifirst  $ilast  1 
trajout   press_snap  ncrestart  time0 0.0 dt 0.0 keepext 
go
EOF


if [ $MULTI_MD -gt 1 ]
then

# Preparing Multi-MD production phase 
   echo "Preparing $MULTI_MD multi MD production directories"
   for ((I=1;I<=MULTI_MD;I++))
   do
   if [ $I -lt 10 ]
   then
         txt="0${I}"
   else
         txt=${I}
   fi
   cp -p -r ../5.PRODUCTION ../5.PRODUCTION_${txt} 
   if [ "$MD_TYPE" == "MD" ]  
   then 
        mv press_snap.${I}       ../5.PRODUCTION_${txt}/md_000.rst
   elif [ "$MD_TYPE" == "GAMD" ]
   then
        mv press_snap.${I}       ../5.PRODUCTION_${txt}/gamd.rst   
   fi
   cp PEEL.info             ../5.PRODUCTION_${txt}/
   cat <<EOF > ../5.PRODUCTION_${txt}/job_${MD_TYPE}_${txt}.sh
#!/bin/bash
export APTAMD=$APTAMD
export NPROCS=$NPROCS
export USE_GPU=$USE_GPU
export PEEL=$PEEL
export PEELONLY=$PEELONLY
export SIEVE=$SIEVE
./run_production_peel.sh ../${TOPOLOGY}  $INIT $NJOBS
EOF
   if [ "$MD_TYPE" == "GAMD" ]; then  echo "./run_gamd_eq.sh ../$TOPOLOGY" >> ../5.PRODUCTION_${txt}/job_${MD_TYPE}_${txt}.sh ; fi 
   echo "./run_production_peel.sh ../${TOPOLOGY}  $INIT $NJOBS"  >> ../5.PRODUCTION_${txt}/job_${MD_TYPE}_${txt}.sh
   chmod 755 ../5.PRODUCTION_${txt}/job_${MD_TYPE}_${txt}.sh

   echo "Submit job_${MD_TYPE}_${txt}.sh in ../5.PRODUCTION_${txt}/"

   done
   rm -r -f ../5.PRODUCTION
   cd ../

else

# Preparing and running production phase 
   if [ "$MD_TYPE" == "MD" ]  
   then 
        mv press_snap       ../5.PRODUCTION/md_000.rst
   elif [ "$MD_TYPE" == "GAMD" ]
   then
        mv press_snap       ../5.PRODUCTION/gamd.rst   
   fi
   cp PEEL.info  ../5.PRODUCTION/
   cat <<EOF > ../5.PRODUCTION/job_${MD_TYPE}.sh
#!/bin/bash
export APTAMD=$APTAMD
export NPROCS=$NPROCS
export USE_GPU=$USE_GPU
export PEEL=$PEEL
export PEELONLY=$PEELONLY
export SIEVE=$SIEVE
EOF
   if [ "$MD_TYPE" == "GAMD" ]; then  echo "./run_gamd_eq.sh ../$TOPOLOGY" >> ../5.PRODUCTION/job_${MD_TYPE}.sh ; fi 
   echo "./run_production_peel.sh ../${TOPOLOGY}  $INIT $NJOBS"  >> ../5.PRODUCTION/job_${MD_TYPE}.sh
   chmod 755 ../5.PRODUCTION/job_${MD_TYPE}.sh
   cd ../5.PRODUCTION
   if [ $DO_PROD == "YES" ] || [ $DO_PROD == "yes" ] 
   then 
      bash ./job_${MD_TYPE}.sh 
   fi 

fi
