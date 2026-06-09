#!/bin/bash

# Autodock 4 calculations

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
if [ "$#" -eq 0 ]; then more $APTAMD/DOC/do_adck.txt ; exit; fi

source $APTAMD/ENV/aptamd_env.sh

# Specific MGLTOOLS to be used
AUTODOCKTOOLS="$MGLTOOLS/MGLToolsPckgs/AutoDockTools"
PYTHONSH="$MGLTOOLS/bin/pythonsh"
if [ ! -e $PYTHONSH ] || [ ! -e $AUTODOCKTOOLS ]
then
   echo "$PYTHONSH and/or $AUTODOCKTOOLS do no exist, but required"
   exit
fi
if [ ! -e $ADCK/autogrid4 ] || [ ! -e $ADCK/autodock4  ]
then
   echo " autogrid4 and/or autodock4 not accessible, but required"
   exit
fi

# Checking options
extension="${1##*.}"
if [ $extension == "src" ]
then
    echo "Sourcing $1"
    source $1
else
# Input files containing the list of receptor pdb files (Absolute path is required) !
 RECEPTOR_PDB_LIST=$1 
 LIGAND_PDB_LIST=$2
 if [  -z "$RECEPTOR_PDB_LIST" ]; then echo 'Usage: do_adck.sh RECEPTOR_PDB_LIST LIGAND_PDB_LIST' ; exit; fi
 if [  -z "$LIGAND_PDB_LIST" ];   then echo 'Usage: do_adck.sh RECEPTOR_PDB_LIST LIGAND_PDB_LIST'; exit; fi
fi

# Num pf processors
if [ -z "$NPROCS" ]
then
   NPROCS=$(cat /proc/cpuinfo | grep -c processor)
   echo "Using NPROCS=$NPROCS available"
else
   echo "Using NPROCS=$NPROCS as predefined "
fi
# 
if [ -z "$TOP_PERCEN" ]    # Percentage of top docking poses to keep when cleaning up directories
then
   TOP_PERCEN=10
fi
if [ -z $KEEPSOLV_REC ]
then
    KEEPSOLV_REC="NO"
elif [ $KEEPSOLV_REC == "YES" ]
then 
    echo 'KEEPSOLV_REC=YES'
    echo 'Na+ ions and H2O molecules preserved in receptor coordinates'
fi 
if [ $KEEPSOLV_REC == "NO" ] ; then echo 'Counterions and H2O molecules are removed from receptor coordinates' ; fi

# For MM relaxation of Autodock poses
if [ -z $RELAXOPT ]
then
    RELAXOPT="GBSA"
fi
if [ -z $CUTOFF ]; then CUTOFF="10.0"; fi  # For relaxation. Probably no necessary to change.
if [ $RELAXOPT == "GBSA" ]
then
   echo "Autodock poses will be relaxed by MM/GBSA minimization prior to rescoring,"
   echo "applying positional harmonic restraints to heavy atoms > $CUTOFF  Ang of the contact region"
elif [ $RELAXOPT == "DIEL" ]
then
   echo "Autodock poses will be partially relaxed by dielectric-distance-dependant MM minimization prior to rescoring"
   echo "Cutoff for atom selection in receptor-ligand complex = $CUTOFF"
else
   echo "RELAXOPT=$RELAXOPT not valed. It must be GBSA or DIEL"
   exit
fi

# AUTODOCK OPTIONS (for rescoring)
if [ -z "$GRID" ]
then
   GRID="0.333"
   echo "Using GRID=${GRID}"
else
   echo "Using GRID=${GRID} as predefined "
fi

export OMP_NUM_THREADS=$NPROCS
# Flow control 
if [ -z "$DO_DOCKING" ]
then
    DO_DOCKING="YES"
else
   if [ $DO_DOCKING == "NO" ]; then echo "Omitting Docking Calculations. Doing Relax-Rescore-Analysis"; fi 
fi
if [ -z "$DO_RELAX" ]
then
    DO_RELAX="YES"
else
   if [ $DO_RELAX == "NO" ] ; then echo "Omitting Relax Calculations. Useful for debugging only"; fi 
fi
if [ -z "$TOP_PERCEN_PRESCORE" ]
then
   TOP_PERCEN_PRESCORE=33
fi
if [ -z "$PRESCORE" ]
then
    PRESCORE="NO"
    echo "Omitting Prescoring Calculations."
else
   if [ $PRESCORE == "NO" ] ; then echo "Omitting Prescoring Calculations."; fi 
   if [ $PRESCORE == "YES" ] ; then echo "Performing Prescoring Calculations."; fi 
   if [ $PRESCORE == "YES" ] ; then echo " Top $TOP_PERCEN_PRESCORE % prescore will be relaxed"; fi 
fi
if [ -z "$DO_RESCORE" ]
then
    DO_RESCORE="YES"
else
   if [ $DO_RESCORE == "NO" ] ; then echo "Omitting Rescore Calculations. Useful for debugging only"; fi 
fi
if [ -z "$DO_GRIDMAP" ]
then
    DO_GRIDMAP="YES"
else
   if [ $DO_GRIDMAP == "NO" ] ; then echo "Omitting Gridmap Calculations. Useful for debugging only"; fi 
fi
if [ -z "$RM_RECDIR" ]
then
    RM_RECDIR="YES"
else
   if [ $RM_RECDIR == "NO" ] ; then echo "Do not remove previs REC_I directories. Useful for debugging only"; fi 
fi

WORKDIR=$PWD

echo "WORKDIR=$PWD"
echo "No other working directory is made"

for file in receptor.pdbqt ligand.pdbqt receptor.gpf receptor_ligand.dpf complex.top ligand.off input_leap_ligand.src ligand.frag receptor.pdb ligand.pdb 
do
  if [ ! -e $file ]
  then
     echo "$file does not exist. This should be available from a "
     echo "previous execution of do_adck_templates.sh" 
     exit
  fi
done

# Clean up clean_up_REC_LIG_data.sh script from prev executions
rm -f clean_up_REC_LIG_data.sh
    
# Getting some dimensions 
nrec=$(cat $RECEPTOR_PDB_LIST | wc -l)
nlig=$(cat $LIGAND_PDB_LIST | wc -l)

nres_rec=$(grep 'ATOM ' receptor.pdbqt | tail -1 | awk '{print $5}')
nres_lig=$(grep 'ATOM  ' ligand.pdbqt | tail -1 | awk '{print $5}')

#  SANDER input for MM/GBSA relaxation 
#  For speed , cut=25.0 and rgbmax=15   Note that restraints apply
cat << EOF > sander_min_gbsa.inp
GBSA Minimization IGB = 5  
&cntrl  
 imin=1, maxcyc=500, ntmin=2, drms=0.002 
 ntb=0, ntf=1, ntc=1, ntpr=100, 
 igb=5, gbsa=1, saltcon=0.150, rgbmax=15.0, 
 intdiel=1.0, extdiel=80.0, 
 cut=25.0, nsnb=5000,  
 ntr=1,
 restraint_wt=25.0,
 restraintmask='DUMMY_RSTMASK',
 
 / 
EOF

cat  <<EOF  >  sander_min_belly.inp 
Distance dependent eps
&cntrl
 imin=1, ncyc=100,  maxcyc=500, ntmin=2, drms=0.02
 ntb=0, ntf=1, ntc=1, ntpr=100,
 cut=100.0, nsnb=5000,  igb=0,
 ibelly=1,
 bellymask='DUMMY_BELLY',
 /
 &ewald
  eedmeth=5,
 /
EOF

IATCEN=$(grep '#IATCEN=' receptor.gpf | awk '{print $2}')
echo "IATCEN=$IATCEN"
IZATM=$(grep '#IZATM=' receptor.gpf | awk '{print $2}')
echo "IZATM=$IZATM"
if [ $IZATM -gt 0 ] && [ ! -e $WORKDIR/receptor.zatm ]
then
   echo "Z-potential defined in receptor.gpf, but no receptor.zatm file is available"
   exit
fi
if [ $IZATM -gt 0 ] && [ ! -e $WORKDIR/ligand.zatm ]
then
   echo "Z-potential defined in receptor.gpf, but no ligand.zatm file is available"
   exit
fi
if [ $IZATM -gt 0 ]
then
  echo "Z-potential defined in receptor.gpf/ligand.zatm: Autodock poses must be relaxed"
  echo "by means of dielectric-distance-dependant MM minimization."
  cp sander_min_belly.inp sander_min.inp 
  BELLYMASK="!@"
  ZATM_REC=$(grep ZATM $WORKDIR/receptor.zatm | awk '{print $2}') 
  I=0
  for IDATM in $ZATM_REC
  do  
    let "I=$I+1"
    if [ $I -eq 1 ]
    then 
       BELLYMASK="${BELLYMASK}${IDATM}" 
    else
       BELLYMASK="${BELLYMASK},${IDATM}" 
    fi
  done

# Maybe a problem here if FF19 solvent is used

  pdb_rec=$(head -1 $RECEPTOR_PDB_LIST)
  nat_rec=$(grep -c  'ATOM  \|HETATM' $pdb_rec)
  ZATM_LIG=$(grep ZATM $WORKDIR/ligand.zatm | awk '{print $2}') 
  for IDATM in $ZATM_LIG
  do  
    let "IAT=$nat_rec + $IDATM"
    BELLYMASK="${BELLYMASK},${IAT}" 
  done
  sed -i "s/DUMMY_BELLY/${BELLYMASK}/" sander_min.inp

elif [ $RELAXOPT == "GBSA" ]
then
   cp sander_min_gbsa.inp sander_min.inp 
   NRESREC=$(grep 'ATOM  \|HETATM' receptor.pdb | awk '{printf("%3s %6i \n", $4,$5)}' |sort -n -k 2 | uniq | wc -l)
   NRESLIG=$(grep 'ATOM  \|HETATM' ligand.pdb   | awk '{printf("%3s %6i \n", $4,$5)}' |sort -n -k 2 | uniq | wc -l)
   let "IRES=$NRESREC+1" 
   let "JRES=$NRESREC+$NRESLIG" 
   RSTMASK="(( (:1-${NRESREC})\&(:${IRES}-${JRES}>@${CUTOFF}) ) | ((:${IRES}-${JRES})\&(:1-${NRESREC}>@${CUTOFF})) )\&(!@H=)"
   sed -i "s/DUMMY_RSTMASK/${RSTMASK}/" sander_min.inp
   
elif [ $RELAXOPT == "DIEL" ]
then
   # Probably not adequate for DNA complexes 
   cp sander_min_belly.inp sander_min.inp 
   NRESREC=$(grep 'ATOM  \|HETATM' receptor.pdb | awk '{printf("%3s %6i \n", $4,$5)}' |sort -n -k 2 | uniq | wc -l)
   NRESLIG=$(grep 'ATOM  \|HETATM' ligand.pdb   | awk '{printf("%3s %6i \n", $4,$5)}' |sort -n -k 2 | uniq | wc -l)
   let "IRES=$NRESREC+1" 
   let "JRES=$NRESREC+$NRESLIG" 
   BELLYMASK="( (:1-${NRESREC})\&(:${IRES}-${JRES}<@${CUTOFF}) ) | ((:${IRES}-${JRES})\&(:1-${NRESREC}<@${CUTOFF})) | (:${IRES}-${JRES}\&@H=)"
   sed -i "s/DUMMY_BELLY/${BELLYMASK}/" sander_min.inp
fi

if [ ${DO_DOCKING} == "YES" ]
then 

irec=0
rm -f TASK.sh 


for lfile in $(cat $RECEPTOR_PDB_LIST)
do
   if [ -e $lfile ]
   then 
      file=$lfile
   elif [ -e $WORKDIR/$lfile ]
   then
      file=$WORKDIR/$lfile
   else
      echo "$lfile listed in $RECEPTOR_PDB_LIST does not exist!"
      continue
   fi 
   let "irec=$irec+1"
   if [ $RM_RECDIR == "YES" ] && [ $DO_GRIDMAP == "YES" ] &&  [ -d REC_${irec} ]; then rm -r -f REC_${irec} ; fi  
   if [ ! -d REC_${irec} ]; then mkdir REC_${irec} ; fi 
   if [ $KEEPSOLV_REC == "NO" ]
   then
       grep -v 'WAT\|HOH' $file | grep -v 'Na+' | grep -v 'K+' | grep -v 'Cl-' > temp.pdb
   else
       grep -v 'K+\|Cl-\|EPW' $file  > temp.pdb
   fi
   $TOOLS/center_geom -i temp.pdb -o REC_${irec}/rec_${irec}.pdb
   rm -f temp.pdb 
   grep -v PDBQT_PQR_MATCH receptor.pdbqt >  REC_${irec}/tmp.pdbqt
   $TOOLS/merge_pdbqt  REC_${irec}/tmp.pdbqt REC_${irec}/rec_${irec}.pdb  REC_${irec}/receptor.pdbqt
   rm -f REC_${irec}/tmp.pdbqt
   cp ${WORKDIR}/receptor.gpf REC_${irec}/
   if [ ${IATCEN} -gt 0 ]
   then
      CEN_LINE=$(grep 'ATOM  \|HETATM' REC_${irec}/rec_${irec}.pdb | sed -n  "${IATCEN},${IATCEN}p")
      CEN_REC=$(echo $CEN_LINE |  awk '{printf("%f %f %f\n",$6,$7,$8)}')
      sed -i "s/DUMMY_IATCEN_COORD/${CEN_REC}/" REC_${irec}/receptor.gpf
      echo "# Grid centered at $CEN_LINE"  >>   REC_${irec}/receptor.gpf
   fi
   if [ ${IZATM} -gt 0 ]
   then
      ZATM_REC=$(grep ZATM $WORKDIR/receptor.zatm | awk '{print $2}') 
      for IDATM in $ZATM_REC
      do  
         grep 'ATOM  \|HETATM' REC_${irec}/rec_${irec}.pdb | sed -n  "${IDATM},${IDATM}p" | \
         awk '{printf("%f %f %f\n",$6,$7,$8)}' >> tmp_zatm.xyz
      done
      CEN_ZATM=$(awk '{x+=$1;y+=$2;z+=$3;printf("%f %f %f\n",x/NR,y/NR,z/NR)}' tmp_zatm.xyz | tail -1)
      echo "ZPOT-CENTER coordinates=${CEN_ZATM}" 
      sed -i "s/DUMMY_ZPOT_CENTER/${CEN_ZATM}/" REC_${irec}/receptor.gpf
      rm -f tmp_zatm.xyz 
   fi
   echo "cd ${WORKDIR}/REC_${irec}; $ADCK/autogrid4 -p receptor.gpf -l receptor.glg" >> TASK.sh 
done 

if [ $DO_GRIDMAP == "YES" ]   # NO only for debugging
then
   echo "Running initial AUTOGRID as parallel task in $PWD"
   cat TASK.sh  | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS
fi

# Double loop over receptor and ligand files 
rm -f TASK.sh 
for ((irec=1;irec<=nrec;irec++))
do
   cd REC_${irec}
   ilig=0
   for lfile in  $(cat $WORKDIR/$LIGAND_PDB_LIST)
   do
     if [ -e $lfile ]
     then
         file=$lfile
     elif [ -e $WORKDIR/$lfile ]
     then
        file=$WORKDIR/$lfile
     else
        echo "$lfile listed in $LIGAND_PDB_LIST does not exist!"
        continue
     fi
     let "ilig=$ilig+1"
     if [ -e LIG_${ilig} ]; then rm -r -f LIG_${ilig} ; fi
     mkdir LIG_${ilig}
     cd LIG_${ilig}
     echo $file
     $TOOLS/center_geom -i $file -o lig_${ilig}.pdb
     $TOOLS/merge_pdbqt ../../ligand.pdbqt lig_${ilig}.pdb ligand.pdbqt
     sed -i '/COMMENT/d' ligand.pdbqt
#    Same map files 
     for maps in $(ls ../receptor*map*)
     do
        ln -s $maps $(basename $maps)
     done
     echo "cd ${WORKDIR}/REC_${irec}/LIG_${ilig}; $ADCK/autodock4  -p ../../receptor_ligand.dpf -l receptor_ligand.dlg"  >> $WORKDIR/TASK.sh
     cd ../
   done
   cd ../
done

cd $WORKDIR

echo "Running AUTODOCK as parallel tasks in $PWD"
cat TASK.sh  | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS

fi    #  END of if DO_DOCKING

cd $WORKDIR

# Triple loop over receptor, ligand and pose files 
# to get preliminay scores before relaxaion 

#  Note that for a given Autodock run: we get only energies and
#  coordinates of the cluster representatives generated by the
#  analysis command in Autodock

rm -f PRESCORE*.dat

for ((irec=1;irec<=nrec;irec++))
do
    
rm -f TASK_extract_${irec}.sh
cat <<EOF > TASK_extract_${irec}.sh 
#!/bin/bash

   cd $WORKDIR/REC_${irec} 

   for ((ilig=1;ilig<=${nlig};ilig++))
   do
      cd LIG_\${ilig}
      declare -a energy=""
      energy=(\$(grep "Final Intermolecular Energy" receptor_ligand.dlg | grep -v DOCKED |  awk '{print \$7}'))
      npose=\${#energy[@]}
      for ((ipose=1;ipose<=npose;ipose++))
      do 
         let "jpose=\$ipose-1"
         epose=\${energy["\$jpose"]}
         echo "complex_${irec}_\${ilig}_\${ipose}  \$epose "  >> ${WORKDIR}/PRESCORE_${irec}_\${ilig}.dat
      done 
      cd ../
   done
   cd ../
EOF

done

cd $WORKDIR
echo "Running extraction tasks in $PWD"
rm -f TASK.sh 
for file in $(ls TASK_extract_*)
do
    chmod 755 $file
    echo "$WORKDIR/${file}" >> TASK.sh
done 
cat TASK.sh  | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS
rm -f TASK_extract_*

cat PRESCORE_*.dat > PRESCORE.dat
sed -i 's/+/ /g' PRESCORE.dat
sort -n -k 2 PRESCORE.dat > tmp; mv tmp PRESCORE.dat
rm -f PRESCORE_*.dat

if [ ${PRESCORE} == "YES"  ]   # Chopping PRESCORE file only when PRESCORE=YES
then 

NPOSES=$(cat PRESCORE.dat | wc -l)
let "NCLEAN=$NPOSES * (100 - $TOP_PERCEN_PRESCORE ) / 100" 
let "NKEEP=$NPOSES - $NCLEAN" 
if [ $NCLEAN -gt 0 ]  
then 
  echo "Only the top $NKEEP poses out of total of  $NPOSES will be relaxed and rescored"
  head -${NKEEP} PRESCORE.dat  > tmp; mv tmp PRESCORE.dat 
fi 

fi


# Check if TORSDOF is active 
TORSDOF=$(grep TORSDOF ligand.pdbqt | awk '{print $NF}') 
echo "TORSDOF=$TORSDOF"
if [ $TORSDOF -gt 0 ]
then 
   echo "TORSDOF > 0 ---> Sorting PDBQT atoms into PQR ordering " 
   iord=($(grep MATCH ligand.pdbqt | sort -n -k 4 | awk '{print $(NF-1)}'))
   nord=${#iord[@]}
fi


if [ $DO_RELAX == "YES" ]               # BEGIN DO_RELAX IF 
then 

# Triple loop over receptor , ligand and pose files to extract and relax poses
for ((irec=1;irec<=nrec;irec++))
do

#  Inner loop tasks are performed by independent scripts that will be run in 
#  parallel mode. Care is needed here to escape the local script variables

rm -f $WORKDIR/TASK_pose_${irec}.sh $WORKDIR/TASK_sander_${irec}.sh

cat <<EOF > $WORKDIR/TASK_pose_${irec}.sh
#!/bin/bash

declare -a iord=""
iord=(${iord[*]})

   cd $WORKDIR/REC_${irec}
   for ((ilig=1;ilig<=${nlig};ilig++))
   do
      cd LIG_\${ilig}
      if [ -e OPT_MM ]; then rm -r -f OPT_MM; fi 
      mkdir OPT_MM
      ligand_atoms=\$(grep "ATOM" ligand.pdbqt | wc -l)
      grep 'ATOM \|HETATM' receptor_ligand.dlg | grep -v 'DOCKED:\|INPUT' | split -d -l \$ligand_atoms - ./OPT_MM/dock_ligand
      cd OPT_MM/
      ln -s ../../rec_${irec}.pdb receptor.pdb
      ipose=0
      rm -f NPOSE.dat
      for file in \$(ls dock_ligand*)
      do 
          let "ipose=\$ipose+1"
          if [ "${PRESCORE}" == "YES" ]
          then 
               icheck=\$(grep -c "complex_${irec}_\${ilig}_\${ipose} " $WORKDIR/PRESCORE.dat )
          else
               icheck=1
          fi
          if [ \${icheck} -eq 1  ]
          then 

          echo "\$ipose" >> NPOSE.dat
          rm -f ligand.pdb 
          sed 's/DOCKED: ATOM/ATOM/' \$file > temp.pdb 
          if [ $TORSDOF -gt 0 ] 
          then 
             declare -a FLINES=""
             iline=0
             while read line
             do
                FLINE["\$iline"]=\$line
                let "iline=\$iline+1"
             done < temp.pdb
             let   "nlines=\$iline"
             if [ \$nlines -ne $nord ]
             then
                 echo " ERROR nlines=\$nlines, but nord=$nord"
                 exit
             fi
             rm -f temp.pdb 
             for ((i=1;i<=${nord};i++))
             do
                 let "j=\$i-1"
                 let "jline=\${iord["\$j"]} - 1 "
                 echo "\${FLINE["\$jline"]}" >> temp.pdb
             done
          fi
          $TOOLS/reorder < temp.pdb > ligand_noter.pdb 
          iprev=0
          rm -f ligand.pdb
          for natfrag in \$(cat $WORKDIR/ligand.frag)
          do
              let "ifirst=\$iprev+1"
              let "ilast=\$iprev+\$natfrag"
              sed -n "\${ifirst},\${ilast}p" ligand_noter.pdb >> ligand.pdb
              echo 'TER' >> ligand.pdb
              let "iprev=\$iprev+\$natfrag"
          done
          rm -f ligand_noter.pdb
          pose=pose_\${ipose}
          $AMBERHOME/bin/tleap -f $WORKDIR/input_leap_ligand.src > mlog 
          mv  -f ligand_edited.pdb ligand.pdb
          cat receptor.pdb  > receptor_ligand.pdb
          echo 'TER' >>                 receptor_ligand.pdb
          cat  ligand.pdb >> receptor_ligand.pdb
          $TOOLS/pdbcrd < receptor_ligand.pdb > \${pose}.crd
          echo "cd $WORKDIR/REC_${irec}/LIG_\${ilig}/OPT_MM;  \
          $AMBERHOME/bin/sander -O -i $WORKDIR/sander_min.inp -p $WORKDIR/complex.top -c \${pose}.crd -ref \${pose}.crd -r \${pose}.rst -o \${pose}.out -inf \${pose}.inf" >>$WORKDIR/TASK_sander_${irec}.sh

          fi # End of icheck IF

      done 
      cd ../../
    done
    cd ../
EOF

done

cd $WORKDIR
echo "Running pose extraction tasks in $PWD"
rm -f TASK.sh 
for file in $(ls TASK_pose_*)
do
    chmod 755 $file
    echo "$WORKDIR/${file}" >> TASK.sh
done 
cat TASK.sh  | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS

cd $WORKDIR
echo "Running SANDER relaxations as parallel tasks in $PWD"
rm -f TASK.sh 
for file in $(ls TASK_sander_*)
do
    cat $file >> TASK.sh 
done 
rm -f TASK_pose_* TASK_sander_*
cat TASK.sh  | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS

fi  # End of DO_RELAX IF

cd $WORKDIR

if [ $DO_RESCORE  == "YES"  ]
then 

# Triple loop over receptor, ligand and pose files to perform AUTODOCK rescore
for ((irec=1;irec<=nrec;irec++))
do

#  Inner loop tasks are performed by independent scripts that will be run in 
#  parallel mode. Care is needed here to escape the local script variables

rm -f  TASK_prep_${irec}.sh TASK_adck_${irec}.sh

cat <<EOF > TASK_prep_${irec}.sh
#!/bin/bash

   cd $WORKDIR/REC_${irec}
   for ((ilig=1;ilig<=${nlig};ilig++))
   do
      cd LIG_\${ilig}/OPT_MM
      npose=\$(cat NPOSE.dat | wc -l)
      for ipose in \$(cat NPOSE.dat) 
      do 
         pose=pose_\${ipose}
         if [ -e \${pose}_RESCORE ]; then rm -r -f \${pose}_RESCORE; fi
         mkdir \${pose}_RESCORE 
         $AMBERHOME/bin/cpptraj $WORKDIR/complex.top <<EOX  > mlog
trajin \${pose}.rst
strip  !:1-${nres_rec}
trajout \${pose}_RESCORE/receptor.pdb pdbatom
go
EOX
         $AMBERHOME/bin/cpptraj $WORKDIR/complex.top <<EOX > mlog
trajin \${pose}.rst
strip  :1-${nres_rec}
trajout \${pose}_RESCORE/ligand.pdb pdbatom
go
EOX
         cd \${pose}_RESCORE
# Grid size is now tuned for the particular pose geometry
         $TOOLS/merge_pdbqt $WORKDIR/receptor.pdbqt receptor.pdb receptor.pdbqt  > mlog
         $TOOLS/merge_pdbqt $WORKDIR/ligand.pdbqt ligand.pdb  temp.pdbqt  > mlog 
         $TOOLS/center_geom -noatm -i ligand.pdb -o ligand.size
         CEN_LIG=\$(grep 'COM  =' ligand.size |  awk '{printf("%f %f %f\n",\$4,\$5,\$6)}')
         DIM_LIG=\$(grep 'XYZ DIM' ligand.size |  awk '{printf("%f %f %f\n",\$5,\$6,\$7)}')
         SIZE=\$(echo \$DIM_LIG $GRID | awk '{ix=1.5*\$1/\$4; iy=1.5*\$2/\$4; iz=1.5*\$3/\$4; printf("%i %i %i \n",ix,iy,iz)}')
         echo  "npts \${SIZE}          # num.grid points in xyz" > receptor.gpf
         grep -v npts  $WORKDIR/REC_${irec}/receptor.gpf | grep -v gridcenter | grep -v 'IATCEN' | grep -v 'Grid centered at ATOM' | grep -v 'IZTAM' >> receptor.gpf
         ed receptor.gpf <<EOX
/map receptor
i
gridcenter \${CEN_LIG}           # xyz-coordinates or auto
.
wq
EOX

         echo 'ROOT' >  ligand.pdbqt
         grep 'ATOM ' temp.pdbqt >> ligand.pdbqt
         echo 'ENDROOT' >> ligand.pdbqt
         echo 'TORSDOF 0' >> ligand.pdbqt
         rm -f temp.pdbqt

         sed -n '1,/move/p' $WORKDIR/receptor_ligand.dpf > receptor_ligand.dpf
         echo 'epdb' >> receptor_ligand.dpf

         echo "cd $WORKDIR/REC_${irec}/LIG_\${ilig}/OPT_MM/\${pose}_RESCORE; \
         $ADCK/autogrid4 -p receptor.gpf -l receptor.glg;    \
         $ADCK/autodock4 -p receptor_ligand.dpf -l receptor_ligand.dlg; \
         rm -f receptor.*.map" >> $WORKDIR/TASK_adck_${irec}.sh
         cd ../
      done
      cd ../../
   done
   cd ../
EOF
done

cd $WORKDIR
echo "Running ADCK preparatory tasks for rescoring in $PWD"
rm -f TASK.sh 
for file in $(ls TASK_prep_*)
do
    chmod 755 $file
    echo "$WORKDIR/${file}" >> TASK.sh
done 
cat TASK.sh  | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS

cd $WORKDIR
echo "Running AUTOGRID/AUTODOCK rescoring as parallel tasks in $PWD"
rm -f TASK.sh 
for file in $(ls TASK_adck_*)
do
    cat $file >> TASK.sh 
done 
cat TASK.sh  | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS
rm -f TASK_prep_* TASK_adck_*

fi   # END OF DO_RESCORE IF 

# Triple loop over receptor, ligand and pose files to get final scores
rm -f SCORING*.dat
if [ -e  PDB_SCORING ] ; then rm -r -f PDB_SCORING ; fi
mkdir PDB_SCORING

for ((irec=1;irec<=nrec;irec++))
do
    
rm -f TASK_extract_${irec}.sh
cat <<EOF > TASK_extract_${irec}.sh 
#!/bin/bash

   cd $WORKDIR/REC_${irec} 

   for ((ilig=1;ilig<=${nlig};ilig++))
   do
      cd LIG_\${ilig}/OPT_MM
      for ipose in \$(cat NPOSE.dat )
      do 
         pose=pose_\${ipose}
         cd \${pose}_RESCORE
         energy=\$(grep "Final Intermolecular Energy" receptor_ligand.dlg | awk '{print \$8}' )
         echo "complex_${irec}_\${ilig}_\${ipose} \$energy" >> $WORKDIR/SCORING_${irec}.dat 
         sed  '1,/Intermolecular Energy Analysis/d' receptor_ligand.dlg  | sed '/Estimated Free Energy of Binding/,\$d'  > $WORKDIR/PDB_SCORING/complex_${irec}_\${ilig}_\${ipose}.dlg 
         grep 'epdb: USER' receptor_ligand.dlg | sed 's/epdb: USER/COMMENT /' > $WORKDIR/PDB_SCORING/complex_${irec}_\${ilig}_\${ipose}.pdb
         cat receptor.pdb ligand.pdb >> $WORKDIR/PDB_SCORING/complex_${irec}_\${ilig}_\${ipose}.pdb
         cd ../
      done
      cd ../../
   done
   cd ../
EOF
done

cd $WORKDIR
echo "Running extraction tasks in $PWD"
rm -f TASK.sh 
for file in $(ls TASK_extract_*)
do
    chmod 755 $file
    echo "$WORKDIR/${file}" >> TASK.sh
done 
cat TASK.sh  | $PARHOME/bin/parallel --silent --no-notice  -t -j$NPROCS
rm -f TASK_extract_*

# Scoring file 
cat SCORING_*.dat > SCORING.dat
sed -i 's/+/ /g' SCORING.dat
sort -n -k 2 SCORING.dat > tmp; mv tmp SCORING.dat
rm -f SCORING_*.dat

cd $WORKDIR
rm -f mlog TASK.sh


echo 'If everything is OK, you should clean up intermediate files:'
echo '                    Run the clean_up_REC_LIG_data.sh script!'
rm -f clean_up_REC_LIG_data.sh
echo  "TOP_PERCEN=${TOP_PERCEN} " >> clean_up_REC_LIG_data.sh
echo  "WORKDIR=$WORKDIR" >> clean_up_REC_LIG_data.sh
echo  'cd $WORKDIR' >> clean_up_REC_LIG_data.sh
for ((irec=1;irec<=nrec;irec++))
do
    echo 'rm -r -f $WORKDIR/REC_'${irec} >> clean_up_REC_LIG_data.sh
done
echo  'NPOSES=$(cat SCORING.dat | wc -l)' >> clean_up_REC_LIG_data.sh
echo  'let "NCLEAN=$NPOSES * (100 - $TOP_PERCEN ) / 100" ' >> clean_up_REC_LIG_data.sh
echo  'let "NKEEP=$NPOSES - $NCLEAN" ' >> clean_up_REC_LIG_data.sh
echo  'if [ $NCLEAN -gt 0 ]  ' >> clean_up_REC_LIG_data.sh
echo  'then ' >> clean_up_REC_LIG_data.sh
echo  '  cd ${WORKDIR}/PDB_SCORING  ' >> clean_up_REC_LIG_data.sh
echo  '  mkdir ../PDB_SCORING_TOP_${TOP_PERCEN}' >> clean_up_REC_LIG_data.sh
echo  '  cp $(cat ../SCORING.dat | head -${NKEEP} | awk  ' \' '{printf("%s.* ",$1)}' \' ')  ../PDB_SCORING_TOP_${TOP_PERCEN} '  >> clean_up_REC_LIG_data.sh
echo  '  cd ${WORKDIR} ' >> clean_up_REC_LIG_data.sh
echo  '  rm -r -f PDB_SCORING' >> clean_up_REC_LIG_data.sh 
echo  '  head -${NKEEP} SCORING.dat  > tmp; mv tmp SCORING.dat '  >> clean_up_REC_LIG_data.sh
echo  'fi '  >> clean_up_REC_LIG_data.sh
chmod 755 clean_up_REC_LIG_data.sh

