#!/bin/bash

# Auxiliary script for do_rmsd_rgyr_inf.sh

source ./environment.sh
  
MDCRD=$1
ID=$2

echo "MDCRD=$1"
echo "ID=$2"

TT=$(date +%N)
TMPDIR=${SCRATCH}/TMPDIR_${TT}

 # All work is done in a temporal directory
mkdir $TMPDIR
cd $TMPDIR

cp ${WORKDIR}/${MOL}_ref.dat .
NSNAP=$(ncdump -h $WORKDIR/${MDCRD}  |grep frame | grep UNLI | awk '{print $6}' | sed 's/(//')

echo "trajin $WORKDIR/${MDCRD} 1 ${NSNAP} ${SIEVE} " > temp.in
echo "reference ${REFPDB}" >> temp.in
echo "rms reference  ${MOL_MASK}  out ${ID}.rmsd  " >> temp.in
echo "radgyr :1-${NRES} out ${ID}.rgyr " >> temp.in
if [ ${DO_SURF} == "YES" ]
then
echo "molsurf :1-${NRES} out ${ID}.surf probe 1.4 radii vdW " >> temp.in
fi
if [ ${DO_INF} == "YES" ]
then
echo "strip !:${MOL_MASK}" >> temp.in
echo "trajout ${ID}.pdb  pdb" >> temp.in
fi

$AMBERHOME/bin/cpptraj $WORKDIR/../../$TOPOLOGY  < temp.in

cp  ${ID}.rmsd ${ID}.rgyr  $WORKDIR/
if [ ${DO_SURF} == "YES" ]; then cp  ${ID}.surf $WORKDIR/ ; fi 

# Correlation plots of each segment (uncomment if required)
# 
# cp  $MFILE_RMSD_RGYR_SURF ${ID}.m
# sed -i "s/DUMMY_RMS/${ID}.rmsd/" ${ID}.m
# sed -i "s/DUMMY_RGYR/${ID}.rgyr/" ${ID}.m
# sed -i "s/DUMMY_SURF/${ID}.surf/" ${ID}.m
# sed -i "s/DUMMY_PNG/${ID}/" ${ID}.m
# mv ${ID}.m ${ID}_rmsd_rgyr_surf.m
# $OCTAVE --no-gui  -q ${ID}_rmsd_rgyr_surf.m 
# cp ${ID}_rmsd_rgyr_surf.png ${ID}_rmsd_rgyr_surf.m $WORKDIR/
#

if [ ${DO_INF} == "YES" ]
then

#Getting NP and BP data. Working with temporal space in RAM memory
 rm -r -f /dev/shm/TMPSHM_${ID}_${TT}
 mkdir /dev/shm/TMPSHM_${ID}_${TT}
 cd /dev/shm/TMPSHM_${ID}_${TT}

 csplit -n 5 -s -f model_  $TMPDIR/${ID}.pdb '/MODEL /' '{*}'
 rm -f model_00000 

 imodel=0
 for model in $(ls model_*)
 do

  let "imodel=$imodel+1"
  $X3DNAHOME/bin/x3dna-dssr --nmr  --no-pair --auxfile=no -i=${model} -o=tmp.out_x3dna  >& tmp.log 

  icheck=$(grep -c 'Summary of structural features' tmp.out_x3dna)
  if [ $icheck -eq 0 ]
  then

  echo "-1"  >> ${ID}_bp_n.dat
  echo "-1"  >> ${ID}_np_n.dat

  else

  grep -n 'List of ' tmp.out_x3dna | grep 'base pairs' | sed 's/:/  /' |awk '{L=$1+$4+1;printf("%i,%ip\n",$1,L)}' > tmp_bp.sed
  sed -n -f tmp_bp.sed tmp.out_x3dna > tmp_bp.dat
  sed -i 's/:D/   D/g' tmp_bp.dat
  grep List tmp_bp.dat  | awk '{print $3}' > tmp_bp_n.dat
  grep -v List tmp_bp.dat | grep -v DSSR   | awk '{print $3, $5}' | sed 's/3\///' | sed 's/5\///' | sed 's/D[A-T]/   /g' > tmp_bp_ij.dat

  #Getting data for non-pairing interactions, either direct base contacts or base stacking 
  grep -n 'List of ' tmp.out_x3dna | grep 'non-pairing interactions' | sed 's/:/  /' |awk '{l=$1+$4;printf("%i,%ip\n",$1,l)}' > tmp_np.sed
  sed -n -f tmp_np.sed tmp.out_x3dna > tmp_np.dat
  sed -i 's/:D/   D/g' tmp_np.dat
  grep 'List of' tmp_np.dat  | awk '{print $3}' > tmp_np_n.dat

  sed -i 's/interBase/ 0 interBase/g' tmp_np.dat
  sed -i 's/stacking:/ 1 stacking/g' tmp_np.dat
  grep -v List tmp_np.dat | awk '{print $3, $5, $6}' | sed 's/3\///' | sed 's/5\///' | sed 's/D[A-T]/   /g' > tmp_np_ij.dat

  cat tmp_bp_n.dat  >> ${ID}_bp_n.dat
  cat tmp_bp_ij.dat >> ${ID}_bp_ij.dat
  cat tmp_np_n.dat  >> ${ID}_np_n.dat
  cat tmp_np_ij.dat >> ${ID}_np_ij.dat

  fi
  rm -f tmp.*

 done

 mv ${ID}_bp_ij.dat ${ID}_bp_n.dat ${ID}_np_ij.dat ${ID}_np_n.dat $TMPDIR/

 cd $TMPDIR/

 rm -r -f /dev/shm/TMPSHM_${ID}_${TT}

#Running octave
cp $MFILE_RMSD_INF  ${ID}.m

sed -i "s/DUMMY_NRES/${NRES}/" ${ID}.m
sed -i "s/DUMMY_BP_N/${ID}_bp_n.dat/" ${ID}.m
sed -i "s/DUMMY_NP_N/${ID}_np_n.dat/" ${ID}.m
sed -i "s/DUMMY_BP_IJ/${ID}_bp_ij.dat/" ${ID}.m
sed -i "s/DUMMY_NP_IJ/${ID}_np_ij.dat/" ${ID}.m
sed -i "s/DUMMY_REF/${MOL}_ref.dat/" ${ID}.m
sed -i "s/DUMMY_FOUT/${ID}.inf/" ${ID}.m
sed -i "s/DUMMY_RMS/${ID}.rmsd/" ${ID}.m
sed -i "s/DUMMY_PNG/${ID}/" ${ID}.m
sed -i "s/DUMMY_DOPLOT/0/" ${ID}.m    # Replace DUMMY_DOPLOT with 1 for plotting

$OCTAVE -q ${ID}.m

# Copy back all data
cp ${ID}_np_ij.dat ${ID}_np_n.dat ${ID}_bp_ij.dat ${ID}_bp_n.dat ${ID}.inf ${ID}.m  $WORKDIR/
if [ -e ${ID}_1.png ]; then cp ${ID}*.png  $WORKDIR/; fi

fi

cd $WORKDIR

rm -r -f $TMPDIR 

