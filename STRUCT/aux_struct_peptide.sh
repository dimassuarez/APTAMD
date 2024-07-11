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
echo "radgyr  :1-${NRES} out ${ID}.rgyr " >> temp.in 
if [ ${DO_SURF} == "YES"  ] 
then 
   echo "molsurf :1-${NRES} out ${ID}.surf probe 1.4 radii vdW " >> temp.in 
fi

$AMBERHOME/bin/cpptraj $WORKDIR/../../$TOPOLOGY  < temp.in 

cp  ${ID}.rmsd ${ID}.rgyr $WORKDIR/
if [ ${DO_SURF} == "YES"  ] 
then 
   cp  ${ID}.surf $WORKDIR/
fi 

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

  rm -f hb_n.dat hb.dat hb_uu.dat
  $AMBERHOME/bin/cpptraj ${WORKDIR}/../../$TOPOLOGY<<EOF
trajin  $REFPDB
trajin  $WORKDIR/${MDCRD} 1 ${NSNAP} ${SIEVE} 
hbond   hbset  out hb_n.dat :1-${NRES}  angle ${ANGLE} dist ${DIST}  printatomnum  avgout hb.dat series uuseries hb_uu.dat 
go
EOF

mv hb_uu.dat ${ID}_hb.dat

#Running octave
cp $MFILE_RMSD_INF  ${ID}.m

sed -i "s/DUMMY_HB/${ID}_hb.dat/" ${ID}.m
sed -i "s/DUMMY_FOUT/${ID}.inf/" ${ID}.m
sed -i "s/DUMMY_RMS/${ID}.rmsd/" ${ID}.m
sed -i "s/DUMMY_DOPLOT/0/" ${ID}.m       # Replace DUMMY_DOPLOT with 1 for plotting
sed -i "s/DUMMY_PNG/${ID}/" ${ID}.m

$OCTAVE --no-gui -q ${ID}.m

# Copy data
cp  ${ID}_hb.dat ${ID}.inf $WORKDIR/

fi

cd $WORKDIR

rm -r -f $TMPDIR 

