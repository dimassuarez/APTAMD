#!/bin/bash

if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
if [ "$#" -eq 0 ]; then more $APTAMD/DOC/do_snapshots.txt ; exit; fi

source $APTAMD/ENV/geuo_env.sh

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

if [ -z "$MD_TRAJ" ]; then echo 'Usage: do_snapshots.sh "MOL1   MOL2 .."   [ MD | GAMD ]  '; exit ; fi
if [ -z "$MD_TYPE" ]; then echo 'Usage: do_snapshots.sh "MOL1   MOL2 .."   [ MD | GAMD ]  '; exit ; fi

if [ -z "$MAX_NA_GLOBAL" ]; then echo 'Considering all Na+ in TOPOLOGY'; ALL_NA="YES"; else ALL_NA="NO"; echo "Considering MAX_NA_GLOBAL=$MAX_NA_GLOBAL" ; fi
if [ -z "$MD_PROD" ]; then MD_PROD="5.PRODUCTION"; else echo "Assuming MD_PROD=$MD_PROD"; fi
if [ -z "$PEEL" ]; then echo 'Considering PEEL=14';  PEEL=16; else echo "Using PEEL=$PEEL"; fi
if [ -z "$SIEVE" ]; then echo 'Considering SIEVE=50';  SIEVE=50; else echo "Using SIEVE=$SIEVE"; fi
if [ -z "$DISLIM" ]; then echo 'Considering DISLIM=3.0';  DISLIM=3.0; else echo "Using DISLIM=$DISLIM"; fi
if [ -z "$SUFFIX_MDCRD" ]; then echo 'Considering SUFFIX_MDCRD=_solutewat.mdcrd';  SUFFIX_MDCRD="_solutewat.mdcrd"; else echo "Using SUFFIX_MDCRD=$SUFFIX_MDCRD"; fi
if [ -z "$PREFIX_MDCRD" ]; then echo 'Considering PREFIX_MDCRD=md_';  PREFIX_MDCRD="md_"; else echo "Using PREFIX_MDCRD=$PREFIX_MDCRD"; fi
if [ -z "$SKIP_PDB" ]; then echo 'Considering SKIP_PDB=NO';  SKIP_PDB="NO"; else echo "Using SKIP_PDB=$SKIP_PDB"; fi
if [ -z "$SNAPSHOTS_DIR" ]
then
   SNAPSHOTS_DIR="SNAPSHOTS"
fi

if [ -n "$PBS_ENVIRONMENT" ] ; then
  NPROCS=$(cat $PBS_NODEFILE | wc -l)
fi

if [ -z "$NPROCS" ]
then
   NPROCS=$(cat /proc/cpuinfo | grep -c processor)
   echo "Using $NPROCS available processors"
else
   echo "Using NPROCS=$NPROCS processors as predefined "
fi
export NPROCS
export OMP_NUM_THREADS=$NPROCS


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
     if [ -z $TOPOLOGY ]; then TOPOLOGY=${MOL}.top; fi
  elif [ ${SUFFIX_MDCRD} == "_solute.mdcrd" ]  && [ -z $TOPOLOGY ]
  then
     echo "Processing solute trj files so that PEEL=0 MAX_NA=0"
     PEEL=0
     MAX_NA=0
     MAX_NA_GLOBAL=0
     if [ -z $TOPOLOGY ]; then TOPOLOGY=${MOL}_solute.top; fi
  elif [ ${SUFFIX_MDCRD} == "_solutewat.mdcrd" ]   && [ -z $TOPOLOGY ]
  then
     if [ -z $TOPOLOGY ]; then TOPOLOGY=${MOL}_solutewat.top; fi
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
  NTOP_SODIUM=$(sed  '/FLAG CHARGE/,$d' $TOPOLOGY | tr " \t" "\n" | grep -c 'Na+')
  NTOP_WAT=$(sed  '1,/FLAG RESIDUE_LABEL/d' $TOPOLOGY | tr " \t" "\n" | grep -c 'WAT')
  echo "$TOPOLOGY  has  NTOP_SODIUM=${NTOP_SODIUM}"
  echo "$TOPOLOGY  has  NTOP_WAT=${NTOP_WAT}"
  if [ "$ALL_NA" == "YES" ] 
  then 
     MAX_NA=${NTOP_SODIUM}
     MAX_NA_GLOBAL=${NTOP_SODIUM} 
     echo "Considering MAX_NA=$MAX_NA listed in $TOPOLOGY "
  else
     if [ $NTOP_SODIUM -gt $MAX_NA_GLOBAL ]
     then 
         echo "Using MAX_NA=$MAX_NA_GLOBAL as defined by MAX_NA_GLOBAL instead of MAX_NA=$NTOP_SODIUM as listed in $TOPOLOGY "
         MAX_NA=$MAX_NA_GLOBAL
     else
         echo "MAX_NA_GLOBAL=$MAX_NA_GLOBAL, but only MAX_NA=$NTOP_SODIUM available as listed in $TOPOLOGY "
         MAX_NA=${NTOP_SODIUM}
     fi
  fi
  if [ $NTOP_WAT  -eq  0 ]
  then 
     echo "Since $TOPOLOGY has no water, then PEEL=0"
     PEEL="0"
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
trajout $SOLUTE_PDB  pdb pdbatom include_ep
go
EOF
    NRES=$(grep 'ATOM  ' $SOLUTE_PDB   | tail -1 | awk '{print $5}')
    NFRAG=$(grep -c 'TER' $SOLUTE_PDB)
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

  if [ ! -e 6.ANALYSIS ]
  then 
     echo "6.ANALYSIS directory not found in ${MOL}_${MD_TYPE}. Making it"
     mkdir 6.ANALYSIS 
  fi
  cd 6.ANALYSIS

  if [ ! -e ${SNAPSHOTS_DIR} ]
  then 
      echo "Creating ${SNAPSHOTS_DIR} directory in ${MOL}_${MD_TYPE}/6.ANALYSIS"
      mkdir ${SNAPSHOTS_DIR}
  else
      echo "Found ${SNAPSHOTS_DIR} directory in ${MOL}_MD/${MD_TYPE}"
  fi
  cd ${SNAPSHOTS_DIR}
  WORKDIR=$PWD 

  # Printing out options
  cat <<EOF > options.txt
MOL=$MOL     
TOPOLOGY=$TOPOLOGY
NRES=$NRES
NFRAG=$NFRAG

SNAPSHOTS_DIR=$SNAPSHOTS_DIR
MAX_NA_GLOBAL=$MAX_NA_GLOBAL
PEEL=$PEEL
SIEVE=$SIEVE
DISLIM=$DISLIM
SKIP_PDB=$SKIP_PDB

MD_TYPE=$MD_TYPE
MD_PROD=$MD_PROD
SUFFIX_MDCRD=$SUFFIX_MDCRD
PREFIX_MDCRD=$PREFIX_MDCRD

NPROCS=$NPROCS
EOF

  if [ ! -e processed_files.txt ] ; then touch processed_files.txt; fi 
  if [ ! -e LISTA ] ; then touch LISTA; fi 

#  Get SNAP PDBs (full atom)
  rm -f TASK.sh trajin.txt lista.txt current_processed_files.txt current_LISTA
  ifile=0
  jfile=0
  INIT=0
  for file in $(ls ../../${MD_PROD}/${PREFIX_MDCRD}???${SUFFIX_MDCRD} )
  do 
      let "ifile=$ifile+1"
      NSNAP=$(ncdump -h $file |grep frame | grep UNLI | awk '{print $6}' | sed 's/(//')
      ncheck=$(grep "$file" processed_files.txt | grep -c "$file  1 ${NSNAP} ${SIEVE} ${PEEL} ${MAX_NA}" )      
      if [ $ncheck -eq 0 ]
      then 
          let "jfile=$jfile+1"
          echo "$file  1 $NSNAP ${SIEVE}"  >> lista.txt 
          echo "trajin $file  1 $NSNAP ${SIEVE} "  >> trajin.txt 
      else
          echo "$file  1 $NSNAP ${SIEVE} already processed."
          let "INIT=$INIT + ${NSNAP}/${SIEVE}"
      fi
      echo "$file  1 $NSNAP ${SIEVE} ${PEEL} ${MAX_NA}"  >> current_processed_files.txt
  done
  nfile=$ifile
  echo "nfile=$nfile"
  if [ $nfile -eq 0 ]
  then
      echo "ls ../../${MD_PROD}/${PREFIX_MDCRD}???${SUFFIX_MDCRD}  FAILED!!"
      echo "Trying ls ../../${MD_PROD}/${PREFIX_MDCRD}??${SUFFIX_MDCRD}"
      ifile=0
      jfile=0
      INIT=0
      for file in $(ls ../../${MD_PROD}/${PREFIX_MDCRD}??${SUFFIX_MDCRD} )
      do 
          let "ifile=$ifile+1"
          NSNAP=$(ncdump -h $file |grep frame | grep UNLI | awk '{print $6}' | sed 's/(//')
          ncheck=$(grep "$file" processed_files.txt | grep -c "$file  1 ${NSNAP} ${SIEVE} ${PEEL} ${MAX_NA}" )      
          if [ $ncheck -eq 0 ]
          then 
              let "jfile=$jfile+1"
              echo "$file  1 $NSNAP ${SIEVE}"  >> lista.txt 
              echo "trajin $file  1 $NSNAP ${SIEVE} "  >> trajin.txt 
          else
              echo "$file  1 $NSNAP ${SIEVE} already processed."
              let "INIT=$INIT + ${NSNAP}/${SIEVE}"
          fi
          echo "$file  1 $NSNAP ${SIEVE} ${PEEL} ${MAX_NA}"  >> current_processed_files.txt
      done
      nfile=$ifile
      if [ $nfile -eq 0 ]
      then
          echo "ls ../../${MD_PROD}/${PREFIX_MDCRD}??${SUFFIX_MDCRD}  FAILED!!"
          echo "No way"
          exit
      fi
  fi
  if [ $jfile -eq 0 ]  # No new file to be processed
  then 
     SKIP_PDB="YES"
  else
     echo "Processing the following files:"
     cat lista.txt
  fi
  if [ $jfile == $ifile ]    # If all mdcrds are again processed then remove previous CLOSEST_*dat files
  then
     rm -f CLOSEST_*dat
  fi

  if [ "$SKIP_PDB" != "YES" ]
  then

  echo "$APTAMD/SCRIPTS/peel_cpptraj.sh ../../$TOPOLOGY   lista.txt  $PEEL ":1-${NRES},Na+"  ${MOL}  ${INIT}"
  $APTAMD/STRUCT/peel_cpptraj.sh ../../$TOPOLOGY lista.txt $PEEL ":1-${NRES},Na+"  ${MOL}   ${INIT}  
  if [ $jfile != $ifile ]; then ls ${MOL}*.pdb.gz | sed 's/.gz//' > current_LISTA; fi 
  ls ${MOL}*.pdb >  LISTA.tmp
  cat LISTA.tmp >> current_LISTA

  # Gzipping PDB files 
  rm -f TASK.sh
  npdb=$(cat LISTA.tmp | wc -l)
  let "n=$npdb/$NPROCS+1"
  split -l${n}  -d LISTA.tmp  tmp_lista_
  for file in $(ls tmp_lista_*)  
  do
     echo 'gzip -f $(cat '${file}')' >> TASK.sh
  done
  cat TASK.sh  | $PARHOME/bin/parallel  -t -j$NPROCS >/dev/null 2>&1
  rm -f tmp_lista_* LISTA.tmp

  else 
     if [ -e LISTA ] ; then cp LISTA current_LISTA; fi
  fi


  if [ ${MAX_NA_GLOBAL} -gt 0 ]
  then

# Get Na+ contacts data for each fragment 

  for ((IFRAG=0;IFRAG<=NFRAG;IFRAG++))
  do


  let "JFRAG=$IFRAG+1"

  rm -f closestout_${JFRAG} closestmols_${JFRAG}.dat closest_${JFRAG}.in

  if [ "$IFRAG" -lt "$NFRAG" ]
  then

  KRES=${IRES["$IFRAG"]} 
  LRES=${JRES["$IFRAG"]} 
  echo "Getting Na+ contacts for fragment $JFRAG = $KRES : $LRES"

  else

  KRES="1"
  LRES="$NRES"
  JFRAG="full"
  echo "Getting Na+ contacts for all solute atoms"

  fi

  if [ $jfile -gt 0 ]
  then 
     cat trajin.txt > closest_${JFRAG}.in 
     cat <<EOF >>  closest_${JFRAG}.in 
closest ${MAX_NA}  :${KRES}-${LRES} solventmask :Na+ closestout closestmols_${JFRAG}.dat
hbond HB-Ion avgout contacts_${JFRAG}.dat solventacceptor :Na+ solventdonor :Na+
go

EOF
     $AMBERHOME/bin/cpptraj.OMP ../../$TOPOLOGY < closest_${JFRAG}.in > closest_${JFRAG}.out 
     if [ -e CLOSEST_${JFRAG}.dat ]
     then 
          sed '1,1d' closestmols_${JFRAG}.dat >> CLOSEST_${JFRAG}.dat
          rm -f closestmols_${JFRAG}.dat 
     else
          mv closestmols_${JFRAG}.dat CLOSEST_${JFRAG}.dat
     fi
  fi


cat <<EOF > na_${JFRAG}.m 
fid_1=fopen('NA_SORTED_${JFRAG}.INFO','w');
fid_2=fopen('NA_DIST_${JFRAG}.INFO','w');
A=load('CLOSEST_${JFRAG}.dat');
nion=${MAX_NA};
DISLIM=${DISLIM};
nsnap=length(A(:,1))/nion;
DIS=A(:,3);
INA=A(:,4);
INA_MAT=reshape(INA,nion,nsnap);
DIS_MAT=reshape(DIS,nion,nsnap);
INA_MAT=INA_MAT';
DIS_MAT=DIS_MAT';
fmt='';
for i=[1:nion]
    fmt=[fmt,'%i,'];
end
    fmt=[fmt,'\n'];
for i=[1:nsnap]
    fprintf(fid_1,fmt,INA_MAT(i,:))
end

for k=[1:1:nion]

for i=[1:nsnap]
   DIS_AVG(i)=mean(DIS_MAT(i,1:k));
   DIS_MEDIAN(i)=median(DIS_MAT(i,1:k));
end

% Accumulated PDF 
wx=0.025;  % Bin size
DIS_FULL=reshape(DIS_MAT(:,1:k),nsnap*k,1);
xmin=min(DIS_FULL)-wx;
xmax=max(DIS_FULL)+wx;
NBINS= round ( (  xmax - xmin ) / wx ) + 1  ;
[ P, X ] = hist( DIS_FULL, NBINS , 1.0) ;
prob=0.0;
for i=[1:length(P)]
    if ( X(i) < DISLIM ) 
        prob=prob+P(i);
    end
end
fprintf(fid_2,' =============================================================== \n ')
fprintf(fid_2,'Na···Solute Distance Distribution for #ion= %i  \n ',k)
fprintf(fid_2,'=============================================================== \n ')
fprintf(fid_2,'Prob  ( DIS < %5.2f ) for %i Na+ ions = %6.4f \n ',DISLIM,k,prob)
fprintf(fid_2,'   PMF  Na+···X  \n ')
for i=[1:length(P)]
   if ( P(i) > 0.001 )
      fprintf(fid_2,' %6.4f  %6.4f \n ',P(i),X(i))
   end
end
%Discretized Entropy
R=8.3145/4.184;
S=0.0;
for i=[1:NBINS]
  if P(i) > 0.0
     S=S-R*P(i)*log(P(i));
  end
end
S=k*S; % To make S extensive 
fprintf(fid_2,'Discretized-full Shannon-Entropy  #ion= %i  #bins= %i  rmin= %10.3f  rmax= %10.3f   S=  %10.5f \n', ...
        k,NBINS,xmin,xmax,S)

NBINS=round((max(DIS_AVG)-min(DIS_AVG))/wx)+1;
xmin=min(DIS_AVG)-wx ;
xmax=max(DIS_AVG)+wx ;
[P,X]=hist(DIS_AVG,NBINS,1.0);
prob=0.0;
for i=[1:length(P)]
    if ( X(i) < DISLIM ) 
        prob=prob+P(i);
    end
end
fprintf(fid_2,' ========================================================================== \n ')
fprintf(fid_2,'Na···Solute AVERAGE Distance Distribution for #ion= %i \n ',k)
fprintf(fid_2,'========================================================================== \n ')
fprintf(fid_2,'Prob  ( AVG-DIS < %5.2f ) for %i Na+ ions = %6.4f \n ',DISLIM,k,prob)
fprintf(fid_2,'   PMF  AVG_Na+···X  \n ')
for i=[1:length(P)]
   if (P(i) > 0.001 )
      fprintf(fid_2,' %6.4f  %6.4f \n ',P(i),X(i))
   end
end
S=0.0d0;
for i=[1:NBINS]
  if P(i) > 0.0
     S=S-R*P(i)*log(P(i));
  end
end
S=k*S;
fprintf(fid_2,'Discretized-AVG  Shannon-Entropy  #ion= %i  #bins= %i  rmin= %10.3f  rmax= %10.3f   S=  %10.5f \n \n', ...
        k,NBINS,xmin,xmax,S)

%Representacion grafica 
clf();
h=figure(1);
subplot(2,1,1)
hist(DIS_AVG,50,1, "facecolor", "r", "edgecolor", "black")
grid on;
W = 8; H = 6;
set(h,'PaperUnits','inches')
set(h,'PaperOrientation','portrait');
set(h,'PaperSize',[H,W])
set(h,'PaperPosition',[0,0,W,H])
xmin= fix(min(DIS_AVG)*0.9);
xmax= round(max(DIS_AVG)*1.1);
wx=(xmax-xmin)/5;
xlim([ xmin  xmax ])
ylim([ 0 0.20])
xlabel(' r ') 
ylabel('frequency')
title([' Frag${JFRAG}---Na+ Dist Closest ',num2str(k),' ions'])
text(xmin+wx/2,0.18,'Average','Fontsize',18)
text(xmin+2*wx,0.18,['P(<',num2str(DISLIM,'%4.2f'),')=',num2str(prob,'%4.2f')],'Fontsize',14)
set(gca,'xtick',[xmin:wx:xmax])
set(gca,'ytick',[0:0.05:0.20])
set(gca,'ticklength',[0.01 0.01])
set(gca,'Fontsize',18)
set(gca,'Fontname','Times')

subplot(2,1,2)
hist(DIS_MEDIAN,50,1, "facecolor", "r", "edgecolor", "black")
grid on;
xmin= fix(min(DIS_MEDIAN)*0.9);
xmax= round(max(DIS_MEDIAN*1.1));
wx=(xmax-xmin)/5;
xlim([ xmin  xmax ])
ylim([ 0 0.20])
xlabel(' r ') 
ylabel('frequency')
text(xmin+wx/2,0.18,'Median','Fontsize',18)
set(gca,'xtick',[xmin:wx:xmax])
set(gca,'ytick',[0:0.05:0.20])
set(gca,'ticklength',[0.01 0.01])
set(gca,'Fontsize',18)
set(gca,'Fontname','Times')

print(h,['hist_${JFRAG}_',num2str(k),'.png'],'-dpng','-color')

end

fclose(fid_1);
fclose(fid_2);

EOF
  $OCTAVE --no-gui -W -q  na_${JFRAG}.m  > mlog
  paste current_LISTA NA_SORTED_${JFRAG}.INFO | sed 's/\t/,/g' >  TMP ; mv -f TMP NA_SORTED_${JFRAG}.INFO

  if [ $JFRAG -eq 1 ] && [ $NFRAG -eq 1 ]; then break; fi

  done 

  fi

  unset TOPOLOGY 

  mv -f  current_LISTA  LISTA
  mv -f  current_processed_files.txt  processed_files.txt

  cd $WORKDIR_TRJ 

done 

