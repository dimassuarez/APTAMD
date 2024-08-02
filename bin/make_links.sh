#/bin/bash
# Doing Links (if not done)
DO_LIST=$(find ../ -name "do_*.sh")
for file in ${DO_LIST}
do
    lkfile=$(basename $file)
    lkfile=${lkfile/.sh/} 
    echo $file $lkfile
    if [ ! -e $lkfile ]
    then
         ln -s $file $lkfile
    fi
done
