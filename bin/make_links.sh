#!/bin/bash
if [ -z "$APTAMD" ]; then echo "APTAMD variable is not defined!" ; exit; fi
# Fix permissions
cd $APTAMD
chmod u=rwx,g=rx,o=rx $(find . -name "*.sh")
# Bin directory
cd $APTAMD/bin
# Do links
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
