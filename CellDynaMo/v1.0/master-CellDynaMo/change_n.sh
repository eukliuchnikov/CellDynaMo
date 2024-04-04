#!/bin/bash

file1="../config/conf.conf"
file2='src/Initializiation/parameters.h'

N=$(grep 'kin_count' $DPATH $file1 | awk '{print $2}')
#echo $N
n=$(grep 'KIN_COUNT' $DPATH $file2 | awk '{print $3}')
lineID=$(grep -in "KIN_COUNT" $file2 | awk -F: '{print $1}')

sed -i -e ''"$lineID"'s/'"$n"'/'"$N"'/g' $DPATH $file2

