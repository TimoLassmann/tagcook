#!/bin/bash

#  make_benchmark_datasets.sh
#  
#
#  Created by timo on 16/06/2011.
#  Copyright 2011 RIKEN Yokohama Institute, Genome Exploration Research Group. All rights reserved.

NUMBARCODE=

RUNNAME=


function usage()
{
cat <<EOF
usage: $0 -b <numbarcode> -m <run_name>
EOF
exit 1;
}

while getopts b:n: opt
do
case ${opt} in
b) NUMBARCODE=${OPTARG};;
n) RUNNAME=${OPTARG};;
*) usage;;
esac
done

if [ "${NUMBARCODE}" = "" ]; then usage; fi
if [ "${RUNNAME}" = "" ]; then usage; fi

foo=$RUNNAME
tmpfile=${foo##*/}
base=${tmpfile%%.*}
echo $base

array=(0.01 0.02 0.03 0.04 0.05 0.06 0.07)
len=${#array[*]}

for (( i=4; i <= 12; i+=2 )); do
j=0

while [ $j -lt $len ]
do
tagdust -simulation $i  -numbarcode $NUMBARCODE  -e ${array[$j]}  -o  $base.$NUMBARCODE.$i.${array[$j]};
let j++
done
done