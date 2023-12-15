#!/bin/bash
# select a random set of $TOTAL $1/*.fa files
# and compare the current sff output with the original
TOTAL=3
if [ $# -eq 0 ]; then
    echo "expecting DIR [num] [julia args]" 2>&1
    exit 1
fi
DIR=$1
shift
if [ ! -d $DIR ]; then
    echo "expecting directory" 2>&1
    exit 1
fi

if [ $# -gt 0 ]; then
    TOTAL=$1
    shift
fi

if [ ! -d testo ]; then
    mkdir testo
fi
rm -rf testo/*

O='\e[0m'
G='\e[1;32m'
R='\e[1;31m'
C='\e[1;36m' # bold cyan
A='\e[1;30m' # gray


# fa files as array
fafiles=($(ls $DIR/*.fa))
# number
n=${#fafiles[@]}
if (($n == 0)); then
    echo -e "$R no files!$O"
    exit 1
fi
# find random indexes
index=($(shuf -i 0-$((n - 1)) -n $TOTAL))
declare -a todo
for idx in ${index[@]}
do
    f=${fafiles[$idx]}
    todo=("${todo[@]}" $f)
done
echo -e "index: ${index[@]}: $TOTAL/${C}$n${O}"
echo "start annotations with: ${todo[@]}"

time -p julia --threads=8 --project=. "$@" chloe.jl -l warn annotate -o testo --numgsrefs 16 "${todo[@]}"

for idx in ${index[@]}
do
    f=${fafiles[$idx]}
    bn=$(basename $f)
    o="${bn%.*}.sff"
    echo -e "diffing[$idx]: $o ${A}(> means line from input)${O}"
    diff testo/$o $DIR/$o
    if [ $? -eq 0 ]; then
        echo -e "$G******** test OK ***********$O"
    else
        echo -e "$R******** test FAILED *******$O"
    fi
done
