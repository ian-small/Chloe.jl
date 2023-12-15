#!/bin/bash
if [ ! -d testo ]; then
    mkdir testo
fi
rm -rf testo/*

O='\e[0m'
G='\e[1;32m'
R='\e[1;31m'
C='\e[1;36m' # bold cyan
A='\e[1;30m' # grey

echo -e "testing ${C}${#}${O} files"
time -p julia --threads=8 --project=. chloe.jl -l warn annotate -o testo --numgsrefs 16 "$@"

for f in "$@"
do
    bn=$(basename $f)
    d=$(dirname $f)
    o="${bn%.*}.sff"
    echo -e "diffing: $o ${A}(> means line from input)${O}"
    diff testo/$o $d/$o
    if [ $? -eq 0 ]; then
        echo -e "$G******** test OK ***********$O"
    else
        echo -e "$R******** test FAILED *******$O"
    fi
done
