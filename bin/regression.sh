#!/bin/bash
if [ ! -d testo ]; then
    mkdir testo
fi

rm -rf testo/*

O='\e[0m'
G='\e[1;32m'
R='\e[1;31m'
A='\e[1;30m' # grey

echo "start annotations..."
time -p julia --threads=8 "$@" chloe.jl -l info annotate -o testo testfa/*.fa
for f in $(ls testo)
do 
    echo -e "diffing: $f ${A}(> means line from original)${O}"
    diff testo/$f testfa/$f
    if [ $? -eq 0 ]; then
        echo -e "$G******** test OK ***********$O"
    else
        echo -e "$R******** test FAILED *******$O"
    fi
done
# rm -rf testo
