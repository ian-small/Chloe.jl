#!/bin/bash
if [ ! -d testo ]; then
    mkdir testo
fi

rm -rf testo/*

O='\e[0m'
G='\e[1;32m'
R='\e[1;31m'

echo "start annotations..."
time -p julia --project=. --threads=8 "$@" chloe.jl -l info annotate -o testo --numgsrefs 16 --numchloerefs 0  testfa/*.fa
for f in $(ls testo)
do 
    echo "diffing $f"
    diff testo/$f testfa/$f
    if [ $? -eq 0 ]; then
        echo -e "$G******** test OK ***********$O"
    else
        echo -e "$R******** test FAILED *******$O"
    fi
done
# rm -rf testo
