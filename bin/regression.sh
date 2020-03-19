#!/bin/sh
mkdir testo
julia src/chloe.jl -o testo testfa/*.fa --level=info
for f in $(ls testo)
do 
    echo "diffing $f"
    diff testo/$f testfa/$f
done
rm -rf testo
