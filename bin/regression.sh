#!/bin/bash
if [ ! -d testo ]; then
    mkdir testo
fi
echo "start annotations..."
JULIA_NUM_THREADS=8 time julia src/chloe.jl -l info annotate -o testo testfa/*.fa
for f in $(ls testo)
do 
    echo "diffing $f"
    diff testo/$f testfa/$f
    if [ $? -eq 0 ]; then
        echo -e "\e[32m******** test OK ***********\e[0m"
    else
        echo -e "\e[31m******** test FAILED *******\e[0m"
    fi
done
rm -rf testo
