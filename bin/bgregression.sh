#!/bin/bash
if [ ! -d testo ]; then
    mkdir testo
fi
rm -rf testo/*
echo "ensure: make run-chloe-broker"
for f in $(ls testfa/*.fa)
do
    echo "annotating $f"
    python bin/chloe.py annotate -o testo $f 
    bn=$(basename $f)
    o="${bn%.*}.sff"
    echo "diffing testfa/$o testo/$o"
    diff testfa/$o testo/$o
    if [ $? -eq 0 ]; then
        echo -e "\e[32m******** test OK ***********\e[0m"
    else
        echo -e "\e[31m******** test FAILED *******\e[0m"
    fi
done
# python bin/chloe.py terminate -a $A
rm -rf testo
