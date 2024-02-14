#!/bin/bash
if [ ! -d testo ]; then
    mkdir testo
fi
rm -rf testo/*
O='\e[0m'
G='\e[1;32m'
R='\e[1;31m'
echo -e "ensure: ${G}make run-chloe-broker$O"
fafiles=$(ls testfa/*.fa)
python bin/chloe.py annotate -o testo --workers=4 $fafiles
for f in $fafiles
do  
    bn=$(basename $f)
    o="${bn%.*}.sff"
    echo "diffing testfa/$o testo/$o"
    diff testfa/$o testo/$o
    if [ $? -eq 0 ]; then
        echo -e "$G******** test OK ***********$O"
    else
        echo -e "$R******** test FAILED *******$O"
    fi
done
# python bin/chloe.py terminate -a $A
# rm -rf testo
