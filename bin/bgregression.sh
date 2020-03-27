#!/bin/sh
if [ ! -d testo ]; then
    mkdir testo
fi
A=tcp://127.0.0.1:9999
JULIA_NUM_THREADS=4 julia src/chloe_svr.jl -a $A --level=warn &
for f in $(ls testfa/*.fa)
do
    echo "annotating $f"
    python bin/chloe.py annotate -a $A -o testo $f 
    bn=$(basename $f)
    o="${bn%.*}.sff"
    echo "diffing testfa/$o testo/$o"
    diff testfa/$o testo/$o
done
python bin/chloe.py terminate -a $A
rm -rf testo
