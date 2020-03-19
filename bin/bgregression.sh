#!/bin/sh
mkdir testo
JULIA_NUM_THREADS=4 julia src/chloe_svr.jl --level=warn &
for f in $(ls testfa/*.fa)
do
    echo "annotating $f"
    python chloe.py annotate $f -o testo
    bn=$(basename $f)
    o="${bn%.*}.sff"
    echo "diffing testfa/$o testo/$o"
    diff testfa/$o testo/$o
done
python chloe.py terminate
rm -rf testo
