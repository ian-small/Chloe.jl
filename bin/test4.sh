#!/bin/bash
if [ ! -d testo ]; then
    mkdir testo
fi
tfiles="testfa/NC_020019.1.fa testfa/NC_020318.1.fa testfa/NC_020152.1.fa testfa/NC_020320.1.fa"
exec python bin/chloe.py annotate -o testo --parallel  $tfiles
