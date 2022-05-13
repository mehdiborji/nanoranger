#!/bin/bash

echo 'cores =' $1
echo 'ref =' $2
echo 'infile=' $3
echo 'outdir =' $4
echo 'sample =' $5

minimap2 -aY --eqx -x map-ont -t $1  --secondary=no --sam-hit-only $2 $3 > $4/$5_trns.sam
