#! /bin/bash

. $(dirname $0)/find_atlas

bname=$(basename $0)
name=${bname%.sh}
out=$name/$name

$atlas --task recal --bam BAM/BAM.bam --fasta BAM/BAM.fasta --recal "intercept[0.1];quality:polynomial[0.8,0.2,0.01]" --chr chr2 --minDeltaLL 1 --fixedSeed 0 --out $out --logFile $out.out 2> $out.err > /dev/null
