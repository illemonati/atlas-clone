#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 222

out="sexEstimation"
$atlas --task sexEstimation --wholeGenome --window 4567 \
	   --bam simulate.bam --fasta simulate.fasta \
	   --fixedSeed 221 --out $out --logFile $out.out 2> $out.eout
