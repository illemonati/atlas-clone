#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 22

out="allelicDepth"
$atlas --task allelicDepth --readUpToDepth 97 \
	   --bam simulate.bam --window 4567 \
	   --fixedSeed 25 --out $out --logFile $out.out 2> $out.eout
