#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 44

out="call"
$atlas --task call --window 4567 \
	   --bam simulate.bam --fasta simulate.fasta \
	   --fixedSeed 42 --out $out --logFile $out.out 2> $out.eout
