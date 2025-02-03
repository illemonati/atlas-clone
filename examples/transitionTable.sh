#! /bin/bash

. $(dirname $0)/find_atlas

. $(dirname $0)/simulate --fixedSeed 321 --seqType paired

out="transitionTable"
$atlas --task transitionTable --bam simulate.bam --fasta simulate.fasta \
	   --fixedSeed 123 --out $out --logFile $out.out 2> $out.eout
