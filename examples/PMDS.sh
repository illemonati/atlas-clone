#! /bin/bash

. $(dirname $0)/find_atlas

for strand in single double; do
	. $(dirname $0)/simulate --pmd "${strand}Strand:Exponential[30,0.1,0.1,0.05]:Exponential[40,0.2,0.3,0.07]" --out $strand
	$atlas --task PMDS --bam $strand.bam --fasta $strand.fasta --fixedSeed 0 --out $strand --logFile $strand.out
done
