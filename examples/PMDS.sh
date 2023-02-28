#! /bin/bash

. $(dirname $0)/find_atlas

for strand in single double; do
	$atlas --task simulate --pmd "${strand}Strand:Exponential[30,0.1,0.1,0.05]:Exponential[40,0.2,0.3,0.07]"  --fixedSeed 0 --logFile simulate.out --out $strand
	$atlas --task PMDS --bam $strand.bam --fasta $strand.fasta --out $strand --logFile $strand.out
done
