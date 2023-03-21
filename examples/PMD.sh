#! /bin/bash

. $(dirname $0)/find_atlas

for strand in single double; do
	$atlas --task simulate --pmd "${strand}Strand:Exponential[30,0.1,0.1,0.05]:Exponential[40,0.2,0.3,0.07]"  --fixedSeed 0 --logFile simulate.out --out $strand

	$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "${strand}Strand:Exponential:Exponential"  --fixedSeed 0 --out ${strand}_exponential --logFile ${strand}_exponential.out
	$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "${strand}Strand:Empiric:Empiric"  --fixedSeed 0 --out ${strand}_empiric --logFile ${strand}_empiric.out

# Parsing output file
	$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "${strand}_exponential_PMD.txt" --fixedSeed 0 --out ${strand}_exponential_reestimate --logFile ${strand}_exponential_reestimate.out
	$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "${strand}_empiric_PMD.txt" ---fixedSeed 0 --out ${strand}_empiric_reestimte --logFile ${strand}_empiric_reestimte.out
done
