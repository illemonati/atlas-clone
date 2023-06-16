#! /bin/bash

. $(dirname $0)/find_atlas

for strand in single double; do
	. $(dirname $0)/simulate --pmd "${strand}Strand:Exponential[30,0.1,0.1,0.05]:Exponential[40,0.2,0.3,0.07]" --out $strand --depth 100

	$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "Exponential"  --fixedSeed 0 --out ${strand}_exponential --logFile ${strand}_exponential.out
	$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "Empiric"  --fixedSeed 0 --out ${strand}_empiric --logFile ${strand}_empiric.out
	#$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "${strand}Strand:[80]:Empiric"  --fixedSeed 0 --out ${strand}_empiric_per --logFile ${strand}_empiric_per.out

	# Parsing output file
	$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "${strand}_exponential_PMD.txt" --fixedSeed 0 --out ${strand}_exponential_reestimate --logFile ${strand}_exponential_reestimate.out
	$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "${strand}_empiric_PMD.txt" ---fixedSeed 0 --out ${strand}_empiric_reestimte --logFile ${strand}_empiric_reestimte.out
	#$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "${strand}_empiric_per_PMD.txt" ---fixedSeed 0 --out ${strand}_empiric_per_reestimte --logFile ${strand}_empiric_per_reestimte.out
done
