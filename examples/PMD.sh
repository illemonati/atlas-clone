#! /bin/bash

pmd5="0.1*exp(-0.1*p)+0.05"
pmd3="0.2*exp(-0.3*p)+0.07"

. $(dirname $0)/find_atlas

for strand in single double; do
	pmd="CT5:$pmd5;CT3:$pmd3"
	[ $strand == "double" ] && pmd="CT5:$pmd5;GA3:$pmd3"

	. $(dirname $0)/simulate --pmd $pmd --out $strand --depth 100 --logFile simulate_$strand.out

	$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "Exponential"  --fixedSeed 0 --out ${strand}_exponential --logFile ${strand}_exponential.out
	$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "Empiric"  --fixedSeed 0 --out ${strand}_empiric --logFile ${strand}_empiric.out
	$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "Exponential[]"  --fixedSeed 0 --out ${strand}_empiric_per --logFile ${strand}_empiric_per.out

	# Parsing output file
	$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "${strand}_exponential_PMD.txt" --fixedSeed 0 --out ${strand}_exponential_reestimate --logFile ${strand}_exponential_reestimate.out
	$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "${strand}_empiric_PMD.txt" ---fixedSeed 0 --out ${strand}_empiric_reestimte --logFile ${strand}_empiric_reestimte.out
	$atlas --task PMD --bam $strand.bam --fasta $strand.fasta --pmd "${strand}_empiric_per_PMD.txt" ---fixedSeed 0 --out ${strand}_empiric_per_reestimte --logFile ${strand}_empiric_per_reestimte.out
done
