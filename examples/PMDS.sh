#! /bin/bash

# `--fixedSeed = N` is needed to have reproducable results in regression test

pmd5="0.1*exp(-0.1*p)+0.05"
pmd3="0.2*exp(-0.3*p)+0.07"

. $(dirname $0)/find_atlas

for strand in single double; do
	pmd="CT5:$pmd5;CT3:$pmd3"
	[ $strand == "double" ] && pmd="CT5:$pmd5;GA3:$pmd3"

	. $(dirname $0)/simulate --pmd $pmd  --fixedSeed 7 \
		--out $strand --logFile simulate_$strand.out

	$atlas --task PMDS --bam $strand.bam --fasta $strand.fasta \
		   --fixedSeed 5 --out $strand --logFile $strand.out 2> $strand.eout
done
