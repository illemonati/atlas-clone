#! /bin/bash

# `--fixedSeed = N` is needed to have reproducable results in regression test

pmdSmall="0.01*exp(-0.1*p)+0.0"
pmdBig="0.2*exp(-0.3*p)+0.07"

. $(dirname $0)/find_atlas

for strand in single double; do
	pmd="CT5:$pmdBig;CT3:$pmdBig"
	[ $strand == "double" ] && pmd="CT5:$pmdSmall"

	. $(dirname $0)/simulate --pmd $pmd  --fixedSeed 3 \
		--out $strand --logFile simulate_$strand.out

	out=filternoPMDS$strand
	$atlas --task filterBAM --bam ${strand}.bam --filterPMDS 0.5 \
		   --fixedSeed 4 --out $out --logFile $out.out 2> $out.eout

	out=PMDS$strand
	$atlas --task PMDS --bam $strand.bam --fasta $strand.fasta --pmd $pmd \
		   --fixedSeed 5 --out $out --logFile $out.out 2> $out.eout

	out=filterPMDS$strand
	$atlas --task filterBAM --bam PMDS${strand}_PMDS.bam --filterPMDS 0.5 \
		   --fixedSeed 6 --out $out --logFile $out.out 2> $out.eout
done
