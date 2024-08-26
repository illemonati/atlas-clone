#! /bin/bash

# `--fixedSeed = N` is needed to have reproducable results in regression test

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --type HW --F 0.1 --fixedSeed 0 \
	--sampleSize 17 --chrLength 11111 --fracPoly 1.0 \
	--alpha 2.0 --beta 2.0 --seqType single --seqCycles 101

for i in {1..17}; do
	samtools view simulate_ind"$i".bam | head -250 | tail -10 | cut -f1 \
			> blacklist_"$i".txt
	u=$(echo "$i*5" | bc)

	out="simple$i"
	$atlas --task BAMDiagnostics --perChromosome --bam simulate_ind$i.bam \
		   --fixedSeed $i --out $out --logFile $out.out 2> $out.eout

	out="complex$i"
	$atlas --task BAMDiagnostics --identifyDuplicates --bam simulate_ind$i.bam \
		   --filterSoftClips --filterMQ 0,$u --blacklist blacklist_$i.txt \
		   --filterReadLength 0,$u --filterFragmentLength 0,$u \
		   --fixedSeed 1$i --out $out --logFile $out.out 2> $out.eout
done
