#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --type HW --F 0.1 \
	--sampleSize 17 --chrLength 11111 --fracPoly 1.0 \
	--alpha 2.0 --beta 2.0 --seqType single --seqCycles 101

for i in {1..17}; do
	samtools view ATLAS_simulations_ind"$i".bam | head -250 | tail -10 | cut -f1 > blacklist_"$i".txt
	u=$(echo "$i*5" | bc)

	$atlas --task BAMDiagnostics --bam ATLAS_simulations_ind$i.bam \
		   --out simple$i --perChromosome --logFile simple$i.out

	$atlas --task BAMDiagnostics --bam ATLAS_simulations_ind$i.bam \
		--filterSoftClips --fixedSeed 0 --filterMQ 0,$u --identifyDuplicates \
		--blacklist blacklist_$i.txt --filterReadLength 0,$u --filterFragmentLength 0,$u \
		--out complex$i --logFile complex$i.out
done
