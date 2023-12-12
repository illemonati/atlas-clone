#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --type HW --F 0.1 \
  --sampleSize 19 --chrLength 1111 --fracPoly 1.0 \
  --alpha 2.0 --beta 2.0 --seqType paired --seqCycles 101

$atlas --task filterBAM --bam ATLAS_simulations_ind1.bam \
	   --fixedSeed 0 --filterSoftClips --out filterBAM_1 --logFile filterBAM_1.out

for i in {2..19}; do
	samtools view ATLAS_simulations_ind$i.bam | head -250 | tail -10 | cut -f1 > blacklist_$i.txt
	u=$(echo "$i*5" | bc)
	$atlas --task filterBAM \
		   --bam ATLAS_simulations_ind$i.bam \
		   --filterMQ 0,$u --blacklist blacklist_$i.txt --filterReadLength 0,$u \
		   --filterFragmentLength 0,$u --filterSoftClips "0.$i" \
		   --fixedSeed 0 --out filterBAM_$i --logFile filterBAM_$i.out
done

