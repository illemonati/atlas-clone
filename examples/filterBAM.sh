#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --type HW --F 0.1 --sampleSize 20 --chrLength 1000 --fracPoly 1.0 --alpha 2.0 --beta 2.0 --seqType single --seqCycles 100 --fixedSeed 0 --logFile simulate.out

$atlas --task filterBAM --bam ATLAS_simulations_ind1.bam --fixedSeed 0 --filterSoftClips --logFile filterBAM_1.out
for i in {2..20}; do
	samtools view ATLAS_simulations_ind$i.bam | head -250 | tail -10 | cut -f1 > blacklist_$i.txt
	u=$(echo "$i*5" | bc)
	$atlas --task filterBAM --bam ATLAS_simulations_ind$i.bam --fixedSeed 0 --filterMQ 0,$u --blacklist blacklist_$i.txt --filterReadLength 0,$u --filterFragmentLength 0,$u --filterSoftClips "0.$i" --logFile filterBAM_$i.out
done

