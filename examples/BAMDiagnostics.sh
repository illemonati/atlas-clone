#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --type HW --F 0.1 --chrLength 1000 --sampleSize 20 --fracPoly 1.0 --alpha 2.0 --beta 2.0 --seqType single --cycles 100 --fragmentLength "gamma(10,0.2)[30,200]" --fixedSeed 0 --logFile simulate.out

for i in {1..20}; do
samtools view ATLAS_simulations_ind"$i".bam | head -250 | tail -10 | cut -f1 > blacklist_"$i".txt
u=$(echo "$i*5" | bc)
$atlas --task BAMDiagnostics --bam ATLAS_simulations_ind$i.bam --filterSoftClips --fixedSeed 0 --filterMQ 0,$u --blacklist blacklist_$i.txt --filterReadLength 0,$u --filterFragmentLength 0,$u --logFile BAMDiagnostics_$i.out
done
