#! /bin/bash

root=`git rev-parse --show-toplevel`
atlas=`find $root -type f -name atlas | tail -n 1`
echo "Using $atlas"

$atlas --task simulate --type HW --F 0.1 --chrLength 1000 --sampleSize 20 --fracPoly 1.0 --alpha 2.0 --beta 2.0 --readLength "single:gamma(10,0.2)[30,101]" --fixedSeed 0

for i in {1..20}; do
samtools view ATLAS_simulations_ind"$i".bam | head -250 | tail -10 | cut -f1 > blacklist_"$i".txt
u=$(echo "$i*5" | bc)
$atlas --task filterBAM --bam ATLAS_simulations_ind"$i".bam --filterSoftClips --fixedSeed 0 --filterMQ 0,$u --blacklist blacklist_"$i".txt --filterReadLength 0,$u --filterFragmentLength 0,$u
done

