#!/bin/bash
root=`git rev-parse --show-toplevel`
atlas=`find $root -type f -name atlas | tail -n 1`
echo "Using $atlas"

$atlas --task simulate --sampleSize 20  --type HW --fixedSeed 0

for i in {1..20}; do
samtools view ATLAS_simulations_ind"$i".bam | head -250 | tail -10 | tr -d "M" | awk '{print $3,$4-1,$6+($4-1)}' | tr " " "\t" > mask_"$i".txt
u=$(echo "$i*0.05" | bc)
$atlas --task writeDepth --bam ATLAS_simulations_ind"$i".bam --fixedSeed 0 --mask mask_"$i".txt  --readUpToDepth "$i"
done
