#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --vcf --type HW --sampleSize 50 --chrLength 1000 --fixedSeed 0 --logFile simulate.out

$atlas inbreeding --vcf ATLAS_simulations.vcf.gz --numBurnin 1 --iterations 100 --fixedSeed 0 --logFile inbreeding.out

