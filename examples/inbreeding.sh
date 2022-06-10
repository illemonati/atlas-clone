#! /bin/bash

root=`git rev-parse --show-toplevel`
atlas=`find $root -type f -name atlas | tail -n 1`
echo "Using $atlas"

$atlas --task simulate --vcf --type HW --sampleSize 100 --chrLength 1000 --fixedSeed 0

$atlas inbreeding --vcf ATLAS_simulations.vcf.gz --numBurnin 2 --iterations 1000 --fixedSeed 0

