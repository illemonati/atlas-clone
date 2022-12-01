#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --vcf --type HW --sampleSize 10 --chrLength 100000 --fixedSeed 0 --logFile simulate.out

$atlas --task calculateF2 --vcf ATLAS_simulations.vcf.gz --filterDepth [2,] --maxMissing 0.1 --minMAF 0.01 --logFile calculateF2.out

