#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 10

$atlas --task calculateF2 --vcf ATLAS_simulations.vcf.gz --filterDepth [2,] --maxMissing 0.1 --minMAF 0.01 --logFile calculateF2.out

