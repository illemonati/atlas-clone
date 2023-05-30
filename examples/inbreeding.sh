#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 50 --chrLength 1000

$atlas inbreeding --vcf ATLAS_simulations.vcf.gz --numBurnin 1 --iterations 100 --fixedSeed 0 --logFile inbreeding.out

