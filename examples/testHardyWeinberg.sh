#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 10 --chrLength 10000 --ploidy 2
$atlas --task testHardyWeinberg --vcf ATLAS_simulations.vcf.gz --out testHardyWeinberg --logFile testHardyWeinberg.out 
