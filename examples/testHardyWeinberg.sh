#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --vcf --type HW --sampleSize 10 --chrLength 10000 --ploidy 2 --fixedSeed 0 --out simulate --logFile simulate.out
$atlas --task testHardyWeinberg --vcf simulate.vcf.gz --out testHardyWeinberg --logFile testHardyWeinberg.out 
