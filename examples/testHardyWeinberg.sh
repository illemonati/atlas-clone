#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --vcf --type HW --sampleSize 10 --chrLength 10000 --ploidy 2 --out simulate
$atlas --task testHardyWeinberg --vcf simulate.vcf.gz --out testHardyWeinberg --logFile testHardyWeinberg.out 
