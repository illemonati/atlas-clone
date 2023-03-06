#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --vcf --type HW --sampleSize 10 --chrLength 1000 --ploidy 2 --fixedSeed 1 --out simulate --logFile simulate.out
$atlas --task VCFFixInt --vcf simulate.vcf.gz --out VCFFixInt --logFile VCFFixInt.out
