#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --vcf --type HW --sampleSize 2 --chrLength 1000 --ploidy 2 --fixedSeed 1 --out simulate --logFile simulate.out

$atlas --task VCFDiagnostics --vcf simulate.vcf.gz --out VCFDiagnostics --logFile VCFDiagnostics.out 
