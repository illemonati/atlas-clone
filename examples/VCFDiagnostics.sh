#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --vcf --type HW --sampleSize 2 --chrLength 2000 --ploidy 2

$atlas --task VCFDiagnostics --vcf ATLAS_simulations.vcf.gz --out VCFDiagnostics --logFile VCFDiagnostics.out 
$atlas --task VCFDiagnostics --fixInt --vcf ATLAS_simulations.vcf.gz --out VCFFixInt --logFile VCFFixInt.out
$atlas --task VCFDiagnostics --writeBED --vcf ATLAS_simulations.vcf.gz --out VCFToInvariantBed --logFile VCFToInvariantBed.out
