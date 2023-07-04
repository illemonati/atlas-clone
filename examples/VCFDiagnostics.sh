#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 2 --chrLength 2000 --ploidy 2

$atlas --task VCFDiagnostics --vcf ATLAS_simulations.vcf.gz --fixedSeed 0 --out VCFDiagnostics --logFile VCFDiagnostics.out 
$atlas --task VCFDiagnostics --fixInt --vcf ATLAS_simulations.vcf.gz --fixedSeed 0 --out VCFFixInt --logFile VCFFixInt.out
$atlas --task VCFDiagnostics --writeBED --vcf ATLAS_simulations.vcf.gz --fixedSeed 0 --out VCFToInvariantBed --logFile VCFToInvariantBed.out
