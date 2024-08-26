#! /bin/bash

# `--fixedSeed = N` is needed to have reproducable results in regression test

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --fixedSeed 3 --sampleSize 2 --chrLength 2321 --ploidy 2

out="VCFDiagnostics"
$atlas --task VCFDiagnostics --vcf simulate.vcf.gz \
	   --fixedSeed 2 --out $out --logFile $out.out 2> $out.eout

out="VCFFixInt"
$atlas --task VCFDiagnostics --fixInt --vcf simulate.vcf.gz \
	   --fixedSeed 1 --out $out --logFile $out.out 2> $out.eout

out="VCFToInvariantBed"
$atlas --task VCFDiagnostics --writeBED --vcf simulate.vcf.gz \
	   --fixedSeed 0 --out $out --logFile $out.out 2> $out.eout
