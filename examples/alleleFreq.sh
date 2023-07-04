#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 2 --chrLength 1000

$atlas --task alleleFreq --vcf ATLAS_simulations.vcf.gz --fixedSeed 0 --out alleleFreq --logFile alleleFreq.out
$atlas --task alleleFreq --likelihoods --vcf ATLAS_simulations.vcf.gz --fixedSeed 0 --out alleleFreqLKs --logFile alleleFreqLKs.out
#$atlas --task alleleFreq --compare --vcf ATLAS_simulations.vcf.gz --samples alleleFreq_alleleFreq.txt.gz --out compare --logFile compare.out
#someone able to set up a working test?
