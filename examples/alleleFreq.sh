#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --vcf --type HW --sampleSize 2 --chrLength 1000 --fixedSeed 0 --logFile simulate.out

$atlas --task alleleFreq --vcf ATLAS_simulations.vcf.gz --out alleleFreq --logFile alleleFreq.out
$atlas --task alleleFreq --likelihoods --vcf ATLAS_simulations.vcf.gz --out alleleFreqLKs --logFile alleleFreqLKs.out
#$atlas --task alleleFreq --compare --vcf ATLAS_simulations.vcf.gz --samples alleleFreq_alleleFreq.txt.gz --out compare --logFile compare.out
#someone able to set up a working test?
