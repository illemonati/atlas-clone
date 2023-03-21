#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --vcf --type HW --sampleSize 2 --fixedSeed 0 --logFile simulate.out

$atlas --task alleleCounts --vcf ATLAS_simulations.vcf.gz --out alleleCounts --logFile alleleCounts.out
$atlas --task alleleCounts --likelihoods --vcf ATLAS_simulations.vcf.gz --out alleleCountsLKs --logFile alleleCountsLKs.out
$atlas --task alleleCounts --transform --countsFile alleleCounts_alleleCounts.txt.gz --out transform --logFile transform.out
