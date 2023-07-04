#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 2

$atlas --task alleleCounts --vcf ATLAS_simulations.vcf.gz --fixedSeed 0 --out alleleCounts --logFile alleleCounts.out
$atlas --task alleleCounts --dosaf --vcf ATLAS_simulations.vcf.gz --fixedSeed 0 --out alleleCountsSAF --logFile alleleCountsSAF.out
$atlas --task alleleCounts --transform --countsFile alleleCounts_alleleCounts.txt.gz --fixedSeed 0 --out transform --logFile transform.out
