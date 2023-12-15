#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf

$atlas --task alleleCounts --vcf ATLAS_simulations.vcf.gz  --outFormat withAlleles \
	   --logFile alleleCounts.out --fixedSeed 0

$atlas --task ancestralAlleles --alleleCounts ATLAS_simulations_alleleCounts.txt.gz \
	   --minorCountMaximum 2 --totalCountMinimum 4 \
	   --logFile ancestralAlleles.out --fixedSeed 0
