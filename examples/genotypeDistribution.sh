#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate

$atlas --task genotypeDistribution --HKY85 --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --out genotypeDistribution --logFile genotypeDistribution.out 
