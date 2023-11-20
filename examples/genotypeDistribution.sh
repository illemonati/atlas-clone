#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --chrLength 100000

$atlas --task genotypeDistribution --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --out genotypeDistribution --logFile genotypeDistribution.out --window 10000 --HKY85 
$atlas --task theta --bam ATLAS_simulations.bam --out theta --logFile theta.out --window 10000
