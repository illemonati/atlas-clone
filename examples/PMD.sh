#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --pmd "doubleStrand:Exponential[30,0.1,0.1,0.05]:Exponential[40,0.2,0.3,0.07]"  --fixedSeed 0 --logFile simulate.out
$atlas --task PMD --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --pmdModels "doubleStrand:Exponential:Exponential"  --fixedSeed 0 --out exponential --logFile exponential.out
$atlas --task PMD --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --pmdModels "doubleStrand:Empiric:Empiric"  --fixedSeed 0 --out empiric --logFile empiric.out
