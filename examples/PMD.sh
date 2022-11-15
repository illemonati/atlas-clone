#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --pmd "doubleStrand:Exponential[30,0.1,0.1,0.05]:Exponential[40,0.2,0.3,0.07]"  --fixedSeed 0 --logFile simulate.out

$atlas --task PMD --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --pmdModels "doubleStrand:Exponential:Exponential"  --fixedSeed 0 --out exponential --logFile exponential.out
$atlas --task PMD --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --pmdModels "doubleStrand:Empiric:Empiric"  --fixedSeed 0 --out empiric --logFile empiric.out

# Parsing output file
$atlas --task PMD --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --pmd "exponential_PMD.txt" --reestimate --fixedSeed 0 --out exponential_reestimate --logFile exponential_reestimate.out
$atlas --task PMD --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --pmd "empiric_PMD.txt" --reestimate --fixedSeed 0 --out empiric_reestimte --logFile empiric_reestimte.out
