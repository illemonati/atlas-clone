#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File with exponential pmd pattern
$atlas simulate --logFile simulate.out --pmd "CT5:0.2*exp(-0.3*p)+0.07;GA3:0.2*exp(-0.25*p)+0.06"

# estimate PMD pattern
$atlas estimateErrors --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --logFile estimateErrors.out --out ee

# Create BAM file with PMDS scores using PMD pattern in 'ee_RGInfo.json'
$atlas PMDS --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --RGInfo ee_RGInfo.json --logFile PMDS.out

# filter file with PMDS score > 1
$atlas filterBAM --filterPMDS 1 --bam ATLAS_simulations_PMDS.bam --logFile filterBAM.out
