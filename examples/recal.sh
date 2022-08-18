#! /bin/bash

root=`git rev-parse --show-toplevel`
atlas=`find $root -type f -name atlas | tail -n 1`
echo "Using $atlas"

model="intercept[0.1];quality:polynomial()[0.8]"

$atlas --task simulate --recal $model --fixedSeed 0
$atlas --task recal --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --recal $model --rerecalibrate --onlyLL --fixedSeed 0 
#$atlas --task recal --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --recal "intercept;quality:polynomial(1)" --rerecalibrate --minDeltaLL 1 --fixedSeed 0
$atlas --task recal --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --recal "intercept;quality:empiric" --rerecalibrate --minDeltaLL 1 --fixedSeed 0
