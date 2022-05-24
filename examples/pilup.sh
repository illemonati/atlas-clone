#! /bin/bash

root=`git rev-parse --show-toplevel`
atlas=`find $root -type f -name atlas | tail -n 1`
echo "Using $atlas"

$atlas --task simulate --fixedSeed 0
$atlas --task pileup --bam ATLAS_simulations.bam --fields depth,alleles --window 100000 --printAll
