#! /bin/bash

root=`git rev-parse --show-toplevel`
atlas=`find $root -type f -name atlas | tail -n 1`
echo "Using $atlas"

$atlas --task simulate --fixedSeed 0
$atlas --task GLF --bam ATLAS_simulations.bam --printAll
