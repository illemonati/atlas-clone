#! /bin/bash

root=`git rev-parse --show-toplevel`
atlas=`find $root -type f -name atlas | tail -n 1`
echo "Using $atlas"

$atlas --task simulate --chrLength 100000
$atlas --task downsample --bam ATLAS_simulations.bam --prob 0.5
$atlas --task downsample --bam ATLAS_simulations.bam --prob 0.5,0.3
$atlas --task downsample --bam ATLAS_simulations.bam --prob 0.5,0.3,0.1
