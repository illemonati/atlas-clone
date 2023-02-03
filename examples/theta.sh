#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --fixedSeed 0 --logFile simulate.out
rm ATLAS_simulations.bam.bai # will automatically recreate it

# theta
$atlas --task theta --bam ATLAS_simulations.bam --printAll --fixedSeed 0 --extraVerbose --out standard --logFile standard.out

# theta with downsample
$atlas --task theta --bam ATLAS_simulations.bam --prob 1,0.5,0.2 --fixedSeed 0 --extraVerbose --out downsample --logFile downsample.out 

# theta genomewide
$atlas --task theta --genomeWide --bootstraps 10 --bam ATLAS_simulations.bam --fixedSeed 0 --extraVerbose --out genomewide --logFile genomewide.out

# theta genomewide with downsample
$atlas --task theta  --prob 1,0.5,0.2 --genomeWide --bootstraps 10 --bam ATLAS_simulations.bam --fixedSeed 0 --extraVerbose --out genomewide_downsample --logFile genomewide_downsample.out
