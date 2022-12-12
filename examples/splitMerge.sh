#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --numReadGroups 10 --fixedSeed 0 --seqType paired --logFile simulate.out
echo "readGroup seqType seqCycles" > rgs.txt
echo "SimReadGroup1 paired 100" >> rgs.txt
echo "SimReadGroup3 paired 100" >> rgs.txt
echo "SimReadGroup5 paired 100" >> rgs.txt
echo "SimReadGroup7 paired 100" >> rgs.txt
echo "SimReadGroup9 paired 100" >> rgs.txt
$atlas --task splitMerge --bam ATLAS_simulations.bam --readGroupSettings rgs.txt --logFile splitMerge.out --mergingMethod middle --fixedSeed 0 --out middleMerge
$atlas --task splitMerge --bam ATLAS_simulations.bam --readGroupSettings rgs.txt --logFile splitMerge.out --mergingMethod firstMate --out firstMerge
$atlas --task splitMerge --bam ATLAS_simulations.bam --readGroupSettings rgs.txt --logFile splitMerge.out --mergingMethod secondMate --out secondMerge
#$atlas --task splitMerge --bam ATLAS_simulations.bam --readGroupSettings rgs.txt --logFile splitMerge.out --mergingMethod highestQuality --fixedSeed 0 --out qualityMerge
