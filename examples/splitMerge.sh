#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --numReadGroups 10 --fixedSeed 0 --seqType paired --logFile simulate.out
echo "readGroup seqType seqCycles" > rgs.txt
echo "SimReadGroup1 paired 100" >> rgs.txt
echo "SimReadGroup3 paired 100" >> rgs.txt
echo "SimReadGroup5 paired 100" >> rgs.txt
echo "SimReadGroup7 paired 100" >> rgs.txt
echo "SimReadGroup9 paired 100" >> rgs.txt
$atlas --task splitMerge --bam ATLAS_simulations.bam --readGroupSettings rgs.txt --logFile middleMerge.out --mergingMethod middle --fixedSeed 0 --out middleMerge
$atlas --task splitMerge --bam ATLAS_simulations.bam --readGroupSettings rgs.txt --logFile firstMerge.out --mergingMethod firstMate --out firstMerge
$atlas --task splitMerge --bam ATLAS_simulations.bam --readGroupSettings rgs.txt --logFile secondMerge.out --mergingMethod secondMate --out secondMerge
$atlas --task splitMerge --bam ATLAS_simulations.bam --readGroupSettings rgs.txt --logFile qualityMerge.out --mergingMethod highestQuality --fixedSeed 0 --out qualityMerge
$atlas --task splitMerge --bam ATLAS_simulations.bam --readGroupSettings rgs.txt --logFile randomMerge.out --mergingMethod randomRead --fixedSeed 0 --out randomMerge
