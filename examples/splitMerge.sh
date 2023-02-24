#! /bin/bash

. $(dirname $0)/find_atlas

# paired end
$atlas --task simulate --numReadGroups 10 --fixedSeed 0 --seqType paired --out paired --logFile paired.out
echo "readGroup seqType seqCycles" > rgs_paired.txt
echo "SimReadGroup1 paired 100" >> rgs_paired.txt
echo "SimReadGroup3 paired 100" >> rgs_paired.txt
echo "SimReadGroup5 paired 100" >> rgs_paired.txt
echo "SimReadGroup7 paired 100" >> rgs_paired.txt
echo "SimReadGroup9 paired 100" >> rgs_paired.txt

$atlas --task splitMerge --bam paired.bam --readGroupSettings rgs_paired.txt --logFile middleMerge.out --mergingMethod middle --fixedSeed 0 --out middleMerge
$atlas --task splitMerge --bam paired.bam --readGroupSettings rgs_paired.txt --logFile firstMerge.out --mergingMethod firstMate --out firstMerge
$atlas --task splitMerge --bam paired.bam --readGroupSettings rgs_paired.txt --logFile secondMerge.out --mergingMethod secondMate --out secondMerge
$atlas --task splitMerge --bam paired.bam --readGroupSettings rgs_paired.txt --logFile qualityMerge.out --mergingMethod highestQuality --fixedSeed 0 --out qualityMerge
$atlas --task splitMerge --bam paired.bam --readGroupSettings rgs_paired.txt --logFile randomMerge.out --mergingMethod randomRead --fixedSeed 0 --out randomMerge

# single end
$atlas --task simulate --numReadGroups 10 --fixedSeed 0 --seqType single --logFile single.out --out single
echo "readGroup seqType seqCycles" > rgs_single.txt
echo "SimReadGroup1 single 100" >> rgs_single.txt
echo "SimReadGroup3 single 100" >> rgs_single.txt
echo "SimReadGroup5 single 100" >> rgs_single.txt
echo "SimReadGroup7 single 100" >> rgs_single.txt
echo "SimReadGroup9 single 100" >> rgs_single.txt

$atlas --task splitMerge --bam single.bam --readGroupSettings rgs_single.txt --logFile split.out --fixedSeed 0 --out split
