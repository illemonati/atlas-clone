#! /bin/bash

. $(dirname $0)/find_atlas

# paired end
. $(dirname $0)/simulate --numReadGroups 10 --seqType paired --out paired --logFile paired.out
echo "readGroup seqType seqCycles" > rgs_paired.txt
echo "SimReadGroup1 paired 100" >> rgs_paired.txt
echo "SimReadGroup3 paired 100" >> rgs_paired.txt
echo "SimReadGroup5 paired 100" >> rgs_paired.txt
echo "SimReadGroup7 paired 100" >> rgs_paired.txt
echo "SimReadGroup9 paired 100" >> rgs_paired.txt

$atlas --task splitMerge  --mergingMethod middle \
	   --bam paired.bam --readGroupSettings rgs_paired.txt\
	   --fixedSeed 0 --out middleMerge --logFile middleMerge.out

$atlas --task splitMerge  --mergingMethod firstMate \
	   --bam paired.bam --readGroupSettings rgs_paired.txt \
	   --out firstMerge --logFile firstMerge.out

$atlas --task splitMerge  --mergingMethod secondMate \
	   --bam paired.bam --readGroupSettings rgs_paired.txt \
	   --out secondMerge --logFile secondMerge.out

$atlas --task splitMerge --mergingMethod highestQuality \
	   --bam paired.bam --readGroupSettings rgs_paired.txt \
	   --fixedSeed 0 --out qualityMerge --logFile qualityMerge.out

$atlas --task splitMerge  --mergingMethod randomRead \
	   --bam paired.bam --readGroupSettings rgs_paired.txt \
	   --fixedSeed 0 --out randomMerge --logFile randomMerge.out

# single end
. $(dirname $0)/simulate --numReadGroups 10 --seqType single --out single --logFile single.out

echo "readGroup seqType seqCycles" > rgs_single.txt
echo "SimReadGroup1 single 100" >> rgs_single.txt
echo "SimReadGroup3 single 100" >> rgs_single.txt
echo "SimReadGroup5 single 100" >> rgs_single.txt
echo "SimReadGroup7 single 100" >> rgs_single.txt
echo "SimReadGroup9 single 100" >> rgs_single.txt

$atlas --task splitMerge --bam single.bam --readGroupSettings rgs_single.txt \
	   --fixedSeed 0 --out split --logFile split.out
