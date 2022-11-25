#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --numReadGroups 10 --fixedSeed 0 --logFile simulate.out
echo "readGroup seqType seqCycles" > rgs.txt
echo "SimReadGroup1 single 75" >> rgs.txt
echo "SimReadGroup3 single 75" >> rgs.txt
echo "SimReadGroup5 single 75" >> rgs.txt
echo "SimReadGroup7 single 75" >> rgs.txt
echo "SimReadGroup9 single 75" >> rgs.txt
$atlas --task splitMerge --bam ATLAS_simulations.bam --readGroupSettings rgs.txt --logFile splitMerge.out
