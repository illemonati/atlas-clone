#! /bin/bash

root=`git rev-parse --show-toplevel`
atlas=`find $root -type f -name atlas | tail -n 1`
echo "Using $atlas"

$atlas --task simulate --numReadGroups 10 --fixedSeed 0
echo "ReadGroup SeqType MaxCycles" > rgs.txt
echo "SimReadGroup1 single 10" >> rgs.txt
echo "SimReadGroup3 single 30" >> rgs.txt
echo "SimReadGroup5 single 50" >> rgs.txt
echo "SimReadGroup7 single 70" >> rgs.txt
echo "SimReadGroup9 single 90" >> rgs.txt
$atlas --task splitMerge --bam ATLAS_simulations.bam --readGroupSettings rgs.txt
