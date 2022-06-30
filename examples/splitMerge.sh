#! /bin/bash

root=`git rev-parse --show-toplevel`
atlas=`find $root -type f -name atlas | tail -n 1`
echo "Using $atlas"

$atlas --task simulate --numReadGroups 10 --fixedSeed 0 --readLength "single:gamma(75,1)[30,150]"
echo "readGroup seqType maxCycles" > rgs.txt
echo "SimReadGroup1 single 75" >> rgs.txt
echo "SimReadGroup3 single 75" >> rgs.txt
echo "SimReadGroup5 single 75" >> rgs.txt
echo "SimReadGroup7 single 75" >> rgs.txt
echo "SimReadGroup9 single 75" >> rgs.txt
$atlas --task splitMerge --bam ATLAS_simulations.bam --readGroupSettings rgs.txt
