#! /bin/bash

. $(dirname $0)/find_atlas

# paired end
. $(dirname $0)/simulate --numReadGroups 10
echo "SimReadGroup1 SimReadGroup2" > rgs.txt
echo "SimReadGroup3 SimReadGroup4 SimReadGroup5 SimReadGroup6 SimReadGroup7 SimReadGroup8 SimReadGroup9 SimReadGroup10" >> rgs.txt

$atlas --task mergeRG --bam ATLAS_simulations.bam --readGroups rgs.txt --logFile mergeRG.out --fixedSeed 0 --out mergeRG
