#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File with 3 readgroups
$atlas simulate --logFile simulate.out --numReadGroups 3

# Create merge-file which merges SimReadGroup2 into SimReadGroup1
echo "receiver donor" > rgs.txt
echo "SimReadGroup1 SimReadGroup2" >> rgs.txt

# Merge readgroups
$atlas mergeRG --bam ATLAS_simulations.bam --logFile mergeRG.out --readGroups rgs.txt
