#! /bin/bash

. $(dirname $0)/find_atlas

. $(dirname $0)/simulate --numReadGroups 10 --fixedSeed 155

echo "SimReadGroup1 SimReadGroup2" > rgs.txt
echo "SimReadGroup3 SimReadGroup4 SimReadGroup5 SimReadGroup6 SimReadGroup7 SimReadGroup8 SimReadGroup9 SimReadGroup10" >> rgs.txt

out="mergeRG"
$atlas --task mergeRG --bam simulate.bam --readGroups rgs.txt \
	   --fixedSeed 159 --out $out --logFile $out.out 2> $out.eout
