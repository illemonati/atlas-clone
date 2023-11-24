#! /bin/bash

. $(dirname $0)/find_atlas

# see python ../tools/plotRecal.py "0.9*T(x) + 0.01*T(x)**3"
recal="intercept[0.0];quality:polynomial[0.9,0.01]"
pmd="CT5:0.2*exp(-0.3*p)+0.01;GA3:0.5*exp(-0.2*p)+0.01"

# Simulate polynomial model
. $(dirname $0)/simulate --recal $recal --pmd $pmd --chrLength 100000 --depth 2 --ploidy 1 --numReadGroups 1

# estimate recal model using default
$atlas --task estimateErrors --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --minDeltaLL 100 --fixedSeed 0 --out default --logFile default.out
$atlas --task estimateErrors --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --recalModel "intercept;quality:polynomial3;position:polynomial3;fragmentLength:polynomial3;mappingQuality:polynomial3;context;" --minDeltaLL 100 --fixedSeed 0 --out poly --logFile poly.out

printf "#%-10s %s\n" "LL" "model" > LL.txt
printf "%.4e %s\n" $(grep "Log Likelihood" default.out | tail -n 1 | awk '{print $6}') "default" >> LL.txt
printf "%.4e %s\n" $(grep "Log Likelihood" poly.out | tail -n 1 | awk '{print $6}') "poly" >> LL.txt
