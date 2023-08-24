#! /bin/bash

. $(dirname $0)/find_atlas

# see python ../tools/plotRecal.py "0.1 + 0.8*x + 0.2*x**2 + 0.01*x**3"
model="intercept[0.1];quality:polynomial[0.8,0.2,0.01]"
pmd="CT5:0.2*exp(-0.3*p)+0.01;CT3:0.1*exp(-0.2*p)+0.01"

# Simulate polynomial model
. $(dirname $0)/simulate --recal $model

printf "#%-10s %s\n" "LL" "model" > LL.txt

# estimate recal model using polynomial
#$atlas --task estimateErrors --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --minDeltaLL 1e2 --fixedSeed 0 --out default --logFile default.out
$atlas --task estimateErrors --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --recalModel "intercept;quality:polynomial3" --pmdModel $pmd --minDeltaLL 1e2 --fixedSeed 0 --out poly --logFile poly.out
