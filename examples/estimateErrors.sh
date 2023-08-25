#! /bin/bash

. $(dirname $0)/find_atlas

# see python ../tools/plotRecal.py "0.1 + 0.8*x + 0.2*x**2 + 0.01*x**3"
recal="intercept[0.1];quality:polynomial[0.8,0.2,0.01]"
pmd="CT5:0.2*exp(-0.3*p)+0.01;CT3:0.1*exp(-0.2*p)+0.01"

# Simulate polynomial model
. $(dirname $0)/simulate --recal $recal --pmd $pmd

$atlas --task estimateErrors --onlyLL --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --recalModel $recal --pmdModel $pmd --fixedSeed 0 --out onlyLL --logFile onlyLL.out

# estimate recal model using default
$atlas --task estimateErrors --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --minDeltaLL 1e2 --fixedSeed 0 --out default --logFile default.out
$atlas --task estimateErrors --onlyLL --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --RGInfo default.json --out default_read --logFile default_read.out

# estimate recal model using polynomial
$atlas --task estimateErrors --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --recalModel "intercept;quality:polynomial3" --pmdModel $pmd --minDeltaLL 1e2 --fixedSeed 0 --out poly --logFile poly.out
$atlas --task estimateErrors --onlyLL --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --RGInfo poly.json --out poly_read --logFile poly_read.out

printf "#%-10s %s\n" "LL" "model" > LL.txt
printf "%.4e %s\n" $(grep "Log Likelihood" onlyLL.out | tail -n 1 | awk '{print $5}') "onlyLL" >> LL.txt
printf "%.4e %s\n" $(grep "Log Likelihood" default.out | tail -n 1 | awk '{print $6}') "default" >> LL.txt
printf "%.4e %s\n" $(grep "Log Likelihood" default_read.out | tail -n 1 | awk '{print $5}') "default_read" >> LL.txt
printf "%.4e %s\n" $(grep "Log Likelihood" poly.out | tail -n 1 | awk '{print $6}') "poly" >> LL.txt
printf "%.4e %s\n" $(grep "Log Likelihood" poly_read.out | tail -n 1 | awk '{print $5}') "poly_read" >> LL.txt
