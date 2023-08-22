#! /bin/bash

. $(dirname $0)/find_atlas

# see python ../tools/plotRecal.py "0.1 + 0.8*x + 0.2*x**2 + 0.01*x**3"
model="intercept[0.1];quality:polynomial[0.8,0.2,0.01]"

# Simulate polynomial model
. $(dirname $0)/simulate --recal $model

printf "#%-10s %s\n" "LL" "model" > LL.txt

# Calculate log Likelihood of model given simulated data
$atlas --task estimateErrors --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --recal $model --onlyLL --fixedSeed 0 --logFile onlyLL.out
printf "%.4e %s\n" $(grep "Log Likelihood" onlyLL.out | tail -n 1 | awk '{print $5}') "simulation" >> LL.txt

# estimate recal model using polynomial
$atlas --task estimateErrors --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --recalModel "intercept;quality:polynomial3" --minDeltaLL 1e2 --fixedSeed 0 --out polynomial --logFile polynomial.out
printf "%.4e %s\n" $(grep "Log Likelihood" polynomial.out | tail -n 1 | awk '{print $6}') "polynomial"  >> LL.txt
$atlas --task recal --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --RGInfo "polynomial.json" --onlyLL --fixedSeed 0 --out polynomial_read_json --logFile polynomial_read_json.out
printf "%.4e %s\n" $(grep "Log Likelihood" polynomial_read_json.out | tail -n 1 | awk '{print $5}') "polynomial_parse"  >> LL.txt

# estimate recal model using empiric
$atlas --task estimateErrors --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --recalModel "intercept;quality;position;context;fragmentLength;mappingQuality;" --minDeltaLL 1e2 --fixedSeed 0 --out empiric --logFile empiric.out
printf "%.4e %s\n" $(grep "Log Likelihood" empiric.out | tail -n 1 | awk '{print $6}') "empiric"  >> LL.txt
$atlas --task recal --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --RGInfo "empiric.json" --onlyLL --fixedSeed 0 --out empiric_read_json --logFile empiric_read_json.out
printf "%.4e %s\n" $(grep "Log Likelihood" empiric_read_json.out | tail -n 1 | awk '{print $5}') "empiric_parse"  >> LL.txt
