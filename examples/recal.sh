#! /bin/bash

. $(dirname $0)/find_atlas

# see python ../tools/plotRecal.py "0.1 + 0.8*x + 0.2*x**2 + 0.01*x**3"
model="intercept[0.1];quality:polynomial[0.8,0.2,0.01]"

# Simululate polynomial model
$atlas --task simulate --recal $model --fixedSeed 0 --logFile simulate.out

printf "#%-10s %s\n" "LL" "model" > LL.txt

# Calculate log Likelihood of model given simulated data
$atlas --task recal --bam ATLAS_simulations.bam --recal $model --onlyLL --fixedSeed 0 --logFile onlyLL.out
printf "%.4e %s\n" $(grep "Log Likelihood" onlyLL.out | tail -n 1 | awk '{print $5}') "simulation" >> LL.txt

# estimate recal model using polynomial
$atlas --task recal --bam ATLAS_simulations.bam --recal "intercept;quality:polynomial3" --minDeltaLL 1e5 --fixedSeed 0 --out polynomial --logFile polynomial.out
printf "%.4e %s\n" $(grep "Log Likelihood" polynomial.out | tail -n 1 | awk '{print $6}') "polynomial"  >> LL.txt
$atlas --task recal --bam ATLAS_simulations.bam --recal "polynomial_recal.txt" --onlyLL --fixedSeed 0 --out polynomial_read --logFile polynomial_read.out
printf "%.4e %s\n" $(grep "Log Likelihood" polynomial_read.out | tail -n 1 | awk '{print $5}') "polynomial_read"  >> LL.txt
$atlas --task recal --bam ATLAS_simulations.bam --RGInfo "polynomial_recal.json" --onlyLL --fixedSeed 0 --out polynomial_read_json --logFile polynomial_read_json.out
printf "%.4e %s\n" $(grep "Log Likelihood" polynomial_read_json.out | tail -n 1 | awk '{print $5}') "polynomial_read_json"  >> LL.txt

# estimate recal model using empiric
$atlas --task recal --bam ATLAS_simulations.bam --recal "intercept;quality;position;context;fragmentLength;mappingQuality;" --minDeltaLL 1e5 --fixedSeed 0 --out empiric --logFile empiric.out
printf "%.4e %s\n" $(grep "Log Likelihood" empiric.out | tail -n 1 | awk '{print $6}') "empiric"  >> LL.txt
$atlas --task recal --bam ATLAS_simulations.bam --recal "empiric_recal.txt" --onlyLL --fixedSeed 0 --out empiric_read --logFile empiric_read.out
printf "%.4e %s\n" $(grep "Log Likelihood" empiric_read.out | tail -n 1 | awk '{print $5}') "empiric_read"  >> LL.txt
$atlas --task recal --bam ATLAS_simulations.bam --RGInfo "empiric_recal.json" --onlyLL --fixedSeed 0 --out empiric_read_json --logFile empiric_read_json.out
printf "%.4e %s\n" $(grep "Log Likelihood" empiric_read_json.out | tail -n 1 | awk '{print $5}') "empiric_read_json"  >> LL.txt
