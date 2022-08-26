#! /bin/bash

. $(dirname $0)/find_atlas

# see python ../tools/plotRecal.py "0.1 + 0.8*x + 0.2*x**2 + 0.01*x**3"
model="intercept[0.1];quality:polynomial()[0.8,0.2,0.01]"

# Simululate polynomial model
$atlas --task simulate --recal1 $model --recal2 $model --fixedSeed 0 --logFile simulate.out

# Calculate log Likelihood of model given simulated data
$atlas --task recal --bam ATLAS_simulations.bam --recal $model --rerecalibrate --onlyLL --fixedSeed 0 --logFile recal_onlyLL.out

# estimate recal model using polynomial
$atlas --task recal --bam ATLAS_simulations.bam --recal "intercept;quality:polynomial(3)" --rerecalibrate --minDeltaLL 10 --fixedSeed 0 --out ATLAS_polynomial --logFile recal_polynomial.out

# estimate recal model using empiric
$atlas --task recal --bam ATLAS_simulations.bam --recal "intercept;quality:empiric" --rerecalibrate --minDeltaLL 10 --fixedSeed 0 --out ATLAS_empiric --logFile recal_empiric.out

# Compare likelihoods
printf "#%-10s %s\n" "LL" "model" > ATLAS_ll.txt
printf "%.4e %s\n" $(grep "Log Likelihood" recal_onlyLL.out | tail -n 1 | awk '{print $5}') "simulation" >> ATLAS_ll.txt
printf "%.4e %s\n" $(grep "Log Likelihood" recal_polynomial.out | tail -n 1 | awk '{print $6}') "polynomial"  >> ATLAS_ll.txt
printf "%.4e %s\n" $(grep "Log Likelihood" recal_empiric.out | tail -n 1 | awk '{print $6}') "empiric"  >> ATLAS_ll.txt
