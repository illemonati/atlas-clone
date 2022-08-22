#! /bin/bash

. $(dirname $0)/find_atlas

# see python ../tools/plotRecal.py "0.1 + 0.8*x + 0.2*x**2 + 0.01*x**3"
model="intercept[0.1];quality:polynomial()[0.8,0.2,0.01]"

# Simululate polynomial model
$atlas --task simulate --recal $model --fixedSeed 0 --logFile simulate.out
# Calculate log Likelihood of model given simulated data
$atlas --task recal --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --recal $model --rerecalibrate --onlyLL --fixedSeed 0 --logFile recal_onlyLL.out
# estimate recal model using polynomial
$atlas --task recal --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --recal "intercept;quality:polynomial(3)" --rerecalibrate --minDeltaLL 1 --fixedSeed 0 --out ATLAS_polynomial --logFile recal_polynomial.out
# estimate recal model using empiric
$atlas --task recal --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --recal "intercept;quality:empiric" --rerecalibrate --minDeltaLL 1 --fixedSeed 0 --out ATLAS_empiric --logFile recal_empiric.out

echo "Log likelihood of simulated data:      " $(grep "Log Likelihood" recal_onlyLL.out | tail -n 1 | awk '{print $5}') > ATLAS_ll.txt
echo "Log likelihood of polynomial estimate: " $(grep "Log Likelihood" recal_polynomial.out | tail -n 1 | awk '{print $6}') >> ATLAS_ll.txt
echo "Log likelihood of empiric estimate:    " $(grep "Log Likelihood" recal_empiric.out | tail -n 1 | awk '{print $6}') >> ATLAS_ll.txt
