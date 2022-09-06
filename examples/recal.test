#! /bin/bash

. $(dirname $0)/find_atlas

# see python ../tools/plotRecal.py "0.1 + 0.8*x + 0.2*x**2 + 0.01*x**3"
r_model="intercept[0.0];quality:polynomial()[0.9,0.01];rho[[-,1.,0.,0.];[0.,-,1.,0.];[0.,0.,-,1.];[1.,0.,0.,-]]"
#r_model="intercept[0.0];quality:polynomial()[0.95,0.01]"
p_model="doubleStrand:Exponential[30,0.3,0.3,0.5]:Exponential[30,0.3,0.3,0.5]" 
# Simululate polynomial model

#$atlas --task simulate --recal1 $r_model --recal2 $r_model --pmd $p_model --fixedSeed 0 --logFile simulate.out
$atlas --task simulate --recal1 $r_model --recal2 $r_model --fixedSeed 0 --logFile simulate.out

#$atlas --task PMD --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --pmdModels "doubleStrand:Exponential:Exponential"  --fixedSeed 0 --logFile PMD.out

# estimate recal model using polynomial
$atlas --task recal --bam ATLAS_simulations.bam --recal "intercept;quality:polynomial(2)" --rerecalibrate --minDeltaLL 10 --fixedSeed 0 --logFile recal_polynomial.out
#$atlas --task recal --bam ATLAS_simulations.bam --recal "intercept;quality:polynomial(2)" --rerecalibrate --pmd $p_model --minDeltaLL 10 --fixedSeed 0 --logFile recal_polynomial.out

