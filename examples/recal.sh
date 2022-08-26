#! /bin/bash

. $(dirname $0)/find_atlas

# see python ../tools/plotRecal.py "0.1 + 0.8*x + 0.2*x**2 + 0.01*x**3"
model="intercept[0.1];quality:polynomial()[0.8,0.2,0.01]"
#model="intercept[0.060883];quality:empiric[0:1.513380,1:0.596903,2:0.032218,3:-0.272727,4:-0.465499,5:-0.607963,6:-0.715386,7:-0.787215,8:-0.826125,9:-0.855006,10:-0.879048,11:-0.855106,12:-0.832865,13:-0.792596,14:-0.766182,15:-0.707180,16:-0.646609,17:-0.558899,18:-0.480407,19:-0.384575,20:-0.287502,21:-0.182669,22:-0.069839,23:0.037841,24:0.151272,25:0.260900,26:0.390357,27:0.515628,28:0.629658,29:0.755226,30:0.879580,31:0.992878,32:1.120109,33:1.221322,34:1.335066,35:1.445328,36:1.551645,37:1.652000,38:1.731477,39:1.813064,40:1.868523,41:1.953203,42:1.988808,43:2.033326,44:2.080246,45:2.099640,46:2.085670,47:2.105547,48:2.078233,49:2.069726,50:2.045986,51:1.911400,52:1.866787,53:1.800953,54:1.672536,55:1.587674,56:1.395788,57:1.226206,58:1.013648,59:0.857134,60:0.576769,61:0.358988,62:-0.033440,63:-0.171867,64:-0.742910,65:-0.726355,66:-1.125323,67:-1.826102,68:-1.866754,69:-2.518098,70:-3.707765,71:-4.217433,72:-2.671585,73:-2.907388,74:-7.603084,75:-9.459973,76:-3.259274,77:-8.651831,78:-6.067265,79:-6.722767,80:0.000000,81:0.000000,82:0.000000,83:0.000000,84:0.000000,85:0.000000]"
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
