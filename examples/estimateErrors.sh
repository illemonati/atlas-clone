#! /bin/bash

. $(dirname $0)/find_atlas

recal="intercept[.0];quality:polynomial[0.95,0.01];position:polynomial[0.02]"
pmd="CT5:0.2*exp(-0.3*p)+0.01;GA3:0.2*exp(-0.3*p)+0.01"

delta=10
k="111"
L="100$k"

. $(dirname $0)/simulate --recal $recal --pmd $pmd  \
	--type "HKY85" --mu 0.55 --thetaG 0.00033 --thetaR 0.015 \
  --chrLength $L,$L --depth 9 --ploidy 2,1 --numReadGroups 3 \
  --baseQuality 'categorical(12:1237270,22:845176,27:1912959,32:21069477,37:34246008,41:339557323)' \
  --fragmentLength 'normal(40,10)[30,101]' \
  --mappingQuality 'normal(60,10)[30,60]'

echo "chr1 0 $L" > bed1.bed
echo "chr2 0 99$k" > bed2.bed

echo "readGroup poolWith" > recal.pool
echo "SimReadGroup2 SimReadGroup1" >> recal.pool

echo "readGroup poolWith" > pmd.pool
echo "SimReadGroup3 SimReadGroup1" >> pmd.pool

recalModel="intercept;quality;position:polynomial1;fragmentLength:polynomial2;mappingQuality;context"
name="estimateErrors"
$atlas --task estimateErrors --minDeltaLL $delta --recalModel $recalModel \
	   --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta \
	   --regions bed1.bed,bed2.bed --ploidy 2,1  --window 45672 \
	   --poolRecal "recal.pool" --poolPMD "pmd.pool" \
	   --fixedSeed 0 --out $name --logFile $name.out

name="fixedMu"
$atlas --task estimateErrors --minDeltaLL $delta --recalModel $recalModel --fixedMu \
	   --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta \
	   --regions bed1.bed,bed2.bed --ploidy 2,1  --window 45672 \
	   --poolRecal "recal.pool" --poolPMD "pmd.pool" \
	   --fixedSeed 0 --out $name --logFile $name.out
