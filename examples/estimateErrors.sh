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
  --baseQuality 'categorical(12:1237,22:845,27:1912,32:21069,37:34246,41:339557)' \
  --fragmentLength 'normal(40,10)[30,101]' \
  --mappingQuality 'normal(60,10)[30,60]' --fixedSeed 80

echo "chr1 0 $L" > bed1.bed
echo "chr2 0 99$k" > bed2.bed

echo "readGroup poolWith" > recal.pool
echo "SimReadGroup2 SimReadGroup1" >> recal.pool

echo "readGroup poolWith" > pmd.pool
echo "SimReadGroup3 SimReadGroup1" >> pmd.pool

out="estimateErrors"
$atlas --task estimateErrors --minDeltaLL $delta --shuffleSites \
	   --bam simulate.bam --fasta simulate.fasta \
	   --poolRecal "all" --poolPMD "all" \
	   --fixedSeed 81 --out $out --logFile $out.out 2> $out.eout

recalModel="intercept;quality;position:polynomial1;fragmentLength:polynomial2;mappingQuality;context"
out="specificModel"
$atlas --task estimateErrors --minDeltaLL $delta --recalModel $recalModel \
	   --bam simulate.bam --fasta simulate.fasta \
	   --regions bed1.bed,bed2.bed --ploidy 2,1 --window 45672 --filePerIteration \
	   --poolRecal "recal.pool" --poolPMD "pmd.pool" \
	   --fixedSeed 82 --out $out --logFile $out.out 2> $out.eout

