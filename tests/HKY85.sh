#! /bin/bash

. $(dirname $0)/find_atlas

recal="intercept[0.1];quality:polynomial[0.9,0.01];position:polynomial[0.01]"
pmd="CT5:0.1*exp(-0.2*p)+0.001;GA3:0.1*exp(-0.2*p)+0.001"

. $(dirname $0)/simulate --chrLength 121234,121098 --ploidy 2,1 --depth 11,9 \
	--type "HKY85" --mu 0.55 --thetaG 0.00033 --thetaR 0.015 --fixedSeed 111 \
	--recal $recal --pmd $pmd --baseQuality "categorical(12,22,27,32,37,41)" \

echo "chr1 0 5" > diplo.bed
echo "chr2 5 10" > haplo.bed
tstep=$(date +%s)
for i in {1..120000}; do
	echo "chr1 ${i}0 ${i}5" >> diplo.bed
	echo "chr2 ${i}5 $((i+1))0" >> haplo.bed
done
echo "Bed creation: $(($(date +%s) - $tstep)) seconds"

probs="0.2"
out=probSites
$atlas --task HKY85 --prob "1.0,$probs"  --minDeltaLL 1 --sample sites \
	   --bam simulate.bam --fasta simulate.fasta --chr chr1 \
	   --recal $recal --pmd $pmd --window 65432 \
	   --fixedSeed 222 --out $out --logFile $out.out 2> $out.eout

out=probReads
$atlas --task HKY85 --depth "50,20,10,5,2" --minDeltaLL 1 --sample reads \
	   --bam simulate.bam --fasta simulate.fasta --chr chr1 \
	   --recal $recal --pmd $pmd --window 65432 \
	   --fixedSeed 222 --out $out --logFile $out.out 2> $out.eout

out=depthSites
$atlas --task HKY85 --genomeWide 2 --depth "10,5,2" --sample sites \
	   --bam simulate.bam --fasta simulate.fasta --chr chr1 \
	   --recal $recal --pmd $pmd --minDeltaLL 1 \
	   --fixedSeed 222 --out $out --logFile $out.out 2> $out.eout

out=upToDepth
$atlas --task HKY85 --minDeltaLL 0.1 --genomeWide 2 \
	   --depth "10,5,2,1" --sample upToDepth \
	   --bam simulate.bam --fasta simulate.fasta \
	   --regions diplo.bed --ploidy 2 --recal $recal --pmd $pmd \
	   --fixedSeed 444 --out $out --logFile $out.out 2> $out.eout

# no downsampling
out=haplo
$atlas --task HKY85 --genomeWide --minDeltaLL 0.1 \
	   --bam simulate.bam --fasta simulate.fasta \
	   --regions haplo.bed --ploidy 1 --recal $recal --pmd $pmd \
	   --fixedSeed 555 --out $out --logFile $out.out 2> $out.eout

out=posterior
$atlas --task HKY85 --posterior "0.55,0.015,0.00033" \
	   --bam simulate.bam --fasta simulate.fasta \
	   --fixedSeed 555 --out $out --logFile $out.out 2> $out.eout
