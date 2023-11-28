#! /bin/bash

. $(dirname $0)/find_atlas

recal="intercept[0.0];quality:polynomial[0.9,0.01]"
pmd="CT5:0.2*exp(-0.3*p)+0.01;GA3:0.5*exp(-0.2*p)+0.01"
recal=""
pmd=""
recalModel="intercept;quality"

delta=1e-2
k="000"

. $(dirname $0)/simulate --recal $recal --pmd $pmd --baseQuality "unif()[0,93]" \
	--chrLength 5000$k,5000$k,100$k,5000$k,100$k --depth 10,10,10,2,2 --ploidy 2,1,1,1,1 --numReadGroups 1

name="10xlong"
$atlas --task estimateErrors --minDeltaLL $delta --recalModel $recalModel \
	--bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --chr chr2 \
	--fixedSeed 0 --out $name --logFile $name.out

name="10xshort"
$atlas --task estimateErrors --minDeltaLL $delta --recalModel $recalModel \
	--bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --chr chr3 \
	--fixedSeed 0 --out $name --logFile $name.out

name="2xlong"
$atlas --task estimateErrors --minDeltaLL $delta --recalModel $recalModel \
	--bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --chr chr4 \
	--fixedSeed 0 --out $name --logFile $name.out

name="2xshort"
$atlas --task estimateErrors --minDeltaLL $delta --recalModel $recalModel \
	--bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --chr chr4 \
	--fixedSeed 0 --out $name --logFile $name.out

nameT="0"
$atlas --task theta --prob 1,0.5,0.2,01,0.05,0.02,0.01 \
	--bam ATLAS_simulations.bam \
	--fixedSeed 0 --out $nameT --logFile $nameT.out > /dev/null &


for name in "10xlong" "10xshort" "2xlong" "2xshort"; do
	echo $name
	nameT="${name}T"
	$atlas --task theta --prob 1,0.5,0.2,01,0.05,0.02,0.01 \
	   --bam ATLAS_simulations.bam --RGInfo ${name}_RGInfo.json \
	   --fixedSeed 0 --out $name --logFile $nameT.out > /dev/null &
done
wait


#$atlas --task estimateErrors --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --recalModel "intercept;quality:polynomial3;position:polynomial3;fragmentLength:polynomial3;mappingQuality:polynomial3;context;" --minDeltaLL 100 --fixedSeed 0 --out poly --logFile poly.out
exit

printf "#%-10s %s\n" "LL" "model" > LL.txt
printf "%.4e %s\n" $(grep "Log Likelihood" default.out | tail -n 1 | awk '{print $6}') "default" >> LL.txt
printf "%.4e %s\n" $(grep "Log Likelihood" onlyQual.out | tail -n 1 | awk '{print $6}') "onlyQual" >> LL.txt
printf "%.4e %s\n" $(grep "Log Likelihood" poly.out | tail -n 1 | awk '{print $6}') "poly" >> LL.txt
