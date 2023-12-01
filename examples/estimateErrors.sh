#! /bin/bash

. $(dirname $0)/find_atlas

recal="intercept[0.0];quality:polynomial[0.9,0.01]"
pmd="CT5:0.2*exp(-0.3*p)+0.01;GA3:0.5*exp(-0.2*p)+0.01"
recal=""
pmd=""
recalModel="intercept;quality"

delta=10
k="000"
L="100$k"

. $(dirname $0)/simulate --recal $recal --pmd $pmd --baseQuality "unif()[0,93]" \
	--chrLength $L,$L --depth 10 --ploidy 2,1 --numReadGroups 1

rModels=("intercept;quality" "intercept;quality:polynomial3;position:polynomial3;fragmentLength:polynomial3;mappingQuality:polynomial3;context;" "intercept;quality;position;fragmentLength;mappingQuality;context")
for i in {0..2}; do
	name="diplo$i"
	recalModel=${rModels[i]}
	$atlas --task estimateErrors --minDeltaLL $delta --recalModel $recalModel \
		--bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --chr chr2 --ploidy 2 \
		--fixedSeed 0 --out $name --logFile $name.out

	name="haplo$i"
	$atlas --task estimateErrors --minDeltaLL $delta --recalModel $recalModel \
		--bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --chr chr2  --ploidy 1 \
		--fixedSeed 0 --out $name --logFile $name.out
done
