. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --type HW --sampleSize 2 --fixedSeed 101

out="GLF1"
$atlas --task GLF --bam simulate_ind1.bam --printAll \
	   --fixedSeed 102 --out $out --logFile $out.out 2> $out.eout

out="GLF2"
$atlas --task GLF --bam simulate_ind2.bam --printAll \
	   --fixedSeed 103 --out $out --logFile $out.out 2> $out.eout

out="geneticDist"
$atlas --task geneticDist --glf GLF1.glf.gz,GLF2.glf.gz \
	   --fixedSeed 104 --out $out --logFile $out.out 2> $out.eout


