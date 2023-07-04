. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --type HW --sampleSize 2

$atlas --task GLF --bam ATLAS_simulations_ind1.bam --printAll --fixedSeed 0 --out GLF1 --logFile GLF1.out
$atlas --task GLF --bam ATLAS_simulations_ind2.bam --printAll --fixedSeed 0 --out GLF2 --logFile GLF2.out
$atlas --task geneticDist --glf GLF1.glf.gz,GLF2.glf.gz --fixedSeed 0 --out geneticDist --logFile geneticDist.out


