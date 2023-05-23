. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --type HW --sampleSize 2

$atlas --task GLF --bam ATLAS_simulations_ind1.bam --printAll --logFile GLF1.out
$atlas --task GLF --bam ATLAS_simulations_ind2.bam --printAll --logFile GLF2.out
$atlas --task geneticDist --glf ATLAS_simulations_ind1.glf.gz,ATLAS_simulations_ind2.glf.gz --iterations 20  --out geneticDist --logFile geneticDist.out


