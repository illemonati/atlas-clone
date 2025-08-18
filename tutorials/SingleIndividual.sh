#! /bin/bash

shopt -s expand_aliases
alias atlas=$(realpath $(dirname "$0")/../build/atlas)

# ## Data from a Single Individual
#
# In this tutorial, you will learn how to interact with ATLAS, how to provide input and interpret output. It will also guide you through a typical workflow when working with ancient or other low-depth data, focusing on the analysis of a single individual. We will then introduce working with multiple individuals in the next tutorial.
#
# ### Simulating some Data
#
# We will begin by simulating sequencing data of a single individual using ATLAS. We will directly simulate a BAM file that contains data from two sequencing runs present in two distinct read groups.
#
# Let us first create a directory specific for this tutorial:
#
# ```
mkdir atlas-tutorial_A
cd atlas-tutorial_A
# ```
# To simulate data from multiple readgroups, ATLAS expects details for those read groups in a read group info file. Let us create such a file as follows:
#
# ```
echo '
{
	"RG_one": {
		"seqType" : "single",
		"mappingQuality" : "fixed(50)",
		"baseQuality" : "unif()[10,30]",
		"seqCycles" : "100"
	},
	"RG_two":{
		"seqType" : "single",
		"mappingQuality" : "normal(50,10)[10,80]",
		"baseQuality" : "poisson(20)[1,40]",
		"seqCycles" : "100"
	}
}' > ReadGroupInfo.json
# ```
# You can look at this file using the command `cat ReadGroupInfo.json`.
#
# It indicates that both read groups were sequenced as single-end (column seqType) and that the second,
# RG_two, has non-fixed mapping qualities (column mappingQuality).
#
# You may now use this file to simulate a BAM file as follows:
# ```
atlas simulate --chrLength 1000000 --depth 50 --RGInfo ReadGroupInfo.json
# ```
# Several output files will be created:
#
# |  |  |
# |--|--|
# | \*.bam |Simulated bam file called ATLAS_simulations.bam. |
# | \*.bam.bai | Index file for the simulated bam file called ATLAS_simulations.bam.bai. |
# | \*.fasta | Simulated reference sequence file called ATLAS_simulations.fasta.|
# | \*.fasta.fai | Index file for the simulated reference sequence file called ATLAS_simulations.fasta.fai. |
# | \*.parameters | File listing all parameters used called simulate.parameters. |
#
# * You may change the prefix form ATLAS_simulations to whatever you like with `--out`.
#
# More information on `simulate` can be found [here](#simulate).
#
#
# ### Diagnosing the Simulated BAM File
#
# You may of course look at the simulated BAM file using samtools. However, BAM files are typically very large, and we may be more interested in summary statistics than the raw sequencing data. You can easily obtain such summary statistics with ATLAS by using the task `BAMDiagnostics` as follows:
#
# ```
atlas BAMDiagnostics --bam ATLAS_simulations.bam
# ```
#
# This is a good time to carefully study the detailed log ATLAS is writing to the command line: it details every step and informs you about every parameter (i.e. argument) used, as well as about additional options that may be used.
#
# * Can you find how you could limit the analysis to reads with high mapping quality only?
#
# ```
atlas BAMDiagnostics --bam ATLAS_simulations.bam --filterBaseQual 1,90
# ```
#
# You may also get a copy of that log written to a file - simply add the parameter `--logFile filename` or append `&> logfile_name.txt` to the above command to indicate the name of the file to which you want the log to be printed. You may also decide to suppress the log with `--silent` when `--logFile filename` is used.
#
# Running BAM diagnostics produces several output files, among which the file `ATLAS_simulations_diagnostics.txt` lists per-read group summary statistics. You may look at it using
# ```
cat ATLAS_simulations_diagnostics.txt
# ```
#
# * Do both readgroups contribute equally to the total depth?
#
# More information on `BAMDiagnostics` can be found [here](#BAMDiagnostics).
#
# ### Filtering Data
#
# ATLAS offers a wide range of filters, which are detailed [here](#Filter) or also in the log.
#
# Try to rerun `BAMDiagnostics` only on the reads mapping in forward strand. Then, check the log for the filter info as well as the filter statistics printed after the BAM file was parsed. These statistics are all provided as a file `ATLAS_simulations_filterSummary.txt`.
#
# ### Calling Genotypes
#
# For most applications, we recommend inferring parameters of interest directly from genotype likelihoods, rather than first calling genotypes and basing subsequent inferences on those. The reason is simply: calling genotypes requires relatively high depth to be accurate enough to ignore their uncertainty downstream.
#
# Nonetheless, genotypes may need to be called for several applications, and ATLAS boosts many different callers for that purpose, starting from random read callers to Bayesian callers using specific priors.
#
# Calling genotypes for a single individual is very easy. Simply type:
# ```
atlas call --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta
# ```
# This will generate a gzipped VCF file `ATLAS_simulations_calls_maximumLikelihood.vcf.gz` with calls the default MLE (Maximum Likelihood Estimation) caller. You can take a look at it using
#
# ```
zcat ATLAS_simulations_calls_maximumLikelihood.vcf.gz | less
# ```
#
# The caller offers plenty of settings, such as which VCF fields shall be printed, if sites without data should be added to the file (default: no), or if genotypes with two alleles different from the reference are allowed (default: yes). All the parameters can be found [here](#call)
#
# If you prefer a different caller, you may set it with `--method`.
#
#
# ### Estimating Heterozygosity
#
# Let us next use ATLAS to infer something interesting: the heterozygosity of the sample along the genome. You can do so in two ways:
#
# * 1) You may call genotypes and count how many are heterozygous. This requires considerable depth, as genotyping errors affect your estimates.
# * 2) You may infer heterozygosity from genotype likelihoods, which account for genotyping uncertainty and should be accurate even at low depth.
#
# Since we simulated data with a high depth (50), let us try both!
#
# Let us first count the number of heterozygous and homozygous calls in the VCF file. You can get the counts and fraction of heterozygous calls using:
#
# ```
zcat ATLAS_simulations_calls_maximumLikelihood.vcf.gz \
	 | awk '$1!~/#/{print $10}' | cut -d':' -f1 \
	 | awk -F '/' '{if($1==$2){++hom}else{++het}}END{print het, hom, het/(het+hom)}'
# ```
# To estimate heterozygosity from genotype likelihoods, we use the task `theta` as follows:
# ```
atlas theta --bam ATLAS_simulations.bam
# ```
# This task estimates $\theta$, the mutation rate under a Felsenstein model that allows for recurrent mutations. The estimate is written to the log and in the file `ATLAS_simulations_theta.txt.gz`. They should be close to 0.001, which is the default theta used in the simulations.
#
# To compare it to the heterozygosity from calls, we need to translate $\theta$ values into expected heterozygosities. ATLAS provides these both in the log as well as in the output file (column `expHet_MLE`).
#
# You should see that these estimates are pretty similar - which is expected given the high depth.
#
# More information on `theta` can be found [here](#theta).
#
# ### Downsampling your Data
#
# Now, let's assess how these two estimates depend on depth. An easy way to do so is to downsample the simulated data and rerun the estimation. You can downsample simulated data with task `downsample` as
#
# ```
atlas downsample --bam ATLAS_simulations.bam --prob 0.5,0.1,0.04
# ```
# Here, the parameter `--prob` specifies the probabilities with which sequencing reads are kept. With `--prob 0.5,0.1,0.04` ATLAS will generate three new BAM files, for which only 50%, 10% or 4% of all reads will be kept, resulting in a depth of about 25, 5 and 2, respectively. You can check with `BAMDiagnostics` if the resulting files `ATLAS_simulations_downsampled_*0.100000*.bam` have the correct depth.
#
# Now repeat the estimation of heterozygosity and $\theta$ for them by both calling and using genotype likelihoods.
#
# * Do they still match?
#
# More information on `downsample` can be found [here](#downsample).
#
# ### Post-Mortem Damage
# #### The problem
#
# The data we simulated above had sequencing errors, but none of the other complications of ancient DNA. One of those complications is post-mortem damage (PMD), which we want to study next.
#
# Let us begin by simulating data with PMD:
# ```
atlas simulate --RGInfo ReadGroupInfo.json --chrLength 1000000 --depth 2 \
	  --pmd "CT5:0.2*exp(-0.3*p)+0.07;GA3:0.1*exp(-0.3*p)+0.2" --out Ancient
# ```
#
# Here, we indicate the PMD-pattern using `--pmd`. The PMD pattern on the 5' end is C->T and follows an exponential curve, the pattern on the 3' end is G->A with slightly different values for the exponential curve. We also used `--out Ancient` to set the output prefix to `Ancient`. The resulting BAM file will thus be named `Ancient.bam`.
#
# This data was again simulated with default $\theta=0.001$. But what is the estimate of $\theta$ from that data? It should be much higher due to PMD being considered true variants.
#
# ```
atlas theta --bam Ancient.bam
# ```
#
# #### The solution
#
# To correct for PMD, we first need to estimate the errors from the data (as we know the true patterns only for simulated data). You can estimate PMD as
#
# ```
atlas estimateErrors --bam Ancient.bam --fasta Ancient.fasta --minDeltaLL 0.1 --out AncientEE
# ```
#
# We set a relatively high deltaLL of 0.1 to make the estimation finish quicker, in reality, you should use a deltaLL of 0.001 or smaller. (`estimateErrors` also estimates the base-quality score recalibration at the same time, see next section.)
#
# These estimates were written to the file `AncientEE_RGInfo.json`.
#
# We can now use the learned pattern to account for PMD when calculating genotype likelihoods with `--RGInfo AncientEE_RGInfo.json`:
#
# ```
atlas theta --bam Ancient.bam --RGInfo AncientEE_RGInfo.json
# ```
#
# Do you now get $\theta$ estimates close to $\theta=0.001$? Note that we only had 2x data for 1Mb with a lot of PMD!
# Note also that the PMD pattern in the json-File works with all tasks in ATLAS, including calls:
#
# ```
atlas call --bam Ancient.bam --fasta Ancient.fasta --RGInfo AncientEE_RGInfo.json
# ```
#
#
# ### Base-Quality Score Recalibration
# #### The Problem
#
# Another issue with most low-depth data, whether ancient or modern, is that the base quality scores provided by the sequencing machines are usually off. That is an issues since the genotype likelihood calculation relies on them.
#
# To illustrate the problem, we can again simulate data, this time with distorted quality scores using
#
# ```
atlas simulate --RGInfo ReadGroupInfo.json --chrLength 1000000 --depth 5 \
	  --recal "intercept[0.1];quality:polynomial[0.9]" --out Distorted
# ```
#
# Here we use a simple distortion where the real quality scores are only 90% of the sequencing machines values. Again, the default $\theta=0.001$ was used to simulate data, but which $\theta$ do you estimate? It should be too high, as the simulations indicate too high base qualities (i.e. too much confidence in the sequencing data).
#
# ```
atlas theta --bam Distorted.bam
# ```
#
# #### The solution
#
# Just as for PMD, we can estimate recal error patterns. As we simulated a polynomial quality-pattern, let's see whether we can reproduce it:
#
# ```
atlas estimateErrors --bam Distorted.bam --fasta Distorted.fasta \
	  --recalModel "intercept;quality:polynomial1" --poolRecal "all" --poolPMD "all" \
	  --out DistortedPoly --minDeltaLL 0.5
# ```
#
# Here we provide the model with `--recalModel "intercept;quality:polynomial1"`, saying we only want to estimate it based on the quality scores and a constant term (intercept). With `poolRecal "all"`, we are pooling the data of the two readgroups togehter, as we did not simulate a different pattern per readgroup. We als use `poolPMD "all"`, atlas will estimate a PMD-pattern, however, it's values will be  close to 0.
#
# You can now use the polynomial recalibration parameters in the json file to rerun theta:
#
# ```
atlas theta --bam Distorted.bam --RGInfo DistortedPoly_RGInfo.json
# ```
#
# With real data, we do not know which covariates have an effect on the recalibration pattern. Per default, Atlas takes into account five covariates: The quality scores (`quality`), the mapping quality (`mappingQuality`), the fragment length (`fragmentLength`), the distance form the 5'-end (`position`) and the previous base (`context`). We can estimate the default recal-pattern as follows:
#
# ```
atlas estimateErrors --bam Distorted.bam --fasta Distorted.fasta \
	  --poolRecal "all" --poolPMD "all" \
	  --out DistortedDefault --minDeltaLL 0.5
# ```
#
# With no `recalModel` given, an empiric estimation will be done for all covariates. The estimates were written to the file `DistortedDefault_RGInfo.json`.
#
# You can also use these recalibration parameters to rerun theta:
#
# ```
atlas theta --bam Distorted.bam --RGInfo DistortedDefault_RGInfo.json
# ```
