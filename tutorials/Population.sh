#! /bin/bash

shopt -s expand_aliases
alias atlas=$(realpath $(dirname "$0")/../build/atlas)

# ## Data from a Population
#
# Let us first create a directory specific for this tutorial:
#
# ```
mkdir atlas-tutorial_B
cd atlas-tutorial_B
# ```
#
# ### Simulating Population
#
# We will start by simulating a population of 11 samples, adding some Postmortem Damage and Sequencing Errors. We will use the following, non-standard, site frequency spectrum:
#
# ```
echo '#SHAPE=<23>
0.96 0.005 0.01 0.0025 0.003 0.002 0.0016 0.0014 0.0013 0.001 0.001 0.0009 0.0008 0.00075 0.0007 0.00065 0.0006 0.0005 0.0006 0.0005 0.0005 0.0005 0.0004' > sfs.txt
# ```
#
# ```
atlas --task simulate --type SFS --sfs sfs.txt \
	   --sampleSize 11 --chrLength 123456 --depth 13 --ploidy 2 \
	  --recal "intercept[0.1];quality:polynomial[0.95,0.01]" \
	  --pmd "CT5:0.09*exp(-0.2*p)+0.01;GA3:0.11*exp(-0.2*p)+0.01"
# ```
#
# ### Estimating the Errors
#
# First, for each bam-file, we need to estimate the PMD and Sequencing errors. This time, we will estimate the default empiric error model (this may take a few minutes):
# ```
for bam in *.bam; do
	atlas estimateErrors --minDeltaLL 0.1 --fasta ATLAS_simulations.fasta --bam $bam
done
# ```
# This creates 11 json-files, called ATLAS_simulations_ind{1..11}_RGInfo.json
#
# ### Genotype Likelihood Files
#
# The bam-file format is a per-read file format, however, for most downstream analysis, and especially when comparing different samples, we want to do this per site (position on chromosome). For this, the Genotype Likelihood File Format is convenient, which stores for each site the likelihoods of the ten genotypes: AA, AC, AG, AT, CC, CG, CT, GG, GT, TT.
#
# Now, we can create the GLF files using the task GLF:
#
# ```
for bam in *.bam; do
	name=${bam%.bam}
	atlas GLF --bam $bam --RGInfo ${name}_RGInfo.json
done
# ```
#
# The GLF-files are called `*.glf.gz`, and the index-files `*.glf.idx`. You can open the index-files with a text editor and verify that they all have the same number of chromosomes and the same length.
#
# You can print a GLF-file using the following command:
# ```
atlas printGLF --glf ATLAS_simulations_ind5.glf.gz | less
# ```
# The values of the likelihoods are 16-bit Phreds `= -1000*log10(P))`, so a value of 301 would mean a log10Probability of -0.301 or a Probability of 0.5.
#
# We will define a variable that holds all the glf-files:
# ```
samples=$(ls -1 *.glf.gz | paste -s -d ',' -)
# ```
#
# ### Major and Minor Alleles
#
# Now that we have the GLF-files, we can do various population tasks on them. Let's start by Estimating major and minor alleles with the task `majorMinor`
#
# ```
atlas majorMinor --glf $samples
# ```
#
# This creates the Variant Call Format file called ATLAS_majorMinor.vcf.gz
#
# ### Site Allele Frequencies
#
# We simulated the population with a given SFS, so let's see if we can estimate it. First, we need to estimate the site allele frequencies using the task `saf`:
#
# ```
atlas saf --glf $samples --fasta ATLAS_simulations.fasta
# ```
#
# This creates saf-files. If you wish to view it, you can use the command `realSFS print` of the angsd tool chain.
#
#
# ### Site Frequency Spectrum
# To infer the site frequency spectrum of our population, we will use the tool winsfs. Please follow the following instructions to install it:
#
# https://github.com/malthesr/winsfs?tab=readme-ov-file#installation
#
# Please set `winsfs` as an alias to the correct destination, which is usually as follows:
#
# ```
alias winsfs="~/.cargo/bin/winsfs"
# ```
#
# We can calculate and view the normalized site frequency spectrum using the same command:
# ```
winsfs saf.saf.idx | winsfs view --normalise
# ```
#
# Is it similar to the one we simulated? What would the SFS look like without correcting the PMD and sequencing errors?
