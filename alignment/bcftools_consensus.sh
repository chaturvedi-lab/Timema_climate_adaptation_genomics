#!/bin/sh
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=samtoolsMpileup
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=samridhi.chaturvedi@gmail.com

    echo ------------------------------------------------------
    echo SLURM: job identifier is $SLURM_JOBID
    echo SLURM: job name is $SLURM_JOB_NAME
    echo ------------------------------------------------------

#SAMTOOLs mpileup version 1.5 options used:
#C = adjust mapping quality; recommended:50, disable:0 [0]
#d = max per-file depth; avoids excessive memory usage [250]
#f = faidx indexed reference sequence file
#q = skip alignments with mapQ smaller than INT [0]
#Q = skip bases with baseQ/BAQ smaller than INT [13]
#g = generate genotype likelihoods in BCF format
#t = --output-tags LIST  optional tags to output:DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR []

#BCFTOOLs call version 1.6 options used
#v = output variant sites only
#c/m = he original calling method (conflicts with -m) or alternative model for multiallelic and rare-variant calling (conflicts with -c)
#p = variant if P(ref|D)<FLOAT with -c [0.5]
#P =  --prior <float-o mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]
#O =  output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' (here it is 'v') 
#o = write output to a file [standard output]

module load samtools
module load bcftools

cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/variantcalling/consensus/

for i in *.gz; do bcftools consensus -f ../../fasta_genome/re_mod_map_timema_06Jun2016_RvNkF702.fasta $i > ${i}.fasta; done




