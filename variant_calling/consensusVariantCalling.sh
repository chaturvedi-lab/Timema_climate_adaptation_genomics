#!/bin/sh
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --account=usubio-kp
#SBATCH --partition=usubio-kp
#SBATCH --job-name=varcall
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
#q = skip consensus_align with mapQ smaller than INT [0]
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
### job1 for bart ########
cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/sambamfiles/bart

samtools mpileup -C 50 -d 250 -f /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/references/filtered2x_variantsTimemaBart.vcf.gz.fasta -q 20 -Q 15 -g -I -t DP,DPR -u -b bartSortedBam.txt -o con_variantsTimemaBart.bcf
bcftools call -v -c -p 0.01 -P 0.001 -O v -o con_variantsTimemaBart.vcf con_variantsTimemaBart.bcf 

### job2 for cali ########
cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/sambamfiles/cali

samtools mpileup -C 50 -d 250 -f /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/references/filtered2x_variantsTimemaCali.vcf.gz.fasta -q 20 -Q 15 -g -I -t DP,DPR -u -b caliSortedBam.txt -o con_variantsTimemaCali.bcf
bcftools call -v -c -p 0.01 -P 0.001 -O v -o con_variantsTimemaCali.vcf con_variantsTimemaCali.bcf 

### job3 for chum ########
cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/sambamfiles/chum

samtools mpileup -C 50 -d 250 -f /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/references/filtered2x_variantsTimemaChum.vcf.gz.fasta -q 20 -Q 15 -g -I -t DP,DPR -u -b chumSortedBam.txt -o con_variantsTimemaChum.bcf
bcftools call -v -c -p 0.01 -P 0.001 -O v -o con_variantsTimemaChum.vcf con_variantsTimemaChum.bcf 

### job4 for cris ########
cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/sambamfiles/cris

samtools mpileup -C 50 -d 250 -f /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/references/filtered2x_variantsTimemaCris.vcf.gz.fasta -q 20 -Q 15 -g -I -t DP,DPR -u -b  crisSortedBam.txt -o con_variantsTimemaCris.bcf
bcftools call -v -c -p 0.01 -P 0.001 -O v -o con_variantsTimemaCris.vcf con_variantsTimemaCris.bcf 

### job5 for knul ########
cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/sambamfiles/knul
54
samtools mpileup -C 50 -d 250 -f /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/references/filtered2x_variantsTimemaKnul.vcf.gz.fasta -q 20 -Q 15 -g -I -t DP,DPR -u -b knulSortedBam.txt -o con_variantsTimemaKnul.bcf
bcftools call -v -c -p 0.01 -P 0.001 -O v -o con_variantsTimemaKnul.vcf con_variantsTimemaKnul.bcf

### job6 for land ########
cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/sambamfiles/land

samtools mpileup -C 50 -d 250 -f /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/references/filtered2x_variantsTimemaLand.vcf.gz.fasta -q 20 -Q 15 -g -I -t DP,DPR -u -b landSortedBam.txt -o con_variantsTimemaLand.bcf
bcftools call -v -c -p 0.01 -P 0.001 -O v -o con_variantsTimemaLand.vcf con_variantsTimemaLand.bcf 

### job7 for podu ########
cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/sambamfiles/podu

samtools mpileup -C 50 -d 250 -f /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/references/filtered2x_variantsTimemaPodu.vcf.gz.fasta -q 20 -Q 15 -g -I -t DP,DPR -u -b poduSortedBam.txt -o con_variantsTimemaPodu.bcf
bcftools call -v -c -p 0.01 -P 0.001 -O v -o con_variantsTimemaPodu.vcf  con_variantsTimemaPodu.bcf

### job 8 for popp ########
cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/sambamfiles/popp

samtools mpileup -C 50 -d 250 -f /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/consensus_align/references/filtered2x_variantsTimemaPopp.vcf.gz.fasta -q 20 -Q 15 -g -I -t DP,DPR -u -b poppSortedBam.txt -o con_variantsTimemaPopp.bcf
bcftools call -v -c -p 0.01 -P 0.001 -O v -o con_variantsTimemaPopp.vcf con_variantsTimemaPopp.bcf  
 



