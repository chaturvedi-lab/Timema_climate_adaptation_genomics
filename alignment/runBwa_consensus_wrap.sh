#!/bin/bash
#SBATCH --job-name=BWA_timema
#SBATCH --time=96:00:00 #walltime
#SBATCH --nodes=1 #number of cluster nodes
#SBATCH --account=usubio-kp #PI account
#SBATCH --partition=usubio-kp #specify computer cluster, other option is kinspeak

module load perl

cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/alignments/fastqfiles/bart/
perl ../../../consensus_align/wrap_bwa_bart.pl *.fastq

cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/alignments/fastqfiles/cali/
perl ../../../consensus_align/wrap_bwa_cali.pl *.fastq

cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/alignments/fastqfiles/cris/
perl ../../../consensus_align/wrap_bwa_cris.pl *.fastq

cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/alignments/fastqfiles/chum/
perl ../../../consensus_align/wrap_bwa_chum.pl *.fastq

cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/alignments/fastqfiles/knul/
perl ../../../consensus_align/wrap_bwa_knul.pl *.fastq

cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/alignments/fastqfiles/land/
perl ../../../consensus_align/wrap_bwa_land.pl *.fastq

cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/alignments/fastqfiles/podu/
perl ../../../consensus_align/wrap_bwa_podu.pl *.fastq

cd /uufs/chpc.utah.edu/common/home/u6007910/projects/timema_adaptation/alignments/fastqfiles/popp/
perl ../../../consensus_align/wrap_bwa_popp.pl *.fastq
