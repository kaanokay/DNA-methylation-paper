#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=example@hi.is # for example uname@hi.is
#SBATCH --partition=mimir  # request node from a specific partition
#SBATCH --nodes=1                 # number of nodes
#SBATCH --ntasks-per-node=20      # 48 cores per node (96 in total)
#SBATCH --mem-per-cpu=3900        # MB RAM per cpu core
#SBATCH --time=14-00:00:00         # run for 4 hours maximum (DD-HH:MM:SS)
#SBATCH --output=slurm_job_output_%j.log
#SBATCH --error=slurm_job_errors_%j.log # Logs if job crashes

ml load SAMtools/1.15-GCC-11.2.0

samtools sort --threads 128 -o Ctr3_1.bam Ctr3_1.sam
samtools sort --threads 128 -o Ctr3_2.bam Ctr3_2.sam
samtools sort --threads 128 -o Ctr3_3.bam Ctr3_3.sam
samtools sort --threads 128 -o Ctr3_4.bam Ctr3_4.sam
samtools sort --threads 128 -o Dnmt1_1.bam Dnmt1_1.sam
samtools sort --threads 128 -o Dnmt1_2.bam Dnmt1_2.sam
samtools sort --threads 128 -o Dnmt1_3.bam Dnmt1_3.sam
samtools sort --threads 128 -o Dnmt1_4.bam Dnmt1_4.sam
samtools sort --threads 128 -o Kmt2a_1.bam Kmt2a_1.sam
samtools sort --threads 128 -o Kmt2a_2.bam Kmt2a_2.sam
samtools sort --threads 128 -o Kmt2a_3.bam Kmt2a_3.sam
samtools sort --threads 128 -o Kmt2a_4.bam Kmt2a_4.sam
exit
