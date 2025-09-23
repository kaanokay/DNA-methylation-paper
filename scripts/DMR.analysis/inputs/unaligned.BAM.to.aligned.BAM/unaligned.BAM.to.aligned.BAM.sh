#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=example@hi.is # for example uname@hi.is
#SBATCH --partition=mimir # request node from a specific partition
#SBATCH --nodes=1 # number of nodes
#SBATCH --ntasks-per-node=32 # 32 cores per node
#SBATCH --mem-per-cpu=3900 # MB RAM per cpu core
#SBATCH --time=14-00:00:00 # run for 14 days maximum (DD-HH:MM:SS)
#SBATCH --output=slurm_job_output_%j.log
#SBATCH --error=slurm_job_errors_%j.log # Logs if job crashes

ml use /hpcapps/lib-bio/modules/all

ml load SAMtools/1.17
ml load minimap2/2.26

samtools fastq -T '*' --threads 128 EM_KO.bam | minimap2 -k17 -ax map-ont --secondary=yes -t 128 -y UCSC_Mus_musculus.GRCm38_genome.fa.gz - | samtools view -Sb - --threads 128 | samtools sort - --threads 128 > EM_KO.bam
exit