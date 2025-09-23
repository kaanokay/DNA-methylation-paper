#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=example@hi.is # for example uname@hi.is
#SBATCH --partition=mimir  # request node from a specific partition
#SBATCH --nodes=1                 # number of nodes
#SBATCH --ntasks-per-node=32      # 48 cores per node (96 in total)
#SBATCH --mem-per-cpu=3900        # MB RAM per cpu core
#SBATCH --time=14-00:00:00         # run for 4 hours maximum (DD-HH:MM:SS)
#SBATCH --output=slurm_job_output_%j.log
#SBATCH --error=slurm_job_errors_%j.log # Logs if job crashes

ml load minimap2/2.24-GCCcore-11.2.0
ml load SAMtools/1.15-GCC-11.2.0

minimap2 -t 128 \
-uf \
-ax splice -y \
UCSC_Mus_musculus.GRCm38_genome.fa \
Library1.2.3.barcode05.merged.fastq | samtools view -h -b -q 40 -F 2304 --threads 128 -Sb - | samtools sort - --threads 128 \
> Library1.2.3.barcode05.bam

minimap2 -t 128 \
-uf \
-ax splice -y \
UCSC_Mus_musculus.GRCm38_genome.fa \
Library1.2.3.barcode08.merged.fastq | samtools view -h -b -q 40 -F 2304 --threads 128 -Sb - | samtools sort - --threads 128 \
> Library1.2.3.barcode08.bam

minimap2 -t 128 \
-uf \
-ax splice -y \
UCSC_Mus_musculus.GRCm38_genome.fa \
Library1.2.3.barcode11.merged.fastq | samtools view -h -b -q 40 -F 2304 --threads 128 -Sb - | samtools sort - --threads 128 \
> Library1.2.3.barcode11.bam
exit
