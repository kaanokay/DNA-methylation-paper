#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=example@hi.is # for example uname@hi.is
#SBATCH --partition=mimir # request node from a specific partition
#SBATCH --nodes=2 # number of nodes
#SBATCH --ntasks-per-node=32 # 32 cores per node
#SBATCH --mem-per-cpu=3900 # MB RAM per cpu core
#SBATCH --time=14-00:00:00 # run for 14 days maximum (DD-HH:MM:SS)
#SBATCH --output=slurm_job_output_%j.log
#SBATCH --error=slurm_job_errors_%j.log # Logs if job crashes

# Load modules

module load SAMtools/1.15-GCC-11.2.0
module load modkit

modkit \
pileup \
EM_KO.bam \
EM_KO.bed \
--ref UCSC_Mus_musculus.GRCm38_genome.fa \
--cpg \
--combine-strands \
--only-tabs \
--threads 128


modkit \
pileup \
control.bam \
control.bed \
--ref UCSC_Mus_musculus.GRCm38_genome.fa \
--cpg \
--combine-strands \
--only-tabs \
--threads 128

exit