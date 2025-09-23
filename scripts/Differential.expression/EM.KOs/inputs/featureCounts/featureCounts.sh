#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=example@hi.is # for example uname@hi.is
#SBATCH --partition=mimir # request node from a specific partition
#SBATCH --nodes=1 # number of nodes [mimir 66 ve 72 arasında bir node atayacak rastgele]
#SBATCH --ntasks-per-node=32 # 32 cores per node
#SBATCH --mem-per-cpu=3900 # MB RAM per cpu core
#SBATCH --time=14-00:00:00 # run for 4 hours maximum (DD-HH:MM:SS)
#SBATCH --output=slurm_job_output_%j.log
#SBATCH --error=slurm_job_errors_%j.log # Logs if job crashes

# To able to run featureCounts from conda envrionment first thing that should be done:

. ~/.bashrc

ml purge

# Load Anaconda3
ml load Anaconda3/2022.05

# Activate SeqKit environment in conda
conda activate featureCounts

featureCounts Ctr3_1.bam Ctr3_2.bam Ctr3_3.bam Ctr3_4.bam Dnmt1_1.bam Dnmt1_2.bam Dnmt1_3.bam Dnmt1_4.bam Kmt2a_1.bam Kmt2a_2.bam Kmt2a_3.bam Kmt2a_4.bam \
-a mm10.refGene.gtf \
 -o featureCounts.EM.KOs.txt \
-p \
-g gene_id \
-Q 50 \
-T 32 \
--verbose
exit