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

# To able to run SeqKit from conda envrionment first thing that should be done:

. ~/.bashrc

ml purge

# Load Anaconda3
ml load Anaconda3/2022.05

# Activate featureCounts

conda activate featureCounts

featureCounts Library1.2.3.barcode02.bam Library1.2.3.barcode03.bam Library1.2.3.barcode05.bam Library1.2.3.barcode06.bam Library1.2.3.barcode07.bam Library1.2.3.barcode08.bam Library1.2.3.barcode09.bam Library1.2.3.barcode10.bam Library1.2.3.barcode11.bam \
-a mm10.refGene.gtf \
 -o featureCounts.neurodevelopment.txt \
-g gene_id \
-L \
-T 64 \
--verbose
exit

