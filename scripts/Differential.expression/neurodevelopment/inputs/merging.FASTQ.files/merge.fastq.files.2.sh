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

# Merge neuronal differentiation NPC replicates from three libraries

cat Library1.barcode05.fastq.gz Library2.barcode05.fastq.gz Library3.barcode05.fastq.gz > Library1.2.3.barcode05.merged.fastq.gz

cat Library1.barcode08.fastq.gz Library2.barcode08.fastq.gz Library3.barcode08.fastq.gz > Library1.2.3.barcode08.merged.fastq.gz

cat Library1.barcode11.fastq.gz Library2.barcode11.fastq.gz Library3.barcode11.fastq.gz > Library1.2.3.barcode11.merged.fastq.gz

# Merge neuronal differentiation Day2 replicates from three libraries

cat Library1.barcode02.fastq.gz Library2.barcode02.fastq.gz Library3.barcode02.fastq.gz > Library1.2.3.barcode02.merged.fastq.gz

cat Library1.barcode06.fastq.gz Library2.barcode06.fastq.gz Library3.barcode06.fastq.gz > Library1.2.3.barcode06.merged.fastq.gz

cat Library1.barcode09.fastq.gz Library2.barcode09.fastq.gz Library3.barcode09.fastq.gz > Library1.2.3.barcode09.merged.fastq.gz

# Merge neuronal differentiation Day6 replicates from three libraries

cat Library1.barcode03.fastq.gz Library2.barcode03.fastq.gz Library3.barcode03.fastq.gz > Library1.2.3.barcode03.merged.fastq.gz

cat Library1.barcode07.fastq.gz Library2.barcode07.fastq.gz Library3.barcode07.fastq.gz > Library1.2.3.barcode07.merged.fastq.gz

cat Library1.barcode10.fastq.gz Library2.barcode10.fastq.gz Library3.barcode10.fastq.gz > Library1.2.3.barcode10.merged.fastq.gz
exit


