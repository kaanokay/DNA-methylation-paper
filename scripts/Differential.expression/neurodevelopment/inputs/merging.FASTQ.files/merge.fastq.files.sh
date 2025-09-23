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

# Merge fastq files in library 1

cat barcode02/*.fastq.gz > Library1.barcode02.fastq.gz

cat barcode03/*.fastq.gz > Library1.barcode03.fastq.gz

cat barcode05/*.fastq.gz > Library1.barcode05.fastq.gz

cat barcode06/*.fastq.gz > Library1.barcode06.fastq.gz

cat barcode07/*.fastq.gz > Library1.barcode07.fastq.gz

cat barcode08/*.fastq.gz > Library1.barcode08.fastq.gz

cat barcode09/*.fastq.gz > Library1.barcode09.fastq.gz

cat barcode10/*.fastq.gz > Library1.barcode10.fastq.gz

cat barcode11/*.fastq.gz > Library1.barcode11.fastq.gz

# Merge fastq files in library 2

cat barcode02/*.fastq.gz > Library2.barcode02.fastq.gz

cat barcode03/*.fastq.gz > Library2.barcode03.fastq.gz

cat barcode05/*.fastq.gz > Library2.barcode05.fastq.gz

cat barcode06/*.fastq.gz > Library2.barcode06.fastq.gz

cat barcode07/*.fastq.gz > Library2.barcode07.fastq.gz

cat barcode08/*.fastq.gz > Library2.barcode08.fastq.gz

cat barcode09/*.fastq.gz > Library2.barcode09.fastq.gz

cat barcode10/*.fastq.gz > Library2.barcode10.fastq.gz

cat barcode11/*.fastq.gz > Library2.barcode11.fastq.gz

# Merge fastq files in library 3

cat barcode02/*.fastq.gz > Library3.barcode02.fastq.gz

cat barcode03/*.fastq.gz > Library3.barcode03.fastq.gz

cat barcode05/*.fastq.gz > Library3.barcode05.fastq.gz

cat barcode06/*.fastq.gz > Library3.barcode06.fastq.gz

cat barcode07/*.fastq.gz > Library3.barcode07.fastq.gz

cat barcode08/*.fastq.gz > Library3.barcode08.fastq.gz

cat barcode09/*.fastq.gz > Library3.barcode09.fastq.gz

cat barcode10/*.fastq.gz > Library3.barcode10.fastq.gz

cat barcode11/*.fastq.gz > Library3.barcode11.fastq.gz

exit