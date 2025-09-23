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

ml load HISAT2/2.2.1-gompi-2021b

hisat2 --dta -p 128 -x UCSC_Mus_musculus.GRCm38 -1 Ctr3_1_EKRN230046783-1A_HG5C2DSX7_L4_1.fq.gz -2 Ctr3_1_EKRN230046783-1A_HG5C2DSX7_L4_2.fq.gz -S Ctr3_1.sam
hisat2 --dta -p 128 -x UCSC_Mus_musculus.GRCm38 -1 Ctr3_2_EKRN230046784-1A_HG5C2DSX7_L4_1.fq.gz -2 Ctr3_2_EKRN230046784-1A_HG5C2DSX7_L4_2.fq.gz -S Ctr3_2.sam
hisat2 --dta -p 128 -x UCSC_Mus_musculus.GRCm38 -1 Ctr3_3_EKRN230046785-1A_HG5C2DSX7_L4_1.fq.gz -2 Ctr3_3_EKRN230046785-1A_HG5C2DSX7_L4_2.fq.gz -S Ctr3_3.sam
hisat2 --dta -p 128 -x UCSC_Mus_musculus.GRCm38 -1 Ctr3_4_EKRN230046786-1A_HG5C2DSX7_L4_1.fq.gz -2 Ctr3_4_EKRN230046786-1A_HG5C2DSX7_L4_2.fq.gz -S Ctr3_4.sam
hisat2 --dta -p 128 -x UCSC_Mus_musculus.GRCm38 -1 Dnmt1_1_EKRN230046775-1A_HG5C2DSX7_L4_1.fq.gz -2 Dnmt1_1_EKRN230046775-1A_HG5C2DSX7_L4_2.fq.gz -S Dnmt1_1.sam
hisat2 --dta -p 128 -x UCSC_Mus_musculus.GRCm38 -1 Dnmt1_2_EKRN230046776-1A_HG5C2DSX7_L4_1.fq.gz -2 Dnmt1_2_EKRN230046776-1A_HG5C2DSX7_L4_2.fq.gz -S Dnmt1_2.sam
hisat2 --dta -p 128 -x UCSC_Mus_musculus.GRCm38 -1 Dnmt1_3_EKRN230046777-1A_HG5C2DSX7_L4_1.fq.gz -2 Dnmt1_3_EKRN230046777-1A_HG5C2DSX7_L4_2.fq.gz -S Dnmt1_3.sam
hisat2 --dta -p 128 -x UCSC_Mus_musculus.GRCm38 -1 Dnmt1_4_EKRN230046778-1A_HG5C2DSX7_L4_1.fq.gz -2 Dnmt1_4_EKRN230046778-1A_HG5C2DSX7_L4_2.fq.gz -S Dnmt1_4.sam
hisat2 --dta -p 128 -x UCSC_Mus_musculus.GRCm38 -1 Kmt2a_1_EKRN230046779-1A_HG5C2DSX7_L4_1.fq.gz -2 Kmt2a_1_EKRN230046779-1A_HG5C2DSX7_L4_2.fq.gz -S Kmt2a_1.sam
hisat2 --dta -p 128 -x UCSC_Mus_musculus.GRCm38 -1 Kmt2a_2_EKRN230046780-1A_HG5C2DSX7_L4_1.fq.gz -2 Kmt2a_2_EKRN230046780-1A_HG5C2DSX7_L4_2.fq.gz -S Kmt2a_2.sam
hisat2 --dta -p 128 -x UCSC_Mus_musculus.GRCm38 -1 Kmt2a_3_EKRN230046781-1A_HG5C2DSX7_L4_1.fq.gz -2 Kmt2a_3_EKRN230046781-1A_HG5C2DSX7_L4_2.fq.gz -S Kmt2a_3.sam
hisat2 --dta -p 128 -x UCSC_Mus_musculus.GRCm38 -1 Kmt2a_4_EKRN230046782-1A_HG5C2DSX7_L4_1.fq.gz -2 Kmt2a_4_EKRN230046782-1A_HG5C2DSX7_L4_2.fq.gz -S Kmt2a_4.sam
exit
