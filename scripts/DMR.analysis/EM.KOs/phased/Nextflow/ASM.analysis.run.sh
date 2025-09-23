#!/bin/bash

ml use /hpcapps/lib-bio/modules/all/
ml load Nextflow/23.10.0

nextflow run \
ONT_methPhase.nf \
-c methPhase.config \
-profile elja_slurm \
--sample_ID EM_KO \
--modBAM_path EM_KO.bam \
--genome1_name B6 \
--genome2_name FVBNJ \
--genome1_path C57BL6J_b38.fa \
--genome2_path FVBNJ_b38_f.fa \
--anchorPath_g1 B6.anchor_info.csv \
--anchorPath_g2 FVBNJ.anchor_info.csv \
--vcf_g1ref FVBNJ-alt_B6-ref.no_indels.phased.vcf.gz \
--vcf_g2ref MUMmer_variants/B6-alt_FVBNJ-ref.no_indels.phased.vcf.gz \
--pairwise_alignment path/to/pairwiseAlignment \
--out_dir path/to/output/directory \
-with-trace

nextflow run \
ONT_methPhase.nf \
-c methPhase.config \
-profile elja_slurm \
--sample_ID Control \
--modBAM_path Control.bam \
--genome1_name B6 \
--genome2_name FVBNJ \
--genome1_path C57BL6J_b38.fa \
--genome2_path FVBNJ_b38_f.fa \
--anchorPath_g1 B6.anchor_info.csv \
--anchorPath_g2 FVBNJ.anchor_info.csv \
--vcf_g1ref FVBNJ-alt_B6-ref.no_indels.phased.vcf.gz \
--vcf_g2ref B6-alt_FVBNJ-ref.no_indels.phased.vcf.gz \
--pairwise_alignment pairwiseAlignment \
--out_dir path/to/output/directory \
-with-trace