#!/bin/bash

python /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/pyscripts/reg_sldsc.py \
    --run_reg \
    --no_intersections \
    --methods OLS Lasso Ridge\
    --ld /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/annotations/EAS/baselineLD_v2.1/BL_combined_ld.csv \
    --sumstat $1 \
    --weight_ld /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/weights/EAS/weights.EAS.hm3_noMHC.csv \
    --leave_chrom $2 \
    --outfile_prefix /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/results/EAS/chr${2}/${3}_BLv2.1 \
    --annot_prefix /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/annotations/EAS/baselineLD_v2.1/BL. \
    --rsid /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/annot_rsid/EAS/rsid. \
    --m550_rsid /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/annot_5_50_rsid/EAS/rsid.
