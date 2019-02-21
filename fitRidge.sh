#!/bin/bash

python /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/pyscripts/reg_sldsc.py \
    --run_reg \
    --no_intersections \
    --methods Ridge \
    --ld /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/annotations/allv2.1/all_combined_ld.csv \
    --sumstat $1 \
    --weight_ld /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/weights/weights.hm3_noMHC.csv \
    --leave_chrom ${2} \
    --outfile_prefix /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/results/chr${2}/${3}_allv2.1 \
    --annot_prefix /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/annotations/allv2.1/all. \
    --rsid /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/annot_rsid/rsid. \
    --m550_rsid /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/annot_5_50_rsid/rsid.
#    --reg_weights $5
