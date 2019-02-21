#!/bin/bash

python /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/pyscripts/reg_sldsc.py \
    --run_reg \
    --no_intersections \
    --methods OLS \
    --ld /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/annotations/one/1000G.EUR.QC.combined_ld.csv \
    --sumstat $1 \
    --weight_ld /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/weights/weights.hm3_noMHC.csv \
    --leave_chrom $2 \
    --outfile_prefix /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/results/Feb2019_eval_wld/chr${2}/${3}_one \
    --annot_prefix /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/annotations/one/one. \
    --rsid /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/annot_rsid/rsid. \
    --m550_rsid /n/groups/price/ran/high-dim-sldsc/real_pheno_pipeline/annot_5_50_rsid/rsid.
#    --reg_weights $5
