#!/bin/bash

python ~/regularized_sldsc/pyscripts_new/reg_sldsc.py \
    --run_reg \
    --no_intersections \
    --${2} \
    --ld ~/data/annotations/690combined_ld.csv \
    --sumstat $1 \
    --weight_ld ~/data/weights/weights.hm3_noMHC.csv \
    --leave_chrom $3 \
    --outfile_prefix ~/regularized_sldsc/results/real_pheno/new/chr${3}/${4}_${2}_690 \
    --annot_prefix /home/rancui/data/annotations/690annots/690annots. \
    --rsid ~/data/annotations/annot_rsid/rsid. \
    --m550_rsid ~/data/annotations/annot_5_50_rsid/rsid.
#    --reg_weights $5
