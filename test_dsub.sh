#!/bin/bash

dsub \
    --provider google-v2 \
    --project encode-uk-biobank \
    --zones "us-east1-b" \
    --logging gs://regularized_sldsc/logging/ \
    --image "gcr.io/encode-uk-biobank/hs-ldsc:latest" \
    --script "run_test_reg.sh"
