#!/bin/bash

c=$1
aria2c -x 16 -s 16 https://genetics.bwh.harvard.edu/downloads/Vova/Roulette/${c}_rate_v5.2_TFBS_correction_all.vcf.bgz
Rscript Rscripts/parse_muts.R ${c}
rm ${c}_rate_v5.2_TFBS_correction_all.vcf.bgz
