#!/usr/bin/env bash
source activate bcftools
input=pop_call_raw.vcf
x=2018

# raw combined vcf as input
vcftools --vcf ${input} --recode --recode-INFO-all --out ${x} 
echo recoded raw vcf
mv ${x}.recode.vcf ${x}.vcf   #31,687 records, 26,128 SNPs

vcffilter -f "AC > 0" ${x}.vcf > ${x}_AC.vcf #31,687 records, 26,128 SNPs
vcffilter -f "QUAL > 30" ${x}_AC.vcf > ${x}_AC_Q30.vcf  #30,486 records, 25,185 SNPs
vcffilter -f "QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ${x}_AC_Q30.vcf > ${x}_AC_Q30_SAF.vcf #28,166 records, 23,310 SNPs
echo filtered fields

vcftools --vcf ${x}_AC_Q30_SAF.vcf --minDP 10 --recode --recode-INFO-all --out ${x}_AC_Q30_SAF_DP.vcf #28,166 records, 23,310 SNPs
echo filtered based on depth

vcftools --vcf ${x}_AC_Q30_SAF_DP.vcf.recode.vcf --remove-indels --recode --recode-INFO-all --out ${x}_AC_Q30_SAF_DP_SNPs.vcf  
mv ${x}_AC_Q30_SAF_DP_SNPs.vcf.recode.vcf ${x}_AC_Q30_SAF_DP_SNPs.vcf #21,166 SNPs, 117 multiallelic sites
echo removed indels

#rm prophage regions saved in phage.bed
vcftools --vcf ${x}_AC_Q30_SAF_DP_SNPs.vcf --exclude-bed bed/phage.bed  --recode --recode-INFO-all --out ${x}_AC_Q30_SAF_DP_SNPs_ph
mv ${x}_AC_Q30_SAF_DP_SNPs_ph.recode.vcf ${x}_AC_Q30_SAF_DP_SNPs_ph.vcf #20,377 SNPs, 111 multiallelic sites
echo removed SNPs in prohage regions

#rm tandem repeat regions
vcftools --vcf ${x}_AC_Q30_SAF_DP_SNPs_ph.vcf --exclude-bed bed/repeat_n20.bed  --recode --recode-INFO-all --out ${x}_AC_Q30_SAF_DP_SNPs_ph_rp
mv ${x}_AC_Q30_SAF_DP_SNPs_ph_rp.recode.vcf ${x}_AC_Q30_SAF_DP_SNPs_ph_rp.vcf #20,359 SNPs
echo removed SNPs in tandem repeat regions

#remove SNPs with any missing data
vcftools --vcf ${x}_AC_Q30_SAF_DP_SNPs_ph_rp.vcf  --max-missing 1 --recode --recode-INFO-all --out ${x}_final_filt.vcf #17,437 SNPs 69 multiallelic sites
mv ${x}_final_filt.vcf.recode.vcf ${x}_final_filt.vcf

#bgzip and index vcf file for further analysis
bgzip -c 2018_final_filt.vcf > 2018_final_filt.vcf.gz
source activate gubbins #need samtools for tabix
tabix 2018_final_filt.vcf.gz 
