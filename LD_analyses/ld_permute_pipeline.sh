#! /bin/bash


cd /plas1/george.sandler/amaranthus/non_syn_vcf/permute/julia_windows/653/
for i in {1..1000}
do
"/plas1/george.sandler/amaranthus/non_syn_vcf/permute/R_shuffle.R" #script to permute assignment of TSR mutations in a VCF

  for ind in 10 11 12 13 16 19 2 20 21 22 23 24 3 6 8 #list of windows to test for significant LD
  do
  
  #combine VCF of region of interest with VCF file of permuted TSR mutations
  cat /plas1/george.sandler/amaranthus/non_syn_vcf/permute/julia_windows/653/Essex_win_${ind}.reg.vcf "/plas1/george.sandler/amaranthus/non_syn_vcf/permute/permute_gts.vcf" > /plas1/george.sandler/amaranthus/non_syn_vcf/permute/julia_windows/Temp_permute_gts_${ind}.vcf

  #plink script to caluclate LD between scarmbled TSR mutations, and SNPs in region of interest
  "/plas1/george.sandler/apps/plink" --make-bed  --vcf /plas1/george.sandler/amaranthus/non_syn_vcf/permute/julia_windows/Temp_permute_gts_${ind}.vcf --out /plas1/george.sandler/amaranthus/non_syn_vcf/permute/plink/amar_plink --allow-extra-chr --double-id

  "/plas1/george.sandler/apps/plink" --bfile /plas1/george.sandler/amaranthus/non_syn_vcf/permute/plink/amar_plink --recode --tab --out /plas1/george.sandler/amaranthus/non_syn_vcf/permute/plink/amar_plink  --allow-extra-chr 

  #TSR mutation must be labelled with ID matched to --ld--snp option
  "/plas1/george.sandler/apps/plink" --file /plas1/george.sandler/amaranthus/non_syn_vcf/permute/plink/amar_plink --r --out /plas1/george.sandler/amaranthus/non_syn_vcf/permute/plink/amar_plink_permuted --allow-extra-chr  --with-freqs --inter-chr --ld-snp m2 



  #calculate null LD statistics in files created by Plink
  "/plas1/george.sandler/amaranthus/non_syn_vcf/permute/R_LD_calc.R" ${ind}
  done
done



