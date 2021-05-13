#LINKAGE DISEQUILIBRIUM
for i in mut1 mut2 mut3
do
plink --file amar_plink --r --out Amaranthus_essex_mis_$i --allow-extra-chr --with-freqs --ld-snp $i --chr Scaffold_11 --keep essex
plink --file amar_plink --r --out Amaranthus_walpole_mis_$i --allow-extra-chr --with-freqs --ld-snp $i --chr Scaffold_11 --keep walpole
plink --file amar_plink --r --out Amaranthus_midwest_mis_$i --allow-extra-chr --with-freqs --ld-snp $i --chr Scaffold_11 --keep midwest
done

#DIVERSITY
vcftools --vcf scaf5_haploid_parallel_wh.vcf --keep essexS --out essex_sus_haploid --site-pi #diversity across susceptible haplotypes
vcftools --vcf scaf5_haploid_parallel_wh.vcf --keep essexR --out essex_res_haploid --site-pi #diversity across resistant haplotypes

vcftools --vcf scaf5_haploid_parallel_wh.vcf --keep walpoleS --out walpole_sus_haploid --site-pi #diversity across susceptible haplotypes
vcftools --vcf scaf5_haploid_parallel_wh.vcf --keep walpoleR --out walpole_res_haploid --site-pi #diversity across resistant haplotypes

vcftools --vcf scaf5_haploid_parallel_wh.vcf --keep midwestS --out midwest_sus_haploid --site-pi #diversity across susceptible haplotypes
vcftools --vcf scaf5_haploid_parallel_wh.vcf --keep midwestR --out midwest_res_haploid --site-pi #diversity across resistant haplotypes

#XPEHH
selscan --xpehh --vcf-ref ../invariant_vcf/midwestS_haploid_var_nomissing_pseudophased_posfiltmap.recode.vcf --vcf ../invariant_vcf/midwestR_haploid_var_nomissing_pseudophased_posfiltmap.recode.vcf --map ../invariant_vcf/phased_midwest.map --out XPEHH_midwest_haploid --threads 10

selscan --xpehh --vcf-ref ../invariant_vcf/walpoleS_haploid_var_nomissing_pseudophased_posfiltmap.recode.vcf --vcf ../invariant_vcf/walpoleR_haploid_var_nomissing_pseudophased_posfiltmap.recode.vcf --map ../invariant_vcf/phased_walpole.map --out XPEHH_walpole_haploid --threads 10

selscan --xpehh --vcf-ref ../invariant_vcf/essexS_haploid_var_nomissing_pseudophased_posfiltmap.recode.vcf --vcf ../invariant_vcf/essexR_haploid_var_nomissing_pseudophased_posfiltmap.recode.vcf --map ../invariant_vcf/phased_essex.map --out XPEHH_essex_haploid --threads 10
