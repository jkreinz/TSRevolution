#phase with SHAPEIT - pop and read back.
########################################

#split by scafs
while read scafs
do
vcftools --vcf pseudo_qual_dust_indel_fixedmissing_nodups.vcf --recode --recode-INFO-all --max-missing .9 --out ${scafs}_missing10p --chr ${scafs} & 
done < scaf_list

mv /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/filtering/Scaffold_*_missing10p.recode.vcf /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing/

#copy map files
cp /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/ldhat/LDhat_workflow/amaranth/*_finalformat_recombrate.txt /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing

#get scaf list
cp /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/filtering/scaf_list /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/allchroms_phasing

#get SNP positions
while read scafs
do
grep -v "#" ${scafs}_missing10p.recode.vcf | awk '{print $2}' > ${scafs}.pos
done < scaf_list

#check for duplicates in map file
while read scafs
do
uniq ${scafs}_finalformat_recombrate.txt > ${scafs}_finalformat_recombrate_uniq.txt
done < scaf_list

#impute recombination rates between SNPs in our VCF
#arg1 = map file, arg2 = listofpos, arg3 = output name
while read scafs
do
Rscript imputemap.R ${scafs}_finalformat_recombrate_uniq.txt ${scafs}.pos ${scafs}_shapeit.map 
done < scaf_list

#remove unmatched position in vcf
while read scafs
do
#awk -v var=${scafs} '{print var "\t" $1 }' ${scafs}_shapeit.map > ${scafs}_map.pos
awk '{if ($2 <= 0) print $1 " " 0.001 " " $3; else print $1 " " $2 " " $3}' ${scafs}_shapeit.map | sed '1 i\pposition rrate gposition' > ${scafs}_shapeit2.map
#vcftools --vcf ${scafs}_missing10p.recode.vcf --recode --recode-INFO-all --positions ${scafs}_map.pos --out ${scafs}_missing_mapmatched &
done < scaf_list

##############################################
#OK, try read aware, population-level phasing
##############################################

#organize bam files and list
#run if phasing multiple scafs at once
ls *dd.bam > bamlist
while read scaf
do
paste <( awk -F'/' '{print $8}' bamlist | sed 's/.dd.bam//g') bamlist |  awk -v scaf="$scaf" '{print $0 "\t" scaf}' > final_bamlist_$scaf
done < scaf_list


#for all scafs if wanted (need to prep map files for all chromosomes)
while read scaf
do
extractPIRs --bam final_bamlist_$scaf \
            --vcf ${scaf}_missing10p.recode.vcf \
            --out /ohta2/julia.kreiner/commongarden_result/shapeit_phasing/myPIRsList_$scaf &
done < scaf_list

while read scaf
do
shapeit -assemble \
        --input-vcf /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/phasing/${scaf}_missing10p.recode.vcf \
        --input-pir /ohta2/julia.kreiner/commongarden_result/shapeit_phasing/myPIRsList_$scaf \
        -O /ohta/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/TSR_work/phasing/byhap/relate/readback_rerun/${scaf}_readback_haplotypes \
        --effective-size 500000 --thread 20

done < scaf_list

