#subset region, 20K SNPs
begin=282560 #starting snp=23872242
end=283560 #ending snp=24801408
sed -n " ${begin},${end}p " ../scaf11_wTSR_flipped.haps > 20Ksnps_wTSR_WH.haps


awk 'BEGIN {OFS = "\t"}
{       printf $3"\t";
        for (i=6; i<=NF; ++i)
                if (toupper($i) == "0"){
                        printf $i = $4 ;
                      }
                else {
                        printf $i = $5 ;
                      }
                printf "\n";
}' 20Ksnps_wTSR_WH.haps > 20Ksnps_wTSR_WH.sites

#run argweaver on phased --> sites file
arg-sample -s 20Ksnps_wTSR_WH.sites -N 500000 -r 7e-9 -m 1.8e-8 --ntimes 20 --maxtime 100e3 -c 1 -n 500 -o reruns/20Ksnps_c1_t20_n500 --overwrite

#convert smc file to bed
for i in {1..500}
do
smc2bed reruns/20Ksnps_c1_t20_n500.${i}.smc.gz > reruns/20Ksnps_c1_t20_n500.${i}.bed.gz
done

#merge MCMC samples
cat 20Ksnps_c1_t50_n500.*.bed | sort -nk2,2 -nk4,4 | bgzip > 20Ksnps_c1_t50_n500.allsamples.sorted.bed.gz
tabix 20Ksnps_c1_t50_n500.allsamples.sorted.bed.gz

#trees for ALS573
arg-summarize -a 20Ksnps_c1_t50_n500.allsamples.sorted.bed.gz --tree --region chr11:24250050-24250051
#trees for ALS653
arg-summarize -a 20Ksnps_c1_t50_n500.allsamples.sorted.bed.gz --tree --region chr11:24249813-24249814
#trees for PPO
arg-summarize -a 20Ksnps_c1_t50_n500.allsamples.sorted.bed.gz --region chr11:23999755-23999756 --time-file 50times2.txt --tree

#summary stats from ARG
#first keep just last 200 samples
zcat 20Ksnps_c1_t50_n500.allsamples.sorted.bed.gz | awk '$4 > 290' | bgzip > 20Ksnps_c1_t50_n500.after300.sorted.bed.gz


#REDOING
#make snps file
cat <(head -n+1 Scaffold_11_missing10p_w653_PPOrealdel_whatshap_shapeit4_WH.sites) <(grep "23999751" Scaffold_11_missing10p_w653_PPOrealdel_whatshap_shapeit4_WH.sites | sed 's/^/chr11\t23999750\t/g') <(grep "24249814" Scaffold_11_missing10p_w653_PPOrealdel_whatshap_shapeit4_WH.sites | sed 's/^/chr11\t24249813\t/g')  <(grep "24250051" Scaffold_11_missing10p_w653_PPOrealdel_whatshap_shapeit4_WH.sites | sed 's/^/chr11\t24250050\t/g') > finalsnps_forstats.bed 

for i in {1..1250}
do
smc2bed 20Ksnps_whatshap_shapeit4_c1_t30_r8.${i}.smc.gz > 20Ksnps_whatshap_shapeit4_c1_t30_r8.${i}.bed.gz
done

#compile all MCMC samples
rm 20Ksnps_whatshap_shapeit4_c1_t30_r8.allsamples.bed
for i in `seq 0 10 1260`
do 
awk -v var="$i" '{print $1 "\t" $2 "\t" $3 "\t" var "\t" $5}'  20Ksnps_whatshap_shapeit4_c1_t30_r8.${i}.bed.gz >> 20Ksnps_whatshap_shapeit4_c1_t30_r8.allsamples.bed
done

sort -nk2,2 -nk4,4 -T ./ 20Ksnps_whatshap_shapeit4_c1_t30_r8.allsamples.bed | bgzip > 20Ksnps_whatshap_shapeit4_c1_t30_r8.allsamples.sorted.bed.gz

#make sure that origin subsets also contain susceptible individuals, or else will error out
#get allele ages by origin by mutation
for i in ALS574_1sus ALS574_2sus ALS574_3sus ALS653_1sus ALS653_2sus PPO_1sus PPO_2sus
do 
arg-summarize -a argweaver/20Ksnps_whatshap_shapeit4_c1_t30_r8.allsamples.sorted.bed.gz -A -T --snp-file finalsnps_forstats.bed.gz --subset vis_haps/${i} -M -Q 0.05,0.95 > vis_haps/${i}.alleleages
done

#arg-summarize -a 20Ksnps_c1_t50_n500.allsamples.sorted.bed.g --snp-file finalsnps_forstats.bed.gz -A -T -Q 0.05,0.95 #allelic age and time to the most recent common ancestor for each locus

