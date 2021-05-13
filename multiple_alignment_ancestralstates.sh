#prep files for LASTZ (need nib and fasta sizes)
while read chrom
do 

mkdir -p /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/AT${chrom}
mkdir -p /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/mga/AT${chrom}


/ohta/apps/kentUtils/bin/faToNib /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/species/palmeri/${chrom}.fa /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/species/palmeri/${chrom}.nib
/ohta/apps/kentUtils/bin/faSize -detailed /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/species/palmeri/${chrom}.fa >/ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/species/palmeri/${chrom}.fa.sizes

done < scaf_list



for spp in hypochondriacus #palmeri #araTha araLyr braDis braRap carPap cucSat fraVes glyMax malDom orySat popTri ricCom solTub sorBic theCac vitVin zeaMay medTru lotJap glySoj
do 
echo "LASTZ on ${chrom} in ${spp}"

#perform multiple alignments, one chromosome at at time (to PALMERI)
cat scaf_list | parallel -j 16 "/gran1/felix.beaudry/apps/lastz-distrib-1.04.00/src/lastz /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/{}.fa /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/species/palmeri/palmeri.fa --step=10 --gapped --nochain --gfextend --strand=both --output=/ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/mga/AT{}/AT_palmeri.{}.axt --format=axt --ambiguous=iupac" #for palmeri

#perform multiple alignments, one chromosome at at time (to HYPOCHONDRIACUS)
cat scaf_list | parallel -j 16 "/gran1/felix.beaudry/apps/lastz-distrib-1.04.00/src/lastz /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/{}.fa /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/species/hypochondriacus/hypochondriacus.fa --step=10 --gapped --nochain --gfextend --strand=both --output=/ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/mga/hypo/ATvAH{}/AT_hypochondriacus.{}.axt --format=axt --ambiguous=iupac"  #for hypo


#chain, get best orthologous alignment, output as MAF
echo "chaining on ${chrom} in ${spp}"

cat scaf_list | parallel -j 16 "/ohta/apps/kentUtils/bin/axtChain -minScore=10000 -linearGap=loose /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/mga/AT{}/AT_palmeri.{}.axt /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/ /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/species/palmeri /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/mga/AT{}/AT_palmeri.{}.chain"

cat scaf_list | parallel -j 16 "/ohta/felix.beaudry/beyond/multiz/selectChains /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/mga/AT{}/AT_palmeri.{}.chain /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/mga/AT{}/AT_palmeri.{}.chain.ortho 9999"

cat scaf_list | parallel -j 16 "/ohta/apps/kentUtils/bin/chainToAxt /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/mga/AT{}/AT_palmeri.{}.chain.ortho /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/species/palmeri /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/mga/AT{}/AT_palmeri.{}.ortho.axt"

cat scaf_list | parallel -j 16 "/ohta/apps/kentUtils/bin/axtToMaf /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/mga/AT{}/AT_palmeri.{}.ortho.axt /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/AT_{}.fa.sizes /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/species/palmeri/palmeri.fa.sizes /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/mga/AT{}/AT_palmeri.{}.ortho.maf"


#sort MAF
cat scaf_list | parallel -j 16  "bash maf-sort.sh /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/mga/AT{}/AT_palmeri.{}.ortho.maf >  /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/mga/AT{}/AT_palmeri.{}.sort.maf"

#rename
cat scaf_list | parallel -j 16 "cat /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/mga/AT{}/AT_palmeri.{}.sort.maf | sed '/^a/{n;s/{}/tuber.{}/}' | sed '/^s/{n;s/Scaffold_/palm.Scaffold_/}' | gzip > /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/renamed/{}.sorted.renamed.maf.gz"

#cd /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/renamed/
#parallel gzip ::: *

#maf to bed
cat scaf_list | parallel -j 16 "python /ohta/julia.kreiner/software/WGAbed/maf_to_bed.py -i /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/renamed/{}.sorted.renamed.maf.gz -r tuber -c {} > /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/renamed/{}.bed"

#keep only unique positions where the allele is from the best scoring block
cat scaf_list | parallel -j 16 "cat /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/renamed/{}.bed | sort -t$'\t' -r -n -k10 | awk '!x[$2]++' | sort -k2n > /ohta/julia.kreiner/palmeri/multiple_alignments/felix_script/renamed/{}.bestmatch.bed"
#Idea is to sort the file in ascending order first (sort -n -k2) and reverse it (-r) on column 2 (which now will be descending order)

#convert to vcf, keeping only sites that differ from ref and are biallelic
sed 's/,/   /g' Scaffold_4.bestmatch.sorted.bed | awk '($4 != $3)'  | awk 'length($4) == 1' | awk '{print $1 "\t" $2 "\t." "\t" $3 "\t" $4 "\t30" "\t." "\t." "\tGT" "\t1/1"}' > ancestral_alleles_palm_SNPs_sorted.vcf

#FINALLY,
#modify reference genome to incorporate ancestral states
java -jar /ohta1/tyler.kent/Software/GATK/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker \
-V ancestral_alleles_palm_SNPs_sorted.vcf \
-R /ohta1/julia.kreiner/waterhemp/data/fixed_assembly/reveal_psuedoassembly/toshare/AMATA_finishedtohypo_renamedfinal.fasta \
-o AMATA_finishedtohypo_renamedfinal_palmacnestral_modheader.fasta
