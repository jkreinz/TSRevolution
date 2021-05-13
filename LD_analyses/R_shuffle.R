#!/usr/bin/env Rscript

#script to randomize assignment of TSR mutations across individuals in a VCF file
library(data.table)

#read file with 3 TSR mutations
df = fread("/plas1/george.sandler/amaranthus/non_syn_vcf/permute/Amaranthus_essex_TSR.vcf", header = F, sep = "\t")


df2 = t(df) #flip dataframe
df3 = df2[10:49,] #isolates fields with genotypes
df5 = df2[1:9,] #isolate information fields
df4 <- df3[sample(nrow(df3)),] #scramble assignment of genotypes
df6 = rbind(df5,df4) #re-combine permuted genotypes with site info fields
df7 = t(df6) #flip back into VCF format

write.table(df7, "/plas1/george.sandler/amaranthus/non_syn_vcf/permute/permute_gts.vcf", sep = "\t", col.names = F, quote = F, row.names = F)
