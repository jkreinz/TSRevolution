#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(data.table)

df = fread("/plas1/george.sandler/amaranthus/non_syn_vcf/permute/plink/amar_plink_permuted.ld")
df = na.omit(df)
md = median(df$R)
mn = mean(df$R)
vr = var(df$R)


line=(noquote(paste(md, mn, vr, args[1])))
write(line,file="/plas1/george.sandler/amaranthus/non_syn_vcf/permute/julia_windows/653_permutations.txt",append=TRUE)