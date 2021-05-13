library(phytools)
library(dispRity)

#import ARG tree
tree<-read.tree("20Ksnps_c1_t50_n500.tree",skip=3)
tree<-tree[3][[1]] #1 == PPO, 2 == 653, 3 == 574

tree<-read.tree("20Ksnps_c1_t50_n500_440_4origins_574.tree")
#tree<-drop.tip(tree,tip = c("7_12-1","15_10-2","E10-2","G2-2")) #these are branches unresolved to the base of the tree

#order data by tips in tree
tips<-as.data.frame(tree$tip.label)
inds<-read.table("~/Documents/all653.txt")
names(inds)<-c("Samples","ALS_653")
row.names(inds)<-inds$Samples #set row names to samples
ordered<-inds[match(tree$tip.label, inds$Samples),] #exlude inds not found in tree

res_allele="A" #A for ALS574 and PPO,T for ALS653
sus_allele="C" #C for ALS574 and ALS653, G for PPO

#find runs of resisistant alleles with clusters containing more than 3 haplotypes
runs<-rle(ordered$ALS_653)
myruns = which(runs$values == res_allele & runs$lengths >=2)

#get end position of runs
runs.lengths.cumsum = cumsum(runs$lengths)
ends = runs.lengths.cumsum[myruns]

#get start position of runs
newindex = ifelse(myruns>1, myruns-1, 0)
starts = runs.lengths.cumsum[newindex] + 1
if (0 %in% newindex) starts = c(1,starts)

#starts<-starts[c(1,2)]
#ends<-ends[c(1,3)]

mut<-list()
#extract runs/haplotypes corresponding with diff independent origins
for (i in 1:length(myruns)) {
  mut[[i]] <- print(ordered$Samples[starts[i]:ends[i]])
}

#for ALS574
mut[[6]]<-c("D11-1","D1-1","D5-1","D12-2","D10-2","D6-2","D2-1","D2-2")
mut[[2]]<-mut[[2]][-c(40:47)]
mut[[2]]<-mut[[2]][-c(40:44)]
mut[[4]]<-c(mut[[4]],mut[[5]][1])
mut[[5]]<-mut[[5]][-c(1)]

#for rerun of ALS653
mut[[2]]<-c(mut[[2]],mut[[3]])
mut[[3]]<-NULL

#get branch corresponding to common ancestor
#count number of lineages that exist at time of common ancestory

##################
#get proportions of res versus susceptible lineages in the current day and at each origin
#################
library(data.table)
res_node<-list(); mut_age<-list(); nt0_sus<-list()
for (i in 1:length(mut)) {
  keeptips<-mut
  keeptips[[i]] <- NULL
  justoneorigin<-drop.tip(tree,tip = unlist(keeptips))
  res_node[[i]]<-findMRCA(justoneorigin,tips=mut[[i]])
  mut_age[[i]]<-tree.age(justoneorigin)[res_node[[i]],]$ages
  n_lineages<-as.data.frame(ltt.plot.coords(justoneorigin,backward = T,tol = 1))
  nt0_sus[[i]]<-last(n_lineages[n_lineages$time <= -mut_age[[i]],]$N)
  
}

#test
pruned.tree<-drop.tip(tree, tree$tip.label[-na.omit(match(mut[[3]], tree$tip.label))])
plot(pruned.tree)

#set variables
nt0_sus 
nt0_res<-as.list(rep(1,length(nt0_sus)))

########################
#coal approach (e.g. speidel et al., 2019; RELATE)
#calculating p-values for the skew in the number of offspring of resistant lineages since their origin (TMRCA)
########################

ks<-unlist(nt0_sus) #number of susceptible lineages when mutation arose + 1
fns<-lengths(mut) #number of current carriers

n_sus <-nrow(inds) - (sum(fns)) 
Ns<-list() 
for(i in 1:length(mut)) {
   Ns[[i]] <- n_sus + fns[i]
}
Ns<-unlist(Ns) #current day sample size, excluding irrelevant resistant clusters

upto <- Ns - ks + 2 #mutation spreads to at least fn haplotypes, and at max, Ns - ks + 2
range_of_spread <- upto - fns

num_list<-list(); denom_list<-list(); pfn_list<-list()
for (i in 1:length(mut)){ #for every independent origin
  num<-list(); denom<-list(); pfn<-list()

    for (k in 1:range_of_spread[i] ) { #for fn to Ns - ks + 2; 
      #calculate the prob that our allele reaches at least our present day freq
    fns2<-fns[i] + (k-1)
    num[[k]] <- (fns2 - 1)*(choose(n = Ns[i] - fns2 -1, k = ks[i] - 3))
    denom[[k]] <- choose(n = Ns[i] - 1, k = ks[i] - 1)
    pfn[[k]] <- num[[k]]/denom[[k]]
    }
  
  num_list[[i]] <- num
  denom_list[[i]] <- denom
  pfn_list[[i]] <- pfn
}


sum(unlist(pfn_list[[1]])) #pvalue for origin 1
sum(unlist(pfn_list[[2]])) #pvalue for origin 2
sum(unlist(pfn_list[[3]])) #pvalue for origin 3
sum(unlist(pfn_list[[4]])) #pvalue for origin 4
sum(unlist(pfn_list[[5]])) #pvalue for origin 5
sum(unlist(pfn_list[[6]])) #pvalue for origin 6

#

#################
#modify to test for recent selection
#calculating p-values for the skew in the number of offspring of resistant lineages over the last 0.2% of the tree

##################
res_node<-list(); mut_age<-list(); nt0_sus<-list();nt0_res<-list()
for (i in 1:length(mut)) {
  pruned.tree<-drop.tip(tree, tree$tip.label[-na.omit(match(mut[[i]], tree$tip.label))])
  n_lineages<-as.data.frame(ltt.plot.coords(pruned.tree,backward = T,tol = 1))
  nt0_res[[i]]<-last(n_lineages[n_lineages$time <= -200,]$N) #2% of the tree
  
  sus_samps<-ordered$Samples[ordered$ALS_653 == sus_allele] #exlude inds not found in tree
  pruned.tree.sus<-drop.tip(tree, tree$tip.label[-na.omit(match(sus_samps, tree$tip.label))])
  n_lineages<-as.data.frame(ltt.plot.coords(pruned.tree.sus,backward = T,tol = 1))
  nt0_sus[[i]]<-last(n_lineages[n_lineages$time <= -200,]$N) #2% of the tree
  
}

nt0_res[[6]] <- 0 #for ALS574 origin 2 since predates our time range (note origins not in the same order as depicted in manuscript)
nt0_res
nt0_sus


## coalescent calcs ##
ks<-unlist(nt0_sus) #number of susceptible lineages when mutation arose + 1
fns<-lengths(mut) #number of current carriers

n_sus <-nrow(inds) - (sum(fns)) 
Ns<-list() 
for(i in 1:length(mut)) {
  Ns[[i]] <- n_sus + fns[i]
}
Ns<-unlist(Ns) #current day sample size, excluding irrelevant resistant clusters

upto <- Ns - ks + 2 #mutation spreads to at least fn haplotypes, and at max, Ns - ks + 2
range_of_spread <- upto - fns


num_list<-list(); denom_list<-list(); pfn_list<-list()
for (i in 1:length(mut)){ #for every independent origin
  num<-list(); denom<-list(); pfn<-list()
  
  for (j in 1:range_of_spread[i] ) { #for fn to Ns - ks + 2; 
    #calculate the prob that our allele reaches at least our present day freq
    fns2<-fns[i] + (j-1)
    num[[j]] <- (choose(n=fns2 - 1,k=(unlist(nt0_res[[i]] - 1))))*(choose(n = Ns[i] - fns2 -1, k = ks[i] - 1))
    denom[[j]] <- choose(n = Ns[i] - 1, k = unlist(nt0_res[[i]]) + ks[i] - 1)
    pfn[[j]] <- num[[j]]/denom[[j]]
  }
  
  num_list[[i]] <- num
  denom_list[[i]] <- denom
  pfn_list[[i]] <- pfn
}
pfn_list

#these are the resulting p values from 6 origins of ALS574, 2 of ALS653, and 3 of PPO
denovo<-c(9.81E-06,
0.000885252,
0.005650493,
0.1401669,
0.6623966,
0.6397279,
5.58E-06,
0.7168271,
0.5886909,
0.2841429,
0.5886909)
denovo_adj<-p.adjust(denovo,method="fdr")
sum(denovo_adj < 0.05)

sgv<-c(5.54E-12,
2.71E-07,
0.1401669,
0.003807775,
0.07407657,
1.81E-18,
7.18E-04,
0.2703239,
0.000178024,
0.2612265)

sgv_adj<-p.adjust(sgv,method="fdr")
sum(sgv < 0.05)
sum(sgv_adj < 0.05)



