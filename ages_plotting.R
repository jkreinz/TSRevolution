<-read.table("~/TSRages_mut11_300plus.txt",header = T)
#ages<-ages[with(ages, order(V3, V10)), ]
#ages$mut<-c(1:nrow(ages))
#write.table(ages,"TSRages.txt",col.names=T,quote=F)

library(ggplot2)
library(scales)
library(PNWColors)


names(ages)[4]<-"TSR locus"
ages$`TSR locus`<-as.factor(ages$`TSR locus`)
star<-pnw_palette(name="Starfish",n=7,type="discrete")

scaleFUN <- function(x) sprintf("%.3f", x)
head(ages)

ggplot(data=ages, aes(mut,allele_age_all/xc,color=`TSR locus`)) +
  geom_point() +
  geom_errorbar(data=ages,aes(ymin=age_05/166.5894, ymax=age95/166.5894, width=.5)) +
  scale_y_continuous(trans = log2_trans(),labels=scaleFUN) +
  annotation_logticks(sides = 'l') +
  theme_bw() +
  ylab("Allele age (Years)") +
  xlab("Mutational Origin") +
  scale_color_manual(values = c(star[4],star[1],star[7])) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11))

ages$denovo<-c(9.81E-06,
          0.000885252,
          0.005650493,
          0.1401669,
          0.6623966,
          0.6397279,
          5.58E-06,
          0.7168271,
          0.2841429,
          0.5886909,
          0.5886909)


ages$sgv<-c(0.2449,
       NA,
       0.9800,
       0.0013,
       0.3289,
       0.3668,
       2.67E-14,
       0.0765,
       0.5470,
       0.2703,
       0.7124)

ps<-ages[,c(1,4,17,18)]
ps_long<-melt(ps,id.vars=c("mut","TSR locus"))
names(ps_long)
names(ps_long$variable)<-"Model of Selection"

ggplot(data=ps_long, aes(mut,-log10(value),color=`TSR locus`,shape=variable)) +
  geom_point(cex=3) +
  scale_y_continuous(trans = log2_trans(),labels=scaleFUN) +
  annotation_logticks(sides = 'l') +
  theme_bw() +
  ylab("-log10(p-value)") +
  xlab("Mutational Origin") +
  scale_color_manual(values = c(star[4],star[1],star[7])) +
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9,10,11)) +
  geom_hline(yintercept = -log10(0.05),lty="dashed")


