rm(list=ls())
library(qtl)
library(ggplot2)
library(dplyr)
library(gridExtra)
ril <- read.cross(format = "csv",".", file = "total_out.input.csv", genotypes = c("a","h","b"),alleles=c("a","b"))
summary(ril)

par(mfrow=c(4,2))
for(i in 1:8){plotPheno(ril, pheno.col=i)}
pairs(jitter(as.matrix(ril$pheno) ), cex=0.6, las=1)
par(mfrow=c(1,1))

ril <- jittermap(ril)
## calc.genoprob calculating genotype frequency; step(cM) determine the density of QTL
ril <- calc.genoprob(ril, step=1, error.prob=0.001)

plotMap(ril)

out.hk<-scanone(ril,pheno.col=1:8,method="hk")
par(mfrow = c(1,1))

out.cim.20 <- cim(ril, pheno.col=5,n.marcovar=10, window=20,method=c("hk"),error.prob=0.0001,map.function=c("kosambi"))
chr <- c(1, 2,3, 4,5, 6,7,8,9,10,11,12,13, 14)
#plot(out.hk, out.cim.20, chr=chr, ylab="LOD score",col=c("blue", "red"), main="window = 20 cM")
#plot(out.cim.20, chr=chr, ylab="LOD score",col=c("blue", "red"), main="window = 20 cM")
operm<-scanone(ril,method = "hk",n.perm = 100,verbose=TRUE)
lod_p <- summary(operm, alpha=c(0.20, 0.05))
sum <- summary(out.cim.20, perms=operm, alpha = 0.05,format = "tabByChr",ci.function = "bayesint", pvalues=TRUE)
#add.threshold(out.cim.20,alpha = 0.05,perms = operm)
#add.cim.covar(out.cim.20, chr=chr, col="yellow")
#lodint(out.cim.20, 7, 1.5) # LOD support intervals
#bayesint(out.cim.20, 3, 0.95) # Bayes credible intervals
mqtl<-makeqtl(ril,chr=sum$lod[,1],pos=sum$lod[,2],what="prob")#make QTL
#rqtl<-refineqtl(ril,pheno.col=1, qtl=mqtl,method="hk",model="normal") #Correction of QTL position
#adqtl<-addqtl(ril,qtl=mqtl,pheno.col=1,method="hk",model="normal") #detecing other QTLs
#intqtl<-addint(ril,pheno.col=1,qtl=mqtl,method="hk",model="normal",pvalues=T) #interaction between QTLs
fqtl<-fitqtl(ril,dropone=T,get.ests=T,model="normal",qtl=mqtl,method="hk") #estimating the pve of each QTL
final_result <- data.frame(sum,fqtl$result.drop,mqtl$pos)
write.csv(final_result,file = 'TKW.cim.csv',quote = F)

##########################################
don <- out.cim.20 %>% 
  group_by(chr) %>% 
  summarise(chr_len=max(pos)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(out.cim.20, ., by=c("chr"="chr")) %>%
  arrange(chr, pos) %>%
  mutate( BPcum=pos+tot)

axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

ggplot(don, aes(x = BPcum, y = lod)) +
  geom_line(size = 1, alpha = 0.8,aes(color=as.factor(chr))) +
  scale_color_manual(values = rep(c("blue", "red"), 14 )) +
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0))+
  geom_hline(yintercept = lod_p[2,1], color = 'black', size = 1.2, linetype = 'twodash')+ 
  theme_minimal() +
  labs(title = "", x = "GW", y = "LOD")+
  theme(axis.text.y = element_text(size = 14,color = 'black'),
        axis.text.x = element_text(size = 14,color = 'black'),
        axis.title.y.left = element_text(size = 14),
        legend.position = "none")

