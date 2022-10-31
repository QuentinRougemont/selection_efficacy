library(ggplot2)
library(dplyr)
library(data.table)
library(magrittr)

setwd("/home/quentin/postdoc/coho/wgs/2021/01.analysis/12pinpis/new/")
#region=read.table("wgs_lat_long.txt",T)  
# region <- read.table("strata.txt",T)
# region <- read.table("strata2.txt",T)
region <- read.table("strata3.txt",T)


pnps_diplo = read.table("DIPLOID/pnps.diploidonly.txt",T)
pnps_trtra = read.table("pnps.tetra_19.01.txt",T)
boot <- read.table("DIPLOID/06.BOOT/boot.pnps.txt",T)
boot$POP<- as.factor(boot$POP)
boot <- dplyr::select(boot, POP, pNpSratio)

boot <- merge(boot, region, by.x="POP", by.y="POP")

myColors <- c("blue","orange",
              "red","green",
              "darkviolet",
              "springgreen4")
names(myColors) <- levels(boot$region)
colScale <- scale_colour_manual(name = "POP",values = myColors)

overallmean = data.frame(mean(pnps_trtra$pNpSratio))
colnames(overallmean) = "mean"
vertical.lines = pnps_trtra$pNpSratio
p<- ggplot(boot, aes(x=pNpSratio)) +
    geom_histogram() +
   labs(x= expression(pi["n"]/pi["s"]))
p <- p + scale_color_brewer(palette = "Dark2") + theme_classic() +
  theme(legend.position="top")

p <- p + geom_vline(xintercept=vertical.lines, color="red",
             linetype="dashed")

p <- p + theme(axis.title.x=element_text(size=20, family="Helvetica",face="bold"),
               axis.text.x=element_text(size=20,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
               axis.title.y=element_text(size=20, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
               axis.text.y=element_text(size=14,family="Helvetica",face="bold"),
               strip.text.x = element_text(size=18))
#add personalized legend!
p

pdf(file="histogram_boot_19_01.pdf", 12,8)
p
dev.off()

a = read.table("pnps.all_pop.all_gene.txt",T)
a$psratio=a$Ps/a$NSS
a$pnratio=a$Pn/(a$Size-a$NSS)
a$pnps=a$pnratio/a$psratio

b=dplyr::select(a,POP, status, pnps)
b <- na.omit(b)
b <-b%>% filter(pnps <5)
ggplot(b, aes(x=status, y=pnps)) + geom_boxplot()

