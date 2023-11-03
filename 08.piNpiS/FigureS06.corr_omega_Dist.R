
library(dplyr)
library(ggplot2)
library(ggpubr)

a = read.table("stats_table" , h = T, sep = "\t")


summary(lm(a$dist_to_ocean ~a$allsites_3sp_omegaNA))
summary(lm(a$dist_to_ocean ~a$allsites_3sp_omegaA))

th_plot <-     theme(axis.title.x=element_text(size=14, family="Helvetica",face="bold"),
                     axis.text.x=element_text(size=14,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
                     axis.title.y=element_text(size=18, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
                     axis.text.y=element_text(size=14,family="Helvetica",face="bold"),
                     strip.text.x = element_text(size=18),
                     panel.grid.major = element_blank()) 



ggplot(a, aes(x = dist_to_ocean, y = allsites_3sp_omegaNA) )+
       geom_point() +
       geom_smooth(method = "lm") +
       th_plot + 
       stat_cor(
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
           label.x = 3, label.y = 0.30) + th_plot