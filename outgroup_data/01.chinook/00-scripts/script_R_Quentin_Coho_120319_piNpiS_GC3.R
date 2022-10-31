# TL - 280318
library("ggplot2")

pnps=read.table("allpops_piNpiS",h=FALSE)

############# Plot pNpS ~ pS ################
ggplot(pnps,aes(x = V6, y =V8))+
  geom_smooth(method="lm", level=0.95, col="darkgrey",fill="lightgrey")+ geom_point(cex=2) + 
  xlab(expression(pi[S]))+ ylab(expression(pi[N]/pi[S]))+
  #scale_linetype_manual(values=c("solid", "dashed"))+
  theme_bw()+
  theme(legend.key.size = unit(2.5,"line"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=30,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=30,angle=90,hjust=.5,vjust=.5,face="italic"))

summary(lm(formula=pnps$V8 ~ pnps$V6)) # NS p~0.2


##### pNpS ~ median GC3 all pops ###
#pnpsgc3=read.table("piNpiS_windowsGC3_medianGC3_allpops",h=TRUE)
pnpsgc3=read.table("all_pop.pNpS.1mb.txt",h=TRUE)

ggplot(pnpsgc3,aes(x = medianGC3, y =pNpS,colour=POP))+
  geom_point(cex=2) + geom_smooth(method="lm", fill=NA)+
  #geom_segment(aes(x = MeanFst, y = cMperMb,
  #                 xend = MeanFst, yend = Fitted, colour=pair))+
  xlab("GC3")+ ylab(expression(pi[N]/pi[S]))+
  #scale_colour_manual(values=group.colors)+
  #scale_linetype_manual(values=c("solid", "dashed"))+
  theme_bw()+
  theme(legend.key.size = unit(2.5,"line"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=30,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=30,angle=90,hjust=.5,vjust=.5,face="italic"))

###### pNpS ~ GC3 per pop#########
library(dplyr)

w <- pnpsgc3 %>%
  group_by(POP) %>% # You can add here additional grouping variables if your real data set enables it
  do(mod = lm(pNpS ~ medianGC3, data = .)) %>%
  mutate(Slope = summary(mod)$coeff[2], 
         pvalue = summary(mod)$coeff[8],
         R2adj = summary(mod)$adj.r.squared) %>%
  select(-mod)

str <-read.table("strata3.txt", T)
w <- merge(w, str)
w$dist_tot <- w$dist_max_km + w$dist_SCO_QR

minmax <- pnpsgc3 %>% group_by(POP) %>% summarise(pnpsMin = min(pNpS), pnpsMax=max(pNpS))


Berners=read.table("Berners_piNpiS.4fold.CDS.sumstats.final.clean.GC3bin4Mb",h=TRUE)
summary(lm(formula=pNpS~medianGC3,data=Berners))

Capilano=read.table("Capilano_piNpiS.4fold.CDS.sumstats.final.clean.GC3bin4Mb",h=TRUE)
summary(lm(formula=pNpS~medianGC3,data=Capilano))

Deschutes=read.table("Deschutes_piNpiS.4fold.CDS.sumstats.final.clean.GC3bin4Mb",h=TRUE)
summary(lm(formula=pNpS~medianGC3,data=Deschutes))

Inch=read.table("Inch_piNpiS.4fold.CDS.sumstats.final.clean.GC3bin4Mb",h=TRUE)
summary(lm(formula=pNpS~medianGC3,data=Inch))

Klamath=read.table("Klamath_piNpiS.4fold.CDS.sumstats.final.clean.GC3bin4Mb",h=TRUE)
summary(lm(formula=pNpS~medianGC3,data=Klamath))

Kwethluk=read.table("Kwethluk_piNpiS.4fold.CDS.sumstats.final.clean.GC3bin4Mb",h=TRUE)
summary(lm(formula=pNpS~medianGC3,data=Kwethluk))

Pallant=read.table("Pallant_piNpiS.4fold.CDS.sumstats.final.clean.GC3bin4Mb",h=TRUE)
summary(lm(formula=pNpS~medianGC3,data=Pallant))

Quilcene=read.table("Quilcene_piNpiS.4fold.CDS.sumstats.final.clean.GC3bin4Mb",h=TRUE)
summary(lm(formula=pNpS~medianGC3,data=Quilcene))

Robertson=read.table("Robertson_piNpiS.4fold.CDS.sumstats.final.clean.GC3bin4Mb",h=TRUE)
summary(lm(formula=pNpS~medianGC3,data=Robertson))

Salmo=read.table("Salmon_piNpiS.4fold.CDS.sumstats.final.clean.GC3bin4Mb",h=TRUE)
summary(lm(formula=pNpS~medianGC3,data=Salmo))

Tsooyess=read.table("Tsooyess_piNpiS.4fold.CDS.sumstats.final.clean.GC3bin4Mb",h=TRUE)
summary(lm(formula=pNpS~medianGC3,data=Tsooyess))

########## 
pnps=read.table("allpops_piNpiS_withstatslinearregression.txt",h=TRUE)
ggplot(pnps,aes(x = slope, y =piNpiS))+
  geom_smooth(method="lm", level=0.95, col="darkgrey",fill="lightgrey")+ geom_point(cex=2) + 
  xlab(expression(pi[S]))+ ylab(expression(pi[N]/pi[S]))+
  #scale_linetype_manual(values=c("solid", "dashed"))+
  theme_bw()+
  theme(legend.key.size = unit(2.5,"line"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=30,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=30,angle=90,hjust=.5,vjust=.5,face="italic"))

summary(lm(formula=pnps$slope ~ pnps$piNpiS)) # NS p=0.06583

#Call:
#  lm(formula = pnps$slope ~ pnps$piNpiS)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.08641 -0.02395 -0.01056  0.03925  0.08378 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)    1.959      1.482   1.322   0.2189  
#pnps$piNpiS  -12.128      5.794  -2.093   0.0658 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.05943 on 9 degrees of freedom
#Multiple R-squared:  0.3274,	Adjusted R-squared:  0.2527 
#F-statistic: 4.382 on 1 and 9 DF,  p-value: 0.06583

ggplot(pnps,aes(x = R2, y =piNpiS))+
  geom_smooth(method="lm", level=0.95, col="darkgrey",fill="lightgrey")+ geom_point(cex=2) + 
  xlab(expression(pi[S]))+ ylab(expression(pi[N]/pi[S]))+
  #scale_linetype_manual(values=c("solid", "dashed"))+
  theme_bw()+
  theme(legend.key.size = unit(2.5,"line"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=30,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=30,angle=90,hjust=.5,vjust=.5,face="italic"))

#summary(lm(formula=pnps$R2 ~ pnps$piNpiS)) # NS p=0.69678
#Call:
#  lm(formula = pnps$R2 ~ pnps$piNpiS)

#Residuals:
#  Min         1Q     Median         3Q        Max 
#-0.0226180 -0.0058911  0.0005888  0.0073003  0.0169901 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)   
#(Intercept)   1.0952     0.3220   3.402  0.00785 **
#  pnps$piNpiS  -0.5064     1.2584  -0.402  0.69678   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.01291 on 9 degrees of freedom
#Multiple R-squared:  0.01767,	Adjusted R-squared:  -0.09147 
#F-statistic: 0.1619 on 1 and 9 DF,  p-value: 0.6968


ggplot(pnps,aes(x = Intercept, y =piNpiS))+
  geom_smooth(method="lm", level=0.95, col="darkgrey",fill="lightgrey")+ geom_point(cex=2) + 
  xlab(expression(pi[S]))+ ylab(expression(pi[N]/pi[S]))+
  #scale_linetype_manual(values=c("solid", "dashed"))+
  theme_bw()+
  theme(legend.key.size = unit(2.5,"line"))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=30,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=30,angle=90,hjust=.5,vjust=.5,face="italic"))

summary(lm(formula=pnps$Intercept ~ pnps$piNpiS)) # S p=0.0476
#Call:
#  lm(formula = pnps$Intercept ~ pnps$piNpiS)

#Residuals:
#  Min       1Q   Median       3Q      Max 
#-0.05659 -0.02620  0.00689  0.01641  0.05789 

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)   -1.287      1.002  -1.283   0.2314  
#pnps$piNpiS    8.980      3.918   2.292   0.0476 *
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.04019 on 9 degrees of freedom
#Multiple R-squared:  0.3685,	Adjusted R-squared:  0.2984 
#F-statistic: 5.252 on 1 and 9 DF,  p-value: 0.04763
