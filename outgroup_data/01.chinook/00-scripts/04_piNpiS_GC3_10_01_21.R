
#script to compute relationship between pnsps and GC3 as well as exploring other possible correlation
#date: 19-01-2021
#author: QR

library(ggplot2)
library(dplyr)

#argv <- commandsArgs(T)
#input <- argv[1] 
setwd("/home/quentin/postdoc/coho/wgs/2021/01.analysis/12pinpis/new/04-plot/")
setwd("/home/quentin/postdoc/coho/wgs/2021/01.analysis/12pinpis/new/DIPLOID/04-plot/")

input <- "pnps_gc3.6000000.txt"
#input <- "pnps_gc3.200000.txt"
input <- "pnps_gc3.500000.txt"
pnpsgc3 <- read.table(input, T)

#perform linear models by pop
models <- pnpsgc3 %>%
  group_by(POP) %>% # You can add here additional grouping variables if your real data set enables it
  do(mod = lm(pNpS ~ medianGC3, data = .)) %>%
  mutate(Slope = summary(mod)$coeff[2], 
         pvalue = summary(mod)$coeff[8],
         R2adj = summary(mod)$adj.r.squared) #%>%
models <- dplyr::select(models, -mod)

str <-read.table("../../strata3.txt", T)
str <-read.table("../strata3.txt", T)
models <- merge(models, str)
#compute the distance to the southernmost site
models$dist_tot <- models$dist_max_km + models$dist_SCO_QR

minmax <- pnpsgc3 %>% group_by(POP) %>% 
  summarise(pnpsMin = min(pNpS), pnpsMax=max(pNpS))
models <- cbind(models,minmax)
#correlation between min max and distance to the source:
summary(lm(models$dist_tot ~ models$pnpsMin)) #signif
summary(lm(models$dist_tot ~ models$pnpsMax)) #not signif
plot(models$dist_tot, models$pnpsMax)
plot(models$dist_tot, models$pnpsMax)
text(models$dist_tot, models$pnpsMax, models$POP)

#correlation between the slope of the relationshpe and distance to the source?
summary(lm(models$Slope ~ models$dist_tot)) #no correlation
plot(models$Slope, models$dist_tot)

plot(models$Slope, models$pnpsMax)
text(models$Slope, models$pnpsMax, models$POP)

plot(models$dist_tot, models$pnpsMax)
text(models$dist_tot, models$pnpsMax, models$POP)

minmaxGC <- pnpsgc3 %>% group_by(POP) %>% 
  summarise(gc3Min = min(medianGC3), gc3Max=max(medianGC3))
models <- cbind(models,minmaxGC)
summary(lm(models$gc3Min ~ models$dist_tot)) #correlation
plot(models$dist_tot , models$gc3Min)
text(models$dist_tot , models$gc3Min, models$POP)
summary(lm(models$gc3Min ~ models$dist_tot + I(models$dist_tot^2) + I(models$dist_tot^3)))
summary(lm(models$gc3Min ~ models$dist_tot + I(models$dist_tot^2) ))
models = models[,-1]

#lot of uncertainty!
myColors <- c("blue","orange",  "red","green",
              "darkviolet", "springgreen4")
models$region<-as.factor(models$region)
names(myColors) <- levels(models$region)
colScale <- scale_colour_manual(name = "region",values = myColors)

p <- ggplot(data=models,
            aes(y=gc3Min, x=dist_tot))
p <- p + stat_smooth(method="lm")
#p <- p + stat_smooth(method = "lm", formula = y ~ x + I(x^2)+ + I(x^3), size = 1)
p <- p + geom_point(aes(color=factor(region)), size=3)
p <- p + colScale
p <- p + theme_bw()
p <- p + labs(x="distance", y = "GC3min")  #+ 
  #expand_limits(y=0.00, x=1) + scale_x_continuous(limits = c(0, 1))
p <- p + theme(axis.title.x=element_text(size=14, family="Helvetica",face="bold"),
               axis.text.x=element_text(size=18,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
               axis.title.y=element_text(size=20, family="Helvetica",face="bold",angle=0, hjust=0, vjust=0.5),
               axis.text.y=element_text(size=14,family="Helvetica",face="bold"),
               strip.text.x = element_text(size=18))
p  <- p + theme(panel.grid.major = element_blank() )#, panel.grid.minor = element_blank())





#merge pnps and slope
pnps <- read.table("../pnps.tetra_19.01.txt",T)
pnps <- read.table("../pnps.overall_19_01_20.txt",T)
pnps <- read.table("../pnps.diploidonly.txt",T )
models <- models[,-1]
m <-left_join(models, pnps, by=c("POP"="myfile"))
summary(lm(m$Slope ~ m$pNpSratio)) #signif 
summary(lm(m$Slope ~ m$pSratio)) #signif 

#### plot correlation between GC3 and piNpiS ###
myColors <- c("blue","orange",
              "red","green",
              "darkviolet",
              "springgreen4")
m$region<-as.factor(m$region)
names(myColors) <- levels(m$region)
colScale <- scale_colour_manual(name = "region",values = myColors)

p <- ggplot(m,aes(x =pSratio, y =Slope))+
  geom_smooth(method="lm", fill="gray87")+
  geom_point(cex=2, (aes(colour=factor(region)))) + 
  xlab(expression(pi[S]))+ 
  ylab(expression(slope (pi[N]/pi[S] ~~~GC3)))+
  #ylab(paste0("slope", expression(pi[N]/pi[S])))+
  colScale+
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
p

pdf(file="slope_pnps_ps.pdf")
p
dev.off()



##### corr median GC3 ~ pNpS
pnpsgc3 <-  merge(pnpsgc3, str)
summary(lm(pnpsgc3$medianGC3 ~ pnpsgc3$pNpS))

pnpsgc3$region<-as.factor(pnpsgc3$region)
names(myColors) <- levels(pnpsgc3$region)
colScale <- scale_colour_manual(name = "region",values = myColors)

p <- ggplot(pnpsgc3,aes(x = medianGC3, y =pNpS ,colour=factor(region)))+
  geom_point(cex=2, (aes(colour=factor(region)))) + 
  geom_smooth(method="lm", fill=NA)+
  xlab("GC3")+ ylab(expression(pi[N]/pi[S]))+
  colScale+
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

pdf(file="correlation_pnps_gc3b_all.6mb.pdf",9,9)
p
dev.off()


m$POP <-as.factor(m$POP)
m$POP <- as.factor(m$POP)
myColors <- c("darkblue","orange","green","darkorange", "red",
              "blue","#003366","darkviolet","#3399FF",
              "#99FF33","chocolate","springgreen4","#0000FF",
              "green")

names(myColors) <- levels(m$POP)
colScale <- scale_colour_manual(name = "region",values = myColors)

p <- ggplot(pnpsgc3,aes(x = medianGC3, y =pNpS,colour=POP))+
  geom_point(cex=2) + geom_smooth(method="lm", fill=NA)+
  xlab("GC3")+ ylab(expression(pi[N]/pi[S]))+colScale+
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

pdf(file="pnps_gc3_all_6mb.pdf",10,8)
p
dev.off()

p1 <- p+ annotate(geom="text", 
                  x=0.55, y=0.14, 
                  label= "paste(italic(R) ^ 2, \" =0.94, p <1e-16***, slope = -0.64\")", 
                  parse = TRUE,
                  color="black", size=5)

pdf(file="pnps_gc3_all_6mb_with_corr.pdf",10,8)
p1
dev.off()


#faire plot pnps slope
#faire plot R2 piNpiS
#faire plot intercept pinpis?
####################################################################################
