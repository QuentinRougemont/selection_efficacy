
library(ggplot2)
library(dplyr)
library(magrittr)
library(cowplot)
library(tidyr)
###################### loading data for rainbow trout ##############################################
### rainbow

rainbow = read.table("outgroup_data/02.rainbow/all_withoutquantiles.CDS.sumstats.4000000_mb", h = T)  #the data for figure 4
rainbow_pop = read.table("outgroup_data/02.rainbow/all.4mb.txt", h = T)      #the data for figure supp (separated by pop)

models <- rainbow_pop %>%
  group_by(POP) %>% # You can add here additional grouping variables if your real data set enables it
  do(mod = lm(pNpS ~ medianGC3, data = .)) %>%
  mutate(Slope = summary(mod)$coeff[2], 
         pvalue = summary(mod)$coeff[8],
         R2adj = summary(mod)$adj.r.squared) %>%
  select(-mod)

models

### chinook
chinook <- read.table("outgroup_data/01.chinook/POP_withoutquantiles.CDS.sumstats.4000000_mb" , h = T)
mod <- lm(pNpS ~ medianGC3, data = chinook) #0.9338
slope = summary(mod)$coeff[2] 
pval = summary(mod)$coeff[8]
R2 = summary(mod)$adj.r.squared
slope
pval
R2

sockeye <- read.table("outgroup_data/03.sockeye//all.CDS.4fold.4mb.txt", h = T)

model_sock <- sockeye %>%
  group_by(POP) %>% # You can add here additional grouping variables if your real data set enables it
  do(mod = lm(pNpS ~ medianGC3, data = .)) %>%
  mutate(Slope = summary(mod)$coeff[2], 
         pvalue = summary(mod)$coeff[8],
         R2adj = summary(mod)$adj.r.squared) #%>%
  #select(-mod)

model_sock

KOKANEE <- read.table("outgroup_data/03.sockeye/KOKANEE_withoutquantiles.4fold.CDS.sumstats.4000000_mb", h = T) #for figure 4 only kokanne from chilo
SOCKEYE <- read.table("outgroup_data/03.sockeye/SOCKEYE_withoutquantiles.4fold.CDS.sumstats.4000000_mb", h = T) #for figure 4 only Socekye from Meadows
summary(lm(KOKANEE$medianGC3 ~ KOKANEE$pNpS))
summary(lm(SOCKEYE$medianGC3 ~ SOCKEYE$pNpS))


pink <- read.table("outgroup_data/04.pink/pink.even.odd.GC3_pnps.txt", h = T)
model_pink <- pink %>%
  group_by(POP) %>% # You can add here additional grouping variables if your real data set enables it
  do(mod = lm(pNpS ~ medianGC3, data = .)) %>%
  mutate(Slope = summary(mod)$coeff[2], 
         pvalue = summary(mod)$coeff[8],
         R2adj = summary(mod)$adj.r.squared) #%>%
#select(-mod)

model_pink

coho <- read.table("08.piNpiS/pnps.allCDS.01_03_21.txt.gz", h = T)
coho <- read.table("08.piNpiS/pnps_gc3.4000000.txt.gz", h = T)
coho <- read.table("08.piNpiS/pnps_gc3.4000000.01.03.txt.gz", h = T)

strata <- read.table("08.piNpiS//strata3.txt", h = T) %>% select(POP, region)
head(coho)

model_coho <- coho %>%
  group_by(POP) %>% # You can add here additional grouping variables if your real data set enables it
  do(mod = lm(pNpS ~ medianGC3, data = .)) %>%
  mutate(Slope = summary(mod)$coeff[2], 
         pvalue = summary(mod)$coeff[8],
         R2adj = summary(mod)$adj.r.squared) #%>%
#select(-mod)
model_coho

## ----------------------------------------------------------------------------##
##correlation between pnps from GC conservative site and GC3:
coho_gc_cons <- read.table("08.piNpiS/allpop.pnpsGCcons.txxt" , h = F) %>% 
    set_colnames(., c("POP","bin", "n", "piN", "piS", "pNpS"))

df <-coho_gc_cons[order(coho_gc_cons$POP, coho_gc_cons$n),]

#on recupere le GC3 #this code needs to be updated #not clean 
w <- rep(sort(rep(seq(1, 14), 2)), 14)
coho$vec=paste0(coho$POP,"_", w)
medgc3 <- aggregate(coho$medianGC3, by=list(coho$vec), mean)
medgc3 <- tidyr::separate(medgc3, Group.1, c("POP","bin"), sep="_")
medgc3pnps <- merge(df, medgc3)

cor(medgc3pnps$pNpS, medgc3pnps$x) #juste to look at correlation between pnps and gc3 broadly

colnames(medgc3pnps)[7] <- "medianGC3"
head(medgc3pnps)

model_cohoGC <- medgc3pnps %>%
  group_by(POP) %>% # You can add here additional grouping variables if your real data set enables it
  do(mod = lm(pNpS ~ medianGC3, data = .)) %>%
  mutate(Slope = summary(mod)$coeff[2], 
         pvalue = summary(mod)$coeff[8],
         R2adj = summary(mod)$adj.r.squared) #%>%
#select(-mod)
model_cohoGC

medgc3pnps <- merge(medgc3pnps, strata)
medgc3pnps$POP <- as.factor(medgc3pnps$POP)


#all data are ready for plotting now

#------------------------------------------------------- PLOT ------------------------------------------------------------#
### define general theme:
th_plot <-     theme(axis.title.x=element_text(size=14, family="Helvetica",face="bold"),
                     axis.text.x=element_text(size=14,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
                     axis.title.y=element_text(size=20, family="Helvetica",face="bold",angle=0, hjust=0, vjust=0.5),
                     axis.text.y=element_text(size=14,family="Helvetica",face="bold"),
                     strip.text.x = element_text(size=18),
                     panel.grid.major = element_blank(),
                     legend.position = "none")

th2 <-   theme_bw() +
  theme(legend.key.size = unit(2.5,"line"),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = 'black', size = 1.25), axis.ticks = element_line(colour = 'black', size = 1.25), 
        axis.text.x = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="black",size=20,angle=0,hjust=.5,vjust=.5,face="plain"),  
        axis.title.x = element_text(colour="black",size=30,angle=0,hjust=.5,vjust=.2,face="italic"),
        axis.title.y = element_text(colour="black",size=30,angle=90,hjust=.5,vjust=.5,face="italic"))

##grey colors for ppt ONLY:
grey <- theme(panel.background = element_rect(fill = "gray20", colour = "lightblue", size = 0.5, linetype = "solid"),
              plot.background = element_rect(fill = "gray20"),
              axis.title.x=element_text(size=22, family="Helvetica",face="bold", color="#CCCCCC"),
              axis.text.x=element_text(size=14, color = "#CCCCCC", family="Helvetica",face = "bold", angle=90, hjust=0, vjust=0.5),
              axis.title.y=element_text(size=22, family="Helvetica", color="#CCCCCC", face="bold",angle=90,  hjust=0.5, vjust=0.5),
              axis.text.y=element_text(size=14,family="Helvetica", color="#CCCCCC", face="bold"),   strip.text.x = element_text(size=18),
              panel.grid.major = element_blank()) #, panel.grid.minor = element_blank())

#------------------------------------------------------------------------------------------------------------------#
#------------------------------------ now plot --------------------------------------------------------------------#
############## rainbow only ############################
coho <- merge(coho, strata)
myColors <- c("blue","orange","red","green","darkviolet","springgreen4")

coho$POP <- as.factor(coho$POP)
coho$region <- as.factor(coho$region)
levels(coho$region)
names(myColors) <- levels(coho$region)

colScale <- scale_colour_manual(name = "population",values = myColors)

head(coho)

plt <- ggplot(coho,aes(x = medianGC3, y = pNpS, colour = region))+
  geom_smooth(method="lm", fill=NA)+
  geom_point(cex=2) + 
  xlab("GC3")+ ylab(expression(pi[N]/pi[S])) +
  colScale+
  th2 +
  guides(col=guide_legend(ncol = 3)) + 
  theme(legend.position = c(0.5,0.75))

plt <-  plt + annotate(geom="text", 
                                     x=0.58, y=0.2, 
                                     label= "paste(italic(R) ^ 2, \" = 0.62 - 0.96 p <2e-16\")", 
                                     parse = TRUE,
                                     color="black", size=6, face = "plain")

plt


head(medgc3pnps)
myColors <- c("blue","orange","green","orange",
              "red","blue","blue","darkviolet",
              "blue","green","orange",
              "springgreen4","blue","green")
names(myColors) <- levels(coho$POP)
colScale <- scale_colour_manual(name = "POP",values = myColors)

plt_gccons <- ggplot(medgc3pnps,aes(x = medianGC3, y =pNpS, colour=POP))+
  geom_smooth(method="lm", fill=NA)+
  geom_point(cex=2) + 
  xlab("GC3")+ ylab(expression(SFSbased * ~ pi[N]/pi[S] ** GCcons)) +
  colScale+
  th2 +
  theme(legend.position = "none")

plt_gconcs <-  plt_gccons + annotate(geom="text", 
                     x=0.58, y=0.2, 
                     label= "paste(italic(R) ^ 2, \" = 0.665 - 0.876 p <2e-16\")", 
                     parse = TRUE,
                     color="black", size=6, face = "plain")
plt_gccons

KOKANEE$POP <- "kokannee"
SOCKEYE$POP <- SOCKEYE$species <- "sockeye" 
KOKANEE$species <- "sockeye"
chinook$POP <- chinook$species <- "chinook"
rainbow$POP <- rainbow$species <- "rainbow"
sockeye$species <- "sockeye"
coho$species <- "coho"
pink$species <- "pink"

#reorder pink and chinook to match other dataset
chinook <- chinook[,c(10,1:9)]
rainbow <- rainbow[,c(10,1:9)]
KOKANEE <- KOKANEE[,c(10,1:9)]
SOCKEYE <- SOCKEYE[,c(10,1:9)]
pink$POP[pink$POP == "EVEN"] <- "pink even"
pink$POP[pink$POP == "ODD"] <- "pink odd"
coho_gc_cons$species <- "coho"

#combine outgroup for plotting
#Figure 4:
outgroup <- rbind(SOCKEYE, KOKANEE, rainbow, chinook, pink)
colnames(outgroup)[1] <- "subspecies"
p <- ggplot(data=outgroup ,aes(y=pNpS, x=medianGC3, colour = subspecies)) +
  stat_smooth(method="lm", se = FALSE) + 
  geom_point(size=3) + 
  theme_bw() + 
  labs(x="GC3", y =  expression(pi["n"]/pi["s"]) )   + 
  th2 + theme(legend.position = c(0.5,0.9)) +
  guides(col=guide_legend(ncol=2)) 
  
p <- p +  annotate(geom="text", 
                    x=0.58, y=0.2, 
                    label= "paste(italic(R) ^ 2, \" = 0.91 - 0.97 p <2e-16\")", 
                    parse = TRUE,
                    color="black", size=6, face = "plain")
p


all <- plot_grid(plt, plt_gccons, p, ncol = 3, labels = "AUTO") #, rel_widths = c(0.5,1))
all

ggsave(filename = "Figure4.svg", all, units = "cm", width = 40, height = 15, dpi = 400)
getwd()
