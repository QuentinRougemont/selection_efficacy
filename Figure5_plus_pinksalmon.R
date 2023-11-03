
library(ggplot2)
library(dplyr)
library(magrittr)
library(cowplot)
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

coho_gc_cons <- read.table("08.piNpiS/allpop.pnpsGCcons.txxt" , h = F) %>% set_colnames(., c("POP","GCclass", "count","medianGC3","pouet","pNpS"))
model_cohoGC <- coho_gc_cons %>%
  group_by(POP) %>% # You can add here additional grouping variables if your real data set enables it
  do(mod = lm(pNpS ~ medianGC3, data = .)) %>%
  mutate(Slope = summary(mod)$coeff[2], 
         pvalue = summary(mod)$coeff[8],
         R2adj = summary(mod)$adj.r.squared) #%>%
#select(-mod)
model_cohoGC



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
dim(chinook)
dim(rainbow)
dim(coho)
dim(sockeye)
dim(pink)
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
#check the data:
head(pink)
head(chinook)
head(coho)
head(sockeye)
head(rainbow)
head(SOCKEYE)
head(KOKANEE)
#check
dim(coho)
dim(chinook)
dim(sockeye)
dim(rainbow)
dim(SOCKEYE)
dim(KOKANEE)
dim(pink)
pink$POP[pink$POP == "EVEN"] <- "pink even"
pink$POP[pink$POP == "ODD"] <- "pink odd"
coho_gc_cons$species <- "coho"
#combine outgroup for plotting
#Figure 4:
outgroup <- rbind(SOCKEYE, KOKANEE, rainbow, chinook, pink)
dim(outgroup)
dim(coho)
dim(coho_gc_cons)
outgroup$status = "outgroup"
coho$status = "coho"
coho_gc_cons$status = "gccons"
all <- rbind(coho, outgroup) %>% select(POP, medianGC3, pNpS, species, status)
dim(all)

coho_gc_cons <- select(coho_gc_cons, POP, medianGC3, pNpS, species, status)
all <-rbind(all, coho_gc_cons)
  
p <- ggplot(data=all,aes(y=pNpS, x=medianGC3, colour = POP)) +
  facet_grid(~status) +
  stat_smooth(method="lm", se = FALSE) + 
  geom_point(size=3) + 
  theme_bw() + 
  labs(x="median GC3", y =  expression(pi["n"]/pi["s"]) )   + 
  th2 + theme(legend.position = c(0.5,0.75))
p




#for supp fig:
outgroup <- rbind(sockeye, rainbow,pink, chinook)

p <- ggplot(data=outgroup,aes(y=pNpS, x=medianGC3, colour = POP)) +
  facet_grid(~species) + 
  stat_smooth(method="lm") + 
  geom_point(size=3) + 
  theme_bw() + 
  labs(x="median GC3", y =  expression(pi["n"]/pi["s"]) )   + 
  th2
p <- p + theme(legend.position = "none")
p

p_grey <- p + grey #for powerpoint only


coho <- merge(coho, strata)
coho$region <- as.factor(coho$region)

myColors <- c("blue","orange","red","green","darkviolet","springgreen4") #outgroups,"black","black","grey")
names(myColors) <- levels(coho$region)
colScale <- scale_colour_manual(name = "region",values = myColors)

p_coho <-  ggplot(data=coho,aes(y=pNpS, x=medianGC3, colour = region)) +
  facet_grid(~species) + 
  stat_smooth(method="lm") + 
  geom_point(size=3) + 
  theme_bw() + colScale +
  labs(x="median GC3", y =  expression(pi["n"]/pi["s"]) )   + 
  th_plot

p_coho_grey <- p_coho + grey


all <- plot_grid(p_coho_grey, p_grey, labels = "AUTO", rel_widths = c(0.5,1))
all

ggsave(filename = "all.svg", all, units = "cm", width = 28, height = 15, dpi = 400)
getwd()


pdf(file="pNps_rainbow_trout.pdf")
p1
dev.off()
