
#script to plot the distribution of genotype per major groups

libs <- c('dplyr','ggplot2', 'factoextra','cowplot')
invisible(lapply(libs, library, character.only = TRUE))

setwd("10.deleterious")

high <- read.table("04_genotype_count/coho_biallelic.renamed.high.ann.matrix.genotype.count" , h =T)
moderate   <- read.table("04_genotype_count/coho_biallelic.renamed.moderate.ann.matrix.genotype.count" , h =T)
synonymous <- read.table("04_genotype_count/coho_biallelic.renamed.synonymous.ann.matrix.genotype.count" , h =T)

high$status = "high"
moderate$status = "missense"
synonymous$status = "synonymous"

all <- rbind(synonymous,moderate,high) 
all.m <- reshape2::melt(all)
all.m$POP <- substr(all.m$IND, start = 1, stop = 3)

region <- read.table("region.txt", h = T)

all.m <- merge(all.m, region,by = "POP")
head(all.m)

############# reorder stuff for plotting #######################################
# Reorder factors
genotypes = c("KLA",
              "DES","QUI","TSO",
              "INC" ,"CAP","ROB",
              "PAL",      
              "SAL",
              "KWE","BER","SNA","MSL","POR")

all.m$POP = with(all.m, factor(POP, levels = genotypes))
#all.m$status = with(all.m, factor(status, levels = c("synonymous","moderate","high")))
all.m$REGION = with(all.m, 
                        factor(REGION, 
                               levels = c("California","Cascadia","BC","HaidaGwaii","Thompson","Alaska")))

all.m$variable2 <- ifelse(all.m$variable=="sum_hom", "homozygous load", 
                         ifelse(all.m$variable == "sum_het",  "heterozygous load", "total load"))

all.m$variable2 = with(all.m, 
                      factor(variable2, 
                             levels = c("total load","heterozygous load", "homozygous load")))

all.m$status <- ifelse(all.m$status=="high", "LoF", 
                       ifelse(all.m$status == "synonymous",  "synonymous",  "missense"))
all.m$status = with(all.m, factor(status, 
                                  levels = c("missense","LoF","synonymous")))

table(all.m$status)
table(all.m$variable)


##############"" ---- plotting the data ---------------------------------------#
devtools::install_github("teunbrand/elementalist")
library(elementalist) #to get rounded border

theme_set(theme_gray())

pgeno1  <- filter(all.m, status == "missense", 
                  variable2 != "heterozygous load") %>% 
  ggplot(., aes(x = POP, y = value)) +
  #geom_violin(aes(x = POP, group = POP, colour = factor(REGION)) ) + 
  geom_boxplot(aes(x = POP, group = POP, colour = factor(REGION)) ) + 
  facet_grid(~variable2, scales = "free_y") +
  #personalize grid:
  theme(strip.text.x = element_text(size = 14, face = "bold"),
        strip.background = element_rect_round(radius =  unit(0.5, "cm")))+
  colScale  +
  #theme_classic() +
  ylab("genotype count") +
  xlab("population") +
  theme(axis.title.x = element_text(size = 12, family = "Helvetica", face ="bold"),
        axis.text.x = element_text(size = 12, family = "Helvetica",
                                   face = "bold", angle = 90, hjust = 0, vjust = 0.5),
        axis.title.y = element_text(size = 12, family = "Helvetica",
                                    face = "bold",angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 12,family = "Helvetica",face ="bold"),
        strip.text.x = element_text(size = 12)) +
  theme(legend.position = "none")

pgeno1 + grey

pdf(file = "geno_boxplot_theme_cowplot_round_noleg_grey.pdf", 10,5)
pgeno1
dev.off()

################################################################################
pgeno2  <- filter(all.m, status == "LoF",
                  variable2 != "heterozygous load") %>% 
  ggplot(., aes(x = POP, y = value)) +
  #geom_violin(aes(x = POP, group = POP, colour = factor(REGION)) ) + 
  geom_boxplot(aes(x = POP, group = POP, colour = factor(REGION)) ) + 
  facet_grid(~variable2, scales = "free_y") +
  #personalize grid:
  theme(strip.text.x = element_text(size = 14, face = "bold"),
        strip.background = element_rect_round(radius =  unit(0.5, "cm")))+
  colScale  +
  #theme_classic() +
  ylab("genotype count") +
  xlab("population") +
  theme(axis.title.x = element_text(size = 12, family = "Helvetica", face ="bold"),
        axis.text.x = element_text(size = 12, family = "Helvetica",
                                   face = "bold", angle = 90, hjust = 0, vjust = 0.5),
        axis.title.y = element_text(size = 12, family = "Helvetica",
                                    face = "bold",angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(size = 12,family = "Helvetica",face ="bold"),
        strip.text.x = element_text(size = 12)) +
  theme(legend.position = "none")

pgeno2 +grey

### combine all graph 
p <- plot_grid(pgeno1, pgeno2, nrow = 2, labels = c("A - Missense", "B - LoF"))
svg("Figure3_version3.svg",width = 12, height = 10)
p
dev.off()

##ppt :
p <- plot_grid(pgeno1 + grey , pgeno2 + grey, nrow = 2, labels = c("A - Missense", "B - LoF"))

svg("Figure3_ppt.svg",width = 12, height = 10)
p
dev.off()





# -------------------------- some tests ---------------------------------------#
miss_t <- filter(all.m, status == "missense", variable2 == "total load") 
miss_h <- filter(all.m, status == "missense", variable2 == "homozygous load") 
lof_h <- filter(all.m, status == "LoF", variable2 == "homozygous load") 
lof_t <- filter(all.m, status == "LoF", variable2 == "total load") 

#ks.test( value ~ POP, data = miss_t)
model <- aov(value ~ POP, data = miss_t)
summary(model)
TukeyHSD(model, conf.level=.95)
plot(TukeyHSD(model, conf.level=.95), las = 2)

model <- aov(value ~ POP, data = lof_t)
summary(model)
TukeyHSD(model, conf.level=.95)
plot(TukeyHSD(model, conf.level=.95), las = 2)

model <- aov(value ~ POP, data = miss_h)
summary(model)
TukeyHSD(model, conf.level=.95)
plot(TukeyHSD(model, conf.level=.95), las = 2)

model <- aov(value ~ POP, data = lof_h)
summary(model)
TukeyHSD(model, conf.level=.95)
plot(TukeyHSD(model, conf.level=.95), las = 2)

strata = read.table("../09.DFE/Table_mainresults_Grapes.txt", h = T) %>%
  select(POP, longitude, latitude, LON_emb, LAT_emb, dist_to_ocean, dist_to_source, pSratio, pNratio)

attach(strata)

head(all.m)
agLoF <- all.m %>% 
  filter(variable2 == "homozygous load" & status == "LoF") %>%
  group_by(POP) %>%
  summarise(mean = mean(value))
agMiss <- all.m %>% 
  filter(variable2 == "homozygous load" & status == "missense") %>%
  group_by(POP) %>%
  summarise(mean = mean(value))

agLoF <- merge(agLoF, strata, by = "POP")
agMiss <- merge(agMiss, strata, by = "POP")

head(agLoF)

summary(lm(mean ~ dist_to_source, data = agLoF ))
summary(lm(mean ~ dist_to_source, data = agMiss ))

###test du weighting the miss_t par pS
head(miss_t)
test <- merge(miss_t, strata, by = "POP")
test$weighted <- test$value / test$pSratio  
head(test)
model <- aov(weighted ~ POP, data = test)
summary(model)
TukeyHSD(model, conf.level=.95)
plot(TukeyHSD(model, conf.level=.95), las = 2)
test %>% group_by(POP) %>% summarise(m = mean(weighted))



