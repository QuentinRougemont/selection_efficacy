

#purpose: construct figure 2 from paper.
# data obtained from pNpS pipeline:
# 1 - https://github.com/QuentinRougemont/PiNPiS
# 2 - smc++: https://github.com/QuentinRougemont/smcpp_input
# 3 - grapes: to add to github
#note: whole genome generated for coho and outgroup using gatk with the pipeline here:
#https://github.com/QuentinRougemont/gatk_haplotype


#------------------------------------------------- check if scripts are installed and load them -----------------------------------------#
if("dplyr" %in% rownames(installed.packages()) == FALSE)
{install.packages("dplyr", repos="https://cloud.r-project.org") }
if("ggplot2" %in% rownames(installed.packages()) == FALSE)
{install.packages("ggplot2", repos="https://cloud.r-project.org") }
if("data.table" %in% rownames(installed.packages()) == FALSE)
{install.packages("data.table", repos="https://cloud.r-project.org") }
if("magrittr" %in% rownames(installed.packages()) == FALSE)
{install.packages("maggritr", repos="https://cloud.r-project.org") }
if("cowplot" %in% rownames(installed.packages()) == FALSE)
{install.packages("cowplot", repos="https://cloud.r-project.org") }

## load libs
libs <- c('dplyr','ggplot2','data.table','magrittr', 'cowplot')
invisible(lapply(libs, library, character.only = TRUE))

#-------------------------------------------------   load the data and metadata -----------------------------------------#

setwd("08.piNpiS")
region <- read.table("strata3.txt",T)                   #various metadata

region$dist_tot = region$dist_max_km + region$dist_SCO_QR #distance to the southern site + distance to water mouth

pnps = read.table("pnps.allCDS.01_03_21.txt.gz",T)        #pnps value

region <- select(region, POP, region, dist_tot, N_recent) #N_recent = smc++ ne over 2e5years
pnps <- merge(pnps, region)

summary(lm(pnps$pNpSratio ~ pnps$dist_tot)) #test correlations


#------------------------------------------------------- PLOT ------------------------------------------------------------#
### define general theme:
th_plot <-     theme(axis.title.x=element_text(size=14, family="Helvetica",face="bold"),
                     axis.text.x=element_text(size=14,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
                     axis.title.y=element_text(size=20, family="Helvetica",face="bold",angle=0, hjust=0, vjust=0.5),
                     axis.text.y=element_text(size=14,family="Helvetica",face="bold"),
                     strip.text.x = element_text(size=18),
                     panel.grid.major = element_blank(),
                     legend.position = "none")

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
pnps$region<-factor(pnps$region)
myColors <- c("blue","orange","red","green","darkviolet","springgreen4") #outgroups,"black","black","grey")
names(myColors) <- levels(pnps$region)
colScale <- scale_colour_manual(name = "region",values = myColors)

p <- ggplot(data=pnps,
            aes(x=dist_tot, y=pNpSratio)) + 
  stat_smooth(method="lm", fill="gray87") + 
  theme_bw() +
  geom_point(aes(colour=factor(region)),size=4) + 
  colScale + 
  labs(x="Distance from the southernmost site (km)", y = expression(pi["n"]/pi["s"]))  +
   scale_x_continuous(limits = c(0, 8000)) +
   scale_y_continuous(limits = c(0.26, 0.30)) #pour NW
   #scale_y_continuous(limits = c(0.20, 0.30)) #pour full

p <- p + th_plot
p1 <- p + annotate(geom="text", 
                  x=3650, y=0.295, 
                  label= "paste(italic(R) ^ 2, \" = .73 p <0.0001***, slope = 1.7e-06\")", 
                  parse = TRUE,
                  color="black", size=5)
#p
p1


###---------------------------------ajouter panel B) omegaNA distance --------------------------------------------------
grape = read.table("../09.DFE/Table_mainresults_Grapes.txt",T)
grape$dist_tot=grape$dist_to_ocean+grape$dist_to_source

# test many possible correlation ------------------------------------------------
attach(grape)
cor(allsites_3sp_omegaNA, dist_to_ocean )
cor(allsites_3sp_omegaNA, dist_to_source )
cor(allsites_3sp_omegaNA, dist_tot )

summary(lm(allsites_3sp_omegaNA ~ dist_to_ocean ))
summary(lm(allsites_3sp_omegaNA ~ dist_to_source ))
summary(lm(allsites_3sp_omegaNA ~ dist_tot ))

summary(lm(allsites_3sp_omegaA ~ dist_to_ocean ))
summary(lm(allsites_3sp_omegaA ~ dist_to_source ))
summary(lm(allsites_3sp_omegaA ~ dist_tot ))


colnames(grape)
summary(lm(allsites_3sp_alpha ~ dist_tot))
summary(lm(allsites_3sp_alpha ~ dist_to_ocean))
summary(lm(allsites_3sp_alpha ~ dist_to_ocean * dist_to_source))
summary(lm(allsites_3sp_omegaNA ~ dist_to_ocean * dist_to_source))
summary(lm(allsites_3sp_omegaNA ~ latitude  * longitude))

summary(lm(allsites_3sp_omegaA ~ dist_tot))

#--------------------------------------------------- DO THE PLOT ---------------------------------------------------------------------
myColors <- c("blue","orange",
              "red","green",
              "darkviolet",
              "springgreen4")
grape$region<-as.factor(grape$region)
names(myColors) <- levels(grape$region)
colScale <- scale_colour_manual(name = "region",values = myColors)

#relation entre piS et omegaA
summary(lm(allsites_3sp_alpha ~ pSratio))
p3 <- ggplot(data=grape,
            aes(x=pSratio, y=allsites_3sp_alpha)) +
  stat_smooth(method="lm", fill="gray87") + 
  theme_bw() +
  geom_point(aes(colour=factor(region)),size=4) +
  colScale + 
  labs(x=expression(pi["s"]), y = expression(alpha)) # +

p3 <- p3 + th_plot

p3b <- p3+ annotate(geom="text", 
                  x=0.001, y=-0.05, 
                  label= "paste(italic(R) ^ 2, \" = .63 p = 0.0004\")", 
                  parse = TRUE,
                  color="black", size=5)
p3b

#--------------------------------------------------- relationship between smc++ Ne and alpha ------------------------------------------ #
#### relation entre Ne et alpha #####

summary(lm(pnps$pSratio ~ Ne_estimate))
summary(lm(pnps$pNpSratio ~ Ne_estimate))
summary(lm(allsites_3sp_alpha ~ Ne_estimate))

#on fait le graphe de la relation:
p4 <- ggplot(data=grape,
            aes(x=Ne_estimate, y=allsites_3sp_alpha)) + 
  stat_smooth(method="lm", fill="gray87") + 
  theme_bw() + 
  geom_point(aes(colour=factor(region)),size=4) + 
  colScale +
  labs(x=expression(N["e"]), y = expression(alpha)) # +
p4 <- p4 + th_plot
p4 <- p4 + annotate(geom="text", 
                  x=30000, y=0.05, 
                  label= "paste(italic(R) ^ 2, \" = .49 p = 0.004\")", 
                  parse = TRUE,
                  color="black", size=5)
p4

#---------------------------------------- relationship omegaNa and distance ---------------------------------------------
summary(lm(allsites_3sp_omegaNA ~ dist_to_ocean * dist_to_source)) #very cool with interaction
summary(lm(allsites_3sp_omegaNA ~ dist_tot))

p2 <- ggplot(data=grape,
            aes(x=dist_to_source, y=allsites_3sp_omegaNA)) + 
  stat_smooth(method="lm", fill="gray87") + 
  theme_bw() + 
  geom_point(aes(colour=factor(region)),size=4) + 
  colScale + 
  labs(x="distance from the southernmost site (km)", y = expression(omega["NA"]))  +
  scale_x_continuous(limits = c(0, 8000)) +
  scale_y_continuous(limits = c(0.15,0.30))

p2 <- p2 + th_plot
p2 <- p2 + annotate(geom="text", 
                    x=4000, y=0.29, 
                    label= "paste(italic(R) ^ 2, \" = .56 p = 0.0011\")", 
                    parse = TRUE,
                    color="black", size=5)

############### combine all "####################################################

pdf(file = "load.pdf", 12,7)
plot_grid(p1,p2,p3b,p4, nrow = 2, labels = "AUTO")
dev.off()
