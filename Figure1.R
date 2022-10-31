

#purpose: script to reproduce our figure 1
#Plot Tajima's D from ANGSD and sampling maps. 
#date: 20-05-22
#Author: QR
#to generate the tajima's D value I used the script available here:
#https://github.com/QuentinRougemont/utility_scripts/tree/master/00.ANGSD

#for each fourteen populations:
#1 - construct saf file: https://github.com/QuentinRougemont/utility_scripts/blob/master/00.ANGSD/00_angsd_bam_to_saf.sh
#2 - saf 2 1d-sfs: https://github.com/QuentinRougemont/utility_scripts/blob/master/00.ANGSD/01_1dsfs.sh
#3 - sfs2theta : https://github.com/QuentinRougemont/utility_scripts/blob/master/00.ANGSD/03_sfs_to_theta.sh
#then reshaped the data in bash to obtain the final output

# ---------------------------- load libs ----------------------------------------------------- #
libs <- c('ggplot2','dplyr','data.table','cowplot','maps','mapplots','gridGraphics', 'rgdal','GISTools')
invisible(lapply(libs, library, character.only = TRUE))

# ----------------------------- load data ---------------------------------------------------- #
tajimas <- fread("zcat 07.tajimasD/ALL.global.tajimas.txt.gz")
strata <- read.table("strata3.txt.gz", h = T)

colnames(tajimas)[1] = "populations"

#add pop data:
taj <- merge(tajimas, strata, by.x="populations", by.y="POP")

#------------------------------------------------------------------------------------------------------------------#
#look at various correlations: 
summary(lm(taj$Tajima ~ taj$dist_SCO_QR)) #distance to the south
summary(lm(taj$Tajima ~ taj$dist_max_km)) #distance to the river mouth
taj$tot = taj$dist_SCO_QR + taj$dist_max_km
summary(lm(taj$Tajima ~ taj$tot))

#test for any interaction between latitude and longitude:
summary(lm(taj$Tajima ~ taj$longitude + taj$latitude + taj$longitude*taj$latitude))

#-------------------- prepare data for plotting ------------------------------------------------------------------#
target <- data.frame(c("KLA","DES","TSO","QUI","CAP","INC","ROB", "PAL","SAL",
            "BER","KWE","SNA","MSL","POR"))
colnames(target) <- "populations"
taj <- left_join(target, taj)

myColors <- c("blue","orange","red","green","darkviolet","springgreen4")
taj$region<-as.factor(taj$region)
names(myColors) <- levels(taj$region)
colScale <- scale_colour_manual(name = "region",values = myColors)

taj$populations <- factor(taj$populations, levels = unique(taj$populations, sort=F))

#Tajimas'D
p <- ggplot(data=a, aes(x=populations,y=Tajima, fill=region)) +
  geom_violin(trim=FALSE) +  scale_fill_manual(values=myColors) + colScale +
  geom_boxplot(width=0.1, fill="white")+
  labs(title = " ", y = "Tajimas'D") +
  theme_minimal()
p <- p + theme(axis.title.x=element_text(size=12, family="Helvetica",face="bold"),
               axis.text.x=element_text(size=10,family="Helvetica",face="bold", angle=45, hjust=0, vjust=0.5),
               axis.title.y=element_text(size=12, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
               axis.text.y=element_text(size=10,family="Helvetica",face="bold"),
               strip.text.x = element_text(size=10),
               plot.title.position = "plot",
               plot.title = element_text(size = 8, family = "Helvetica", face = "bold")) #,   #NEW parameter. Apply for subtitle too.
#panel.grid.major = element_blank() )#, panel.grid.minor = element_blank())
p <- p + coord_flip() + theme(legend.position = "none")


#------------------------------------ construct the map ---------------------------------------- #

coord <- unique(dplyr::select(taj, populations,longitude,latitude, region))          #take the coordinate
#riversData <- readOGR("00.raw_data/down236_latlon.shp")  # load the shapefile
riversData <- readOGR("down236_latlon.shp")  # load the shapefile
waterColor <- 'cadetblue'

#note: It is possible to do the map with ggplot2, but I never liked the display so in the end get back to base R.
par(mar=c(3, 3, 2, 2))
plot(coord[,c(2,3)],
     axes=T,
     asp = 1,
     lty=2, 
     lwd=1,
     bg = waterColor, 
     col="#ffffff",
     #main = "A",
     xlab = "Longitude",
     ylab = "Latitude")
map(add = T, col = "grey90", fill = TRUE)    
maps::map.scale(x=-140, y=45, ratio=FALSE, relwidth=0.2)  
north.arrow(xb=-147.5, yb=45, len=0.35, lab="N")  

plot(riversData, col=waterColor, add=T) # plot rivers

points(coord[,2],coord[,3], 
       col=c("red","green","green","green","orange","orange","orange",
             "darkviolet","darkgreen",c(rep("blue",5)) ),
       cex= 2, pch=16)
text(coord[,2],coord[,3], labels=coord[,1], cex=1.1, pos=3, font = 4)
title(xlab = "Longitude", ylab = "Latitude")
#title(main = "A", adj = 0)

legend(
       legend = c("Alaska","BC","California","Cascadia","HaidaGwai","Thompson"),
       col = c("blue","orange","red","green","darkviolet","darkgreen"),
       box.lwd = 0,
       pch = 16,
       x = -160, y = 55)

map.axes()

pMap <- recordPlot() #this capture the plot for combination with ggplot2 graph


#------------------------ Creation panel C -------------------------------------
countries <- c("Canada","USA")
maps <- map_data("world", region = countries)

sub.maps <- filter(maps, long < -10 & lat < 75 & lat >25) %>%
  mutate(region=replace(region, region=="USA", "Canada"))

inset2 <- ggplot(sub.maps, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group,  fill = region))+
  geom_point(data = coord, aes(x=longitude, y = latitude, color = region)) +
  scale_fill_grey()+
  colScale +
  theme_void()+
  labs(x="Longitude", y = "Latitude", main = "C")+
  theme(legend.position = "none") +
  theme(axis.title.x=element_text(size=8, family="Helvetica",face="bold"),
                 axis.text.x=element_text(size=8,family="Helvetica",face="bold", angle=45, hjust=0, vjust=0.5),
                 axis.title.y=element_text(size=8, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
                 axis.text.y=element_text(size=8,family="Helvetica",face="bold"),
                 strip.text.x = element_text(size=8),
                 plot.title.position = "plot",
                 plot.title = element_text(size = 8, family = "Helvetica", face = "bold"))  #NEW parameter. Apply for subtitle too.
inset2

#combine the two panel
p1p2 <- ggdraw(pMap) +
  draw_plot(inset2, .7, .60, .25, .25)

p1p2

#add tajimas'D panel now:
all <- plot_grid(p1p2, p, rel_widths = c(2,1), ncol = 2, labels = "AUTO")
all

ggsave(filename = "Fig1.Carte_BaseR_ggplot2.svg", all, units = "cm", dpi = 400, height = 20, width = 34)

#then I've slightly move the river label in inkscape since ggrepel or other stuff was not appropriate here.
