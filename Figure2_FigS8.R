### authors: QR
### date: 2023
### purpose: Reproduce figure2 of our plos genetics paper

library(ape)
library(phytools)
library(stats)
library(dplyr)
library(magrittr)
library(ggplot2)
library(adephylo)
library(ggpubr)
library(mapdata)
library(cowplot)
library(grid)
library(gridGraphics)

#-------------------------------------------------------------------------------#
rm(list = ls())

#-------------------------------------------------------------------------------#
#usual plotting stuff for ggplot2:
th_plot <-     theme(axis.title.x=element_text(size=14, family="Helvetica",face="bold"),
                     axis.text.x=element_text(size=14,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
                     axis.title.y=element_text(size=18, family="Helvetica",face="bold",angle=90, hjust=0.5, vjust=0.5),
                     axis.text.y=element_text(size=14,family="Helvetica",face="bold"),
                     strip.text.x = element_text(size=18),
                     panel.grid.major = element_blank()) 


#-------------------------------------------------------------------------------#
### dowload all the data
tablefst <- "10.fst/pairwise.fst.txt"
pop_map = "00.raw_data/strata.txt"

grape = read.table("09.DFE/grape", h = T) #data with all metadata

data = read.table(tablefst, h = T)
pop_map = (read.table(pop_map))
pop_map = data.frame(unique(pop_map[,2]))

#-------------------------------------------------------------------------------#
mat_correlation = as.data.frame(matrix(999,nrow(pop_map),nrow(pop_map)) )

colnames(mat_correlation) = pop_map[,1]
rownames(mat_correlation) = pop_map[,1]


# buid the fst correlation matrix
print("building the fst correlation matrix...")
for (i in 1:length(data[,1])){
  mat_correlation[data[i,1],data[i,2]] = data[i,3]   
  mat_correlation[data[i,2],data[i,1]] = data[i,3]
  for (colonne in 1:ncol(mat_correlation)){
    for (ligne in 1:nrow(mat_correlation)){
      if (colonne == ligne){
        mat_correlation[colonne,ligne] = 0
      }
    }
  }
  mat_correlation = round(mat_correlation, digits = 5)  
}  

fst.matrix <- as.matrix(mat_correlation)

#-------------------------------------------------------------------------------#
## Fst tree ##
tree <- ape::nj(X = fst.matrix) # fst.matrix as a matrix

# for annotating bootstraps values on the tree:
bootstrap.value <- ape::boot.phylo(
    phy = tree, 
    x = fst.matrix, 
    FUN = function(x) ape::nj(x), 
    block = 1, 
    B = 10000, 
    trees = FALSE, 
    rooted = FALSE
    )
 # to get percentage values
bootstrap.value <- round((bootstrap.value/10000)*100, 0)
bootstrap.value
# to include in the tree
tree$node.label <- bootstrap.value 


#-------------------------------------------------------------------------------#
## other type of tree ##
root.tree <- ape::as.phylo(stats::hclust(stats::dist(fst.matrix), method = "average"))
bootstrap.value <- ape::boot.phylo(
	phy = tree, 
	x = fst.matrix, 
	FUN = function(xx) ape::as.phylo(stats::hclust(stats::dist(xx), method = "average")) , 
	block = 1, 
	B = 10000, 
	trees = FALSE, 
	rooted = TRUE)
	 
bootstrap.value <- round((bootstrap.value/10000)*100, 0)
bootstrap.value
tree$node.label <- bootstrap.value

write.tree(tree, "fst.tree")
#tr = read.tree("fst.tree")

#-------------------------------------------------------------------------------#
### récuperer les couleurs
tiplab <- as.data.frame(tree$tip.label)
tiplab$color <- c("blue","orange","green","orange","red","blue","blue","blue","darkviolet",
                  "green","orange","blue","darkgreen","green")
rownames(tiplab) <- tiplab[,1]
colors <- setNames(tiplab[,2],rownames(tiplab))


svl <- dplyr::select(grape, OP, pNpSratio)
rownames(svl) <- svl$OP
svl <- setNames(svl[,2],rownames(svl))


#-------------------------------------------------------------------------------#
## USE PHYTOOLS fonction to create a tree of Fst colored according to Ds value:
obj <- contMap(tree, svl, plot = FALSE)
obj <-setMap(obj,invert=TRUE)

pdf(file = "Fig2B_Fst_tree.pdf", width = 7, height = 7) #, dpi = 600)
#par(fg="transparent")
plot(obj,fsize=0.6,ylim=c(-1,length(tree$tip.label)),
     leg.txt="piN/piS")
lastPP<-get("last_plot.phylo",env=.PlotPhyloEnv)
dev.off()

plot(obj,fsize=0.6,ylim=c(-1,length(tree$tip.label)),
     leg.txt="piN/piS")
lastPP<-get("last_plot.phylo",env=.PlotPhyloEnv)

Fig2B <- recordPlot()


#### bonus we do the same with treemix #########################################
treemix <- read.tree("03.treemix//consensus_tree.tre")
obj1 <- contMap(treemix, svl, plot = FALSE)
obj1 <-setMap(obj1,invert=TRUE)

pdf(file = "Fig2B_treemix_tree.pdf", width = 7, height = 7) #, dpi = 600)
#par(fg="transparent")
plot(obj1,fsize=0.6,ylim=c(-1,length(treemix$tip.label)),
     leg.txt="piN/piS")

lastPP<-get("last_plot.phylo",env=.PlotPhyloEnv)
dev.off()

#note; in this tree it is interesting to see the huge drift on the salmon river, but low pN/pS.
#removal of some class of deleterious mutations?

################################################################################
## we used the distRoot package from adephylo ##
d <- data.frame(distRoot(tree)) %>% set_colnames(., c("len"))
d$pop <- rownames(d)

#prepare pNps data , make a dataframe:
colnames(grape)[1] <- "pop"
svl1 <-  dplyr::select(grape, pop, pNpSratio) %>% set_colnames(., c("pop", "pnps"))

all2 <- merge(d, svl1, by ="pop")
summary(lm(all2$len ~ all2$pnps))
all2$pop

all2$Region <- c("Alaska","BC","Cascadia","BC","Califonria",rep("Alaska",2), "HaidaGwaii", "Alaska",
                 "Cascadia","BC","Thompson","Alaska","Cascadia")
all2$Region <- as.factor(all2$Region)
levels(all2$Region)

myColors <- c("blue","orange","red","green","darkviolet","springgreen4")
names(myColors) <- levels(all2$Region)
colScale <- scale_colour_manual(name = "Region",values = myColors)


plot <- ggplot(all2, aes(x = len, y = svl)) +
  geom_smooth(method = "lm", se = F) +
  geom_point(aes(colour = factor(Region)), size = 2) +
  xlab("tree branch length") + ylab(expression(pi[N]/pi[S])) +
  colScale + theme_bw() +
  th_plot +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 0.03, label.y = 0.285) 

Fig2C <- plot +  theme(legend.position = "none")

#-------------------------------------------------------------------------------#
pdf(file = "Fig2Cplot_pnps_branch.pdf", 10, 10)
Fig2C
dev.off()

#-------------------------------------------------------------------------------
###### Make a MAP with basic phytools function ################################
lat.long <- dplyr::select(grape, pop, latitude, longitude)
rownames(lat.long) <- lat.long$pop
lat.long <- lat.long[,-1]
lat.long <- as.matrix(lat.long[c(1:11,13,12,14),])

root.tree <- ape::as.phylo(stats::hclust(stats::dist(fst.matrix), method = "average"))


obj <- phylo.to.map(root.tree, as.matrix(lat.long), 
                    database="worldHires",
                    regions=c("USA","Canada"),
                    plot=FALSE, ylim = c(35,75), xlim=c(-250,-100))


#plot(obj,cex.points=c(0,1.5), colors = colors, ftype="off", rotate = FALSE, cex.points=c(0,1.5))

pdf(file="panelA_fig2.pdf", 14, 8)
plot(obj,cex.points=c(0,1.5), colors = colors, ftype="i", rotate = FALSE, cex.points=c(0,1.5))
dev.off()

plot(obj,cex.points=c(0,1.5), colors = colors, ftype="i", rotate = FALSE, cex.points=c(0,1.5))
Fig2A <- recordPlot()

#-------------------------------------------------------------------------------#
#Fig 2D:  Same with distance to the southernmost sites:
colnames(all2)
all2$dist_to_south <- all2$dist_to_ocean + all2$dist_to_source

Fig2D <- ggplot(all2, aes(x = dist_to_south, y = svl)) +
  geom_smooth(method = "lm", se = F) + #use se = to get the CI!
  geom_point(aes(colour = factor(Region)), size = 2) +
  xlab("distance to the southernmost site") + ylab(expression(pi[N]/pi[S])) +
  colScale + theme_bw() +
  th_plot +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 0.03, label.y = 0.285) 

Fig2D <- Fig2D +  theme(legend.position = "none")
Fig2D


F2A <- plot_grid(Fig2A, labels = "A") 

#this won't work. I didn't spend time figuring out why and simply merged everything together in inksapce
F2B <- plot_grid(Fig2B)
F2B

#we create fake empty plot to add panel B manually in inkscape
F2CD <- plot_grid(NULL, Fig2C, NULL, Fig2D, labels = c("","C","","D"), ncol = 2)
F2CD

pdf(file = "Fig2A_CD.pdf", 12,15)
plot_grid(F2A, F2CD, ncol = 1, rel_heights = c(1.5,1))
dev.off()

#the stuff will then be merged in inkscape



################ extra stuff below #############################################

#-------------------------------------------------------------------------------#
# Fig S8:
all2 <- merge(all2 , grape , by = "pop")

# omega NA:
plotNA <- ggplot(all2, aes(x = len, y = allsites_3sp_omegaNA)) +
  geom_smooth(method = "lm", se = F) +
  geom_point(aes(colour = factor(Region)), size = 2) +
  xlab("tree branch length") + ylab(expression(omega["NA"])) +
  labs(tag = "A") +
  colScale + theme_bw() +
  th_plot + theme(legend.position = "none") +
  theme(plot.tag = element_text(face = "bold")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 0.03, label.y = 0.285 , size = 6) 


#-------------------------------------------------------------------------------#
# omega A
plotA <- ggplot(all2, aes(x = len, y = allsites_3sp_omegaA)) +
  geom_smooth(method = "lm", se = F) +
  geom_point(aes(colour = factor(Region)), size = 2) +
  xlab("tree branch length") + ylab(expression(omega["A"])) +
  labs(tag = "B") +
  colScale + theme_bw() +
  th_plot + theme(legend.position = "none") +
  theme(plot.tag = element_text(face = "bold")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 0.03, label.y = 0.1, size = 6) 

#-------------------------------------------------------------------------------#
# alpha: 
plotAlpha <- ggplot(all2, aes(x = len, y = allsites_3sp_alpha)) +
  geom_smooth(method = "lm", se = F) +
  geom_point(aes(colour = factor(Region)), size = 2) +
  xlab("tree branch length") + ylab(expression(alpha[""])) +
  labs(tag = "C") +
  colScale + theme_bw() + th_plot + theme(plot.tag = element_text(face = "bold")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 0.10, label.y = 0.2, size = 6) 

#plotAlpha

pdf(file = "FigureS07corrigé.pdf" , 18, 6   )
ggarrange(plotNA, plotA, plotAlpha, nrow = 1, ncol = 3, common.legend =  T, legend = "bottom")
dev.off()

#-------------------------------------------------------------------------------#
# Faire les mêmes graphes avec Treemix pour les insérer dans la Figure 4 ainsi que le graph de pn/ps arbre.
dmix <- data.frame(distRoot(treemix)) %>% set_colnames(., c("treemix"))
dmix$pop <- rownames(dmix)
head(all2)
tmix <- merge(all2, dmix, by = "pop") %>% 
  select(Region, pop, treemix, pnps, 
         allsites_3sp_alpha, allsites_3sp_omegaA, allsites_3sp_omegaNA )

colnames(tmix)
tmix
rm(colScale)

myColors <- c("blue","orange","red","green","darkviolet","springgreen4")
names(myColors) <- levels(tmix$Region)
colScale <- scale_colour_manual(name = "Region",values = myColors)

plotpnps <- ggplot(tmix, aes(x = treemix, y = pnps)) +
  geom_smooth(method = "lm", se = F) +
  #geom_text(label = tmix$pop) +
  geom_point(aes(colour = factor(Region)), size = 2) +
  xlab("tree branch length") + ylab(expression(pi["N"]/pi["S"])) +
  labs(tag = "A") +
  colScale + theme_bw() +
  th_plot + theme(legend.position = "none") +
  theme(plot.tag = element_text(face = "bold")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 0.005, label.y = 0.285 , size = 6) 

plotpnps

plotNA <- ggplot(tmix, aes(x = treemix, y = allsites_3sp_omegaNA)) +
  geom_smooth(method = "lm", se = F) +
  geom_point(aes(colour = factor(Region)), size = 2) +
  xlab("tree branch length") + ylab(expression(omega["NA"])) +
  labs(tag = "B") +
  colScale + theme_bw() +
  th_plot + theme(legend.position = "none") +
  theme(plot.tag = element_text(face = "bold")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 0.005, label.y = 0.285 , size = 6) 

plotNA

plotA <- ggplot(tmix, aes(x = treemix, y = allsites_3sp_omegaA)) +
  geom_smooth(method = "lm", se = F) +
  geom_point(aes(colour = factor(Region)), size = 2) +
  xlab("tree branch length") + ylab(expression(omega["A"])) +
  labs(tag = "C") +
  colScale + theme_bw() +
  th_plot + theme(legend.position = "none") +
  theme(plot.tag = element_text(face = "bold")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 0.005, label.y = 0.075, size = 6) 

plotA

plotAlpha <- ggplot(tmix, aes(x = treemix, y = allsites_3sp_alpha)) +
  geom_smooth(method = "lm", se = F) +
  geom_point(aes(colour = factor(Region)), size = 2) +
  xlab("tree branch length") + ylab(expression(alpha[""])) +
  labs(tag = "D") +
  colScale + theme_bw() + th_plot + theme(plot.tag = element_text(face = "bold")) +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 0.02, label.y = 0.2, size = 6) 

plotAlpha

pdf(file = "FigureS04b.pdf" , 20, 6   )
ggarrange(plotpnps ,plotNA, plotA, plotAlpha, nrow = 1, ncol = 4, common.legend =  T, legend = "bottom")
dev.off()

#-------------------------------------------------------------------------------#


######################### NOW LOOK AT CORRELATION WITH BRANCH LENGTH ###########
colnames(grape)[1] <- "pop"
essai <- merge(grape, all)
colnames(essai)
summary(lm(essai$len ~ essai$allsites_3sp_alpha))
summary(lm(essai$len ~ essai$allsites_3sp_omegaA))
summary(lm(essai$len ~ essai$allsites_3sp_omegaNA))
summary(lm(essai$len ~ essai$allsites_chinook_omegaNA))

essai$region <- as.factor(essai$region)
(essai$region)

myColors <- c("blue","orange","red","green","darkviolet","springgreen4")
names(myColors) <- levels(essai$region)
colScale <- scale_colour_manual(name = "Region",values = myColors)


plotOmegaNA <- ggplot(essai, aes(x = len, y = allsites_3sp_omegaNA)) +
  geom_smooth(method = "lm", se = F) +
  geom_point(aes(colour = factor(region)), size = 2) +
  xlab("tree branch length") + ylab(expression(pi[N]/pi[S])) +
  colScale + theme_bw() +
  th_plot +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 0.03, label.y = 0.285) 

plotOmegaNA

pdf(file = "plot_omegaNA_branch_leg.pdf", 6, 5)
plot
dev.off()
write.table(essai,"grape_with_tree_branch_length.txt" , quote = F , row.names = F, sep = "\t")
