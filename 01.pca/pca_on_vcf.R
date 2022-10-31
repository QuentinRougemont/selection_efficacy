
#DATE: 10-05-18 
#PURPOSE: script to perform PCA on vcffiles
#AUTHOR: Q. Rougemont
#INPUT: vcf file (compressed or not) , strata file

## common checks
if("dplyr" %in% rownames(installed.packages()) == FALSE)
{install.packages("dplyr", repos="https://cloud.r-project.org") }
if("vcfR" %in% rownames(installed.packages()) == FALSE)
{install.packages("vcfR", repos="https://cloud.r-project.org") }
if("ade4" %in% rownames(installed.packages()) == FALSE)
{install.packages("ade4", repos="https://cloud.r-project.org") }
if("adegenet" %in% rownames(installed.packages()) == FALSE)
{install.packages("adegenet", repos="https://cloud.r-project.org") }
if("factoextra" %in% rownames(installed.packages()) == FALSE)
{install.packages("factoextra", repos="https://cloud.r-project.org") }

## load libs
libs <- c('dplyr','vcfR','ade4','adegenet', 
    'factoextra','vegan','ggsci','magrittr')
invisible(lapply(libs, library, character.only = TRUE))

## load and transform data
vcf<-read.vcfR(vcf, verbose=F)                                                                                                                                                                                     
genpop <- vcfR2genind(vcf)   
X <- scaleGen(genpop, NA.method=c("mean")) #imputation                       

pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=100)

## singificance of axis
eig.val <- get_eigenvalue(pca1) #first we get the eigen value in an easy form
eig.val$eig <- eig.val$variance.percent/100 #percent
expected <- bstick(length(eig.val$eig) ) 
signif <- eig.val$eig > expected #get signifcicant axis
#using this we can choose a number of axis1 
signif[1:20]
# 
# ## look at eigen vals:                                                       
pdf(file="eigen_value.pca.salmon.pdf")                                       
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))          
dev.off() 

############################################################################ 
##the most basic and awfull plot                                             
pdf(file="pca.li.pdf")                                                       
s.label(pca1$li,xax=1, yax=2)                                                
title("PCA full salmon\naxes 1-2")                                           
add.scatter.eig(pca1$eig[1:20], 3,1,2)                                       
dev.off()                                                                    
####### Save a few stuff of interest ########################################
#sauvegarde des %d'inertie capturée etc                                      
inertie<-inertia.dudi(pca1, row.inertia = TRUE, col.inertia = TRUE)          
sink("pca.inertie.25.10.txt") 
print(inertie)                
sink()                        
# 
res.var <- get_pca_var(pca1)  
sink("coord.txt")             
print(res.var$coord)          
sink()                        
# 
res.ind <- get_pca_ind(pca1)  
sink("res.ind.coord.txt")         
print(res.ind$coord)          
sink()                        
sink("res.ind.cos2.txt")         
print(res.ind$cos2)          
sink()  
                       
sink("res.ind.contrib.txt")         
print(res.ind$contrib)          
sink()  


# 
#eigen value to file:         
sink("res.eigen.txt")         
print(eig.val)                
sink()                        
# 
######## Prepare data for PLOT #########################################
#### NICE PLOT ###############
strata <- read.table("strata.txt",T) #pop level strata

#getting the strata:
library(tidyr)
ind  <- data.frame(rownames(pca1$li)) %>% set_colnames(.,'ind')
ind1 <- ind%>% separate(ind,c("POP","ID") ,sep="_")
ind  <- cbind(ind,ind1)
#we merge to have regional ID, importantly, we don't sort to preserve individual orders:
strata <- merge(ind,strata,by="POP", sort=F)   

p <- fviz_pca_ind(pca1, label="none", pointsize = 0.0) +
      geom_text(aes(label=strata$POP, colour=factor(strata$POP)),
      size = 2 )
p <- p + scale_color_igv() + theme_minimal() 
pdf(file="pca_from_vcffile_12.pdf")
p
dev.off()
#axe 3-4

p <- fviz_pca_ind(pca1, axes = c(3,4), label="none", pointsize = 0.0) +
      geom_text(aes(label=strata$POP, colour=factor(strata$POP)),
      size = 2 )
p <- p + scale_color_igv() + theme_minimal() 
pdf(file="pca_from_vcffile_axe34.pdf")
p
dev.off()

p <- fviz_pca_ind(pca1, axes = c(5,6), label="none", pointsize = 0.0) +
      geom_text(aes(label=strata$POP, colour=factor(strata$POP)),
      size = 2 )
p <- p + scale_color_igv() + theme_minimal() 
pdf(file="pca_from_vcffile_axe56.pdf")
p
dev.off()

########
col <- c("blue","orange","red","green","darkviolet","darkgreen")

#now plot
pca12 <- fviz_pca_ind(pca1,
             geom.ind = "point", # show points only (nbut not "text")
             pointshape = 19,
             col.ind = strata$region, 
             palette=col,
             addEllipses = TRUE, # Concentration ellipses
             #ellipse.type = "confidence",
             legend.title = "Region",
             #axes = c(2, 3))
             )

pdf(file="pca12_august.pdf",7.5,7.5)
pca12
dev.off()
             


### Plot without a prior
pdf(file="pca.full.1_2.pdf", 12,12)      
colorplot(pca1$li[c(1,2)],               
        pca1$li,                         
        transp=TRUE,                     
        cex=3,                           
        xlab="PC 1",                     
        ylab="PC 2")                     
abline(v=0,h=0,col="grey", lty=2)        
dev.off()  

pdf(file="pca.full.3_4.pdf", 12,12)    
colorplot(pca1$li[c(3,4)],               
        pca1$li,                         
        transp=TRUE,                     
        cex=3,                           
        xlab="PC 3",                     
        ylab="PC 4")                     
abline(v=0,h=0,col="grey", lty=2)        
dev.off()                              
pdf(file="pca.full.5_6.pdf", 12,12)    
colorplot(pca1$li[c(5,6)],               
        pca1$li,                         
        transp=TRUE,                     
        cex=3,                           
        xlab="PC 5",                     
        ylab="PC 6")                     
abline(v=0,h=0,col="grey", lty=2)        
dev.off() 


#### Contribution des SNPs; utiles pour Charles:
# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
