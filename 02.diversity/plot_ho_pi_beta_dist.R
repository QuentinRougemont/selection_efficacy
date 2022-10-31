library(ggplot2)
library(dplyr)
library(data.table)
library(magrittr)

region <- read.table("strata3.txt.gz",T)

hierf <- fread("zcat hierfstat.ho.gz ")
hierfs <- fread("zcat hierfstat.hs.gz ")

h   <-  data.frame((colMeans(hierf[,-1], na.rm=T))) %>% set_colnames(.,c("ho"))
hs  <-  data.frame((colMeans(hierfs[,-1], na.rm=T))) %>% set_colnames(.,c("hs"))
h$POP <-rownames(h)
hs$POP <-rownames(hs)

beta <- read.table("beta_iovl_ci2.gz",T)[,c(1:2)]
hsb <- merge(beta, hs, by="POP")
hsbd <- merge(hsb, region, by="POP")
hsbd <- merge(hsbd, h, by="POP")
hsbd$dist_tot = hsbd$dist_max_km + hsbd$dist_SCO_QR

summary(lm(hsbd$hs ~ hsbd$dist_source))
summary(lm(hsbd$ho ~ hsbd$dist_source))
summary(lm(hsbd$b ~ hsbd$dist_source))

summary(lm(hsbd$hs ~ hsbd$dist_tot))
summary(lm(hsbd$ho ~ hsbd$dist_tot))
summary(lm(hsbd$b ~ hsbd$dist_tot))

# sans thompson
summary(lm(hsbd[-12,16] ~ hsbd[-12,17]))
#hs:
summary(lm(hsbd[-12,3] ~ hsbd[-12,17]))
#bst
summary(lm(hsbd[-12,2] ~ hsbd[-12,17]))

#div_dist$region<-factor(div_dist$region)
#faire le plot

#general theme:
th <- theme(axis.title.x=element_text(size=14, family="Helvetica",face="bold"),
      axis.text.x=element_text(size=14,family="Helvetica",face="bold", angle=90, hjust=0, vjust=0.5),
      axis.title.y=element_text(size=20, family="Helvetica",face="bold",angle=0, hjust=0, vjust=0.5),
      axis.text.y=element_text(size=14,family="Helvetica",face="bold"),
      strip.text.x = element_text(size=18),
      panel.grid.major = element_blank() )#, panel.grid.minor = element_blank())


myColors <- c("blue","orange","red","green","darkviolet","springgreen4")
hsbd$region<-as.factor(hsbd$region)
names(myColors) <- levels(hsbd$region)
colScale <- scale_colour_manual(name = "region",values = myColors)

p <- ggplot(data=hsbd,aes(x=dist_tot, y=hs)) +
  stat_smooth(method="lm", fill="gray87") + 
  theme_bw() + 
  geom_vline(colour="darkgrey",linetype=5, size=3,  data = NULL, xintercept = 1750) + #approximate location of the ice sheet + 
  geom_point(aes(colour=factor(region)),size=4) + 
  colScale
p <- p + labs(x="Distance from southernmost site (km)", y = "Hs")  + 
	expand_limits(y=0.005, x=8000) + 
  scale_x_continuous(limits = c(0, 8000))
p <- p + th

#ajouter les pval
p1 <- p+ annotate(geom="text", 
      x=2800, y=0.028, 
      label= "paste(italic(R) ^ 2, \" = .73 p <0.0001***, slope = -7.56e-06\")", 
      parse = TRUE,
      color="#CCCCCC", size=5)

pdf(file="Ho_dist.pdf",10,8)
p1
dev.off()

p1 <- p1 + annotate(geom="text", 
      x=3300, y=0.017, 
      label= "paste(italic(R) ^ 2, \" without Thompson = .87 p < 0.0001***, slope = -7.9e-06\")", 
      parse = TRUE,
      color="#CCCCCC", size=5)

pdf(file="Ho_dist_correlation_thompson.pdf",10,8)
p1
dev.off()

###same graphe with Bst:
p <- ggplot(data=hsbd,aes(x=dist_tot, y=beta)) + 
  stat_smooth(method="lm", fill="gray86") + 
  theme_bw() + 
  geom_point(aes(colour=factor(region)),size=4) + 
  colScale + 
  stat_smooth(method="lm") + 
  labs(x="Distance from southernmost site (km)", y = expression(beta))  + 
   expand_limits(y=0.005, x=8000) + 
   scale_x_continuous(limits = c(0, 8000)) + 
  th
#p  <- p + theme(legend.position="none")

#ajouter les pval
p1 <- p+ annotate(geom="text", 
                  x=1800, y=0.5, 
                  label= "paste(italic(R) ^ 2, \" = .74 p <0.0001***, slope = 7.05e-05\")", 
                  parse = TRUE,
                  color="black", size=5)
pdf(file="Bst_dist.pdf",10,8)
p1
dev.off()

p1 <- p1 + annotate(geom="text", 
                    x=2500, y=0.4, 
                    label= "paste(italic(R) ^ 2, \" without Thompson = .87 p < 0.0001***, slope = 7.4e-05\")", 
                    parse = TRUE,
                    color="black", size=5)

pdf(file="Bst_dist_correlation_thompson.pdf",10,8)
p1
dev.off()
