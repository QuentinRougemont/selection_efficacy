

yName <- c('LD')
ylab <- bquote(.(yName[1]) ~ r^2)

pdf("Fig2.pdf")
read.table("Fig.BER")->EBER;
plot(EBER[,1]/1000,EBER[,2],type="l",
     col="lightblue",main="LD decay",
     xlab="Distance(Kb)",
     xlim=c(0,300),
     ylim=c(0,0.850658199974543),
     ylab=ylab, # expression(LDr^{2}),
     bty="n",
     lwd=2)
read.table("Fig.POR")->EPOR;
lines(EPOR[,1]/1000,EPOR[,2],col="deepskyblue4",lwd=2)
read.table("Fig.KWE")->EKWE;
lines(EKWE[,1]/1000,EKWE[,2],col="darkblue",lwd=2)
read.table("Fig.MSL")->EMSL;
lines(EMSL[,1]/1000,EMSL[,2],col="blue",lwd=2)
read.table("Fig.SNA")->ESNA;
lines(ESNA[,1]/1000,ESNA[,2],col="LightSkyBlue",lwd=2)
read.table("Fig.CAP")->ECAP;
lines(ECAP[,1]/1000,ECAP[,2],col="orange",lwd=2)
read.table("Fig.INC")->EINC;
lines(EINC[,1]/1000,EINC[,2],col="darkorange",lwd=2)
read.table("Fig.ROB")->EROB;
lines(EROB[,1]/1000,EROB[,2],col="darkorange3",lwd=2)
read.table("Fig.KLA")->EKLA;
lines(EKLA[,1]/1000,EKLA[,2],col="red",lwd=2)
read.table("Fig.DES")->EDES;
lines(EDES[,1]/1000,EDES[,2],col="green",lwd=2)
read.table("Fig.QUI")->EQUI;
lines(EQUI[,1]/1000,EQUI[,2],col="darkolivegreen1",lwd=2)
read.table("Fig.TSO")->ETSO;
lines(ETSO[,1]/1000,ETSO[,2],col="chartreuse",lwd=2)
read.table("Fig.PAL")->EPAL;
lines(EPAL[,1]/1000,EPAL[,2],col="darkviolet",lwd=2)
read.table("Fig.SAL")->ESAL;
lines(ESAL[,1]/1000,ESAL[,2],col="springgreen4",lwd=2)

legend("topright",
       c("Berners R","Porcupine R","Kwethluck","Mile Slough","Snake R",
          "Capilano R","Inch Cr","Robertson Cr",
         "Klamath R",
         "Deschutes R",   "Quilcene R","TsooYess",
         "Pallant",
          "Salmon R"),
       col=c("lightblue","deepskyblue4","darkblue","blue","LightSkyBlue",
             "orange","darkorange","darkorange3",
             "red",
             "green","darkolivegreen1","chartreuse",
             "darkviolet",
             "springgreen4"),
       cex=1,lty=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1),bty="n",lwd=2);
dev.off()

