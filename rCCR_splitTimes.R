##################################################
## Project: split Time 
## Date: Tue 19 Oct 23:18:25 2020
## Author: Carlos S. Reyna-Blanco
##################################################
if (!require("RColorBrewer")) install.packages("RColorBrewer"); library("RColorBrewer")

#setwd("/home/reynac/Dropbox/Backup/MSMC2/phased/v2/Results/cCR/")
setwd("/Users/creyna/Documents/Dropbox/Backup/MSMC2/phased/v2/Results/cCR/")
load("ccR.RData") ##run relaticeCCR_MSMC2.R to create that R data

df.ccR$Pop1 <- factor(df.ccR$Pop1, levels=levels(factor(df.ccR$Pop1))[c(7,6,3,4,1,12,5,2,9,11,8,10)])
df.ccR$Pop2 <- factor(df.ccR$Pop2, levels=levels(factor(df.ccR$Pop2))[c(5,4,7,2,12,6,8,1,10,3,11,9)]) 
df.ccR <- df.ccR[with(df.ccR, order(Pop1,Pop2)),]

df.info <- data.frame(Pop2 = unique(df.ccR$Pop2), Shape = c( rep(19,1),1,rep(19,2),8,rep(19,2),rep(15,2), 1, 0, 19),
                      Color = c( "#9E9AC8" , "blue", "green" ,
                                 "#238B45", "black", "dodgerblue", "#CCCC00", "#FB6A4A",  "#FF69B4",  "brown", "orange3","#6A51A3" ) )
cols2 = rev(brewer.pal(11,'Spectral'))

##a)
mat <- with(df.ccR, matrix(0, ncol=length(unique(Pop2)), nrow=length(unique(Pop1))))
c=0
for (i in 1:length(unique(df.ccR$Pop1))){
  for (j in 1:length(unique(df.ccR$Pop2))){
    if(i == j){
      mat[i,j] <- NA
    } else{
      c=c+1
      mat[i,j] <- df.ccR$DivergenceTime[c]
    }
  }
}
row.names(mat) <- unique(df.ccR$Pop1)
colnames(mat) <- unique(df.ccR$Pop1)
m <- mat
rotate <- function(x) t(apply(x, 2, rev))
m <- rotate(m)
m2 <- rotate(rotate(m))


#jpeg("splitTimerCCR.jpeg", width=752, height=1000, res=90)
pdf("splitTimerCCR.pdf", height = 10, width = 7.8)
par(mar=c(1, 9, 9, 0.4) )
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),
       widths=c(3,1), heights=c(1.25,1))

image(1:ncol(m), 1:nrow(m), m, col = cols2, axes = F, xlab = "", ylab = "" , useRaster = F)
mtext("a", side=3, line=1, adj=-0.1, cex=1, col="black", outer=F, font = 2 )
mtext (at=1:ncol(m), rev(colnames(m)) , cex=0.83, side = 3, las=2, col = c(rep("black",2),"firebrick",rep("black",7),rep("firebrick",2)), line=0.15)
mtext (at=1:nrow(m), rev(rownames(m)) , cex=0.83, side = 2, las=1, col = c(rep("firebrick",2),rep("black",7), "firebrick", rep("black",2)), line=0.15)
for (x in 1:ncol(m2))
  for (y in 1:nrow(m2))
    text(x, y, round(m2[y,x],0), cex = 0.64, font=2)

###b)
df.ccR <- merge(df.ccR,df.info, by="Pop2")
df.ccR <- df.ccR[with(df.ccR, order(Pop1,Pop2)),]

par(mar=c(3.6, 11, 3, 2)  ) 
plot(df.ccR$DivergenceTime, rev(df.ccR$Pop1), log="x",  type="n",
     xlab="", ylab="", cex.axis=0.9, yaxt="n", xaxt="n")
abline(v=c(20000,30000,40000, 50000),lty=2,col="black")
mtext("b", side=3, line=1, adj=-0.24, cex=1, col="black", outer=F, font = 2 )
points(df.ccR$DivergenceTime, rev(df.ccR$Pop1),  type="p", pch=df.ccR$Shape, col=as.character(df.ccR$Color) )
axis(2, at=12:1, labels=rep("",12),col = NA, col.ticks = 1,lwd.ticks = 0.9, tck=-0.01, line=0)
axis(1, c(20000,30000,40000, 50000), col = NA, col.ticks = 1,lwd.ticks = 0.9, tck=-0.01, line=-0.8)
mtext(side=2, line=0.6, at=12:1, unique(df.ccR$Pop1) , 
      col=c(rep("black",2),"firebrick",rep("black",7),rep("firebrick",2)), font=1,cex=0.9, las=2)
mtext(side=2, line=9.4, "Population" , col="black", font=2, cex=1.2, las=3)
mtext(side=1, line=2, "Split time estimated at CCR of 0.5 (Years ago)", col="black", font=2,cex=1.2)

par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(4.3, 0, 1, 0), new=TRUE)
plot(0, type='n', bty='n', xaxt='n', yaxt='n', xlab = "")
legend("bottomright", cex=0.9, legend=unique(df.ccR$Pop1), y.intersp = 0.7, 
       pch=c( rep(19,2),1,rep(19,2),8,rep(19,2),rep(15,2), 1, 0), 
       col= c( "#6A51A3", "#9E9AC8" , "blue", "green" ,
               "#238B45", "black", "dodgerblue", "#CCCC00", "#FB6A4A",  "#FF69B4",  "brown", "orange3" ),
       text.col = c(rep("black",2),"firebrick",rep("black",7),rep("firebrick",2)),
       ncol=1, xpd = TRUE, horiz = F, inset = c(0, 0), bty='n' )

dev.off()
