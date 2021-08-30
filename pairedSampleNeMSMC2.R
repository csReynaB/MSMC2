##################################################
## Project: cRR - MSMC2 analysis 
## Date: Sun 20 Sep 20:42:16 2020
## Author: Carlos S. Reyna-Blanco
##################################################
if (!require("data.table")) install.packages("data.table"); library("data.table")
if (!require("RColorBrewer")) install.packages("RColorBrewer"); library("RColorBrewer")

#setwd("/home/reynac/Dropbox/Backup/MSMC2/phased/v2/Results/Ne/")
setwd("C://Users/creyn/Dropbox/Backup/MSMC2/withKK1/pair/Ne/")

createTable <- function(suffix){
  
  files <- list.files(pattern=suffix)
  samples <- lapply(files, read.table, header=T)
  msmc  <- data.table()
  count <- 1
  for (name in files){
    sample <- unlist(strsplit(name, '\\.'))[1]
    model  <- unlist(strsplit(name, '\\.'))[3]
    msmc <- rbind(msmc, data.table(time_index=samples[[count]][,1], left_time_boundary=samples[[count]][ ,2], 
                                   right_time_boundary=samples[[count]][ ,3], lambda=samples[[count]][ ,4], Sample=sample, Model=model)) #sites
    count <- count + 1
  }
  return(msmc)
}

#########################
df.msmc <- createTable(suffix="final.txt")  #read files by providing the suffix
info <- read.table("samplelist", header = T) #file with samples to be plotted and extra information
cols <- list("Africa" = "gray",
             "WEurope" = "orange",
             "EAsia" = "turquoise",
             "SAmerica" = "chartreuse2",
         
             "ZagrosRegion" = "brown",
             "EasternMarmara" = c(	"#CCCC00"),
             "NorthernGreece" = c("dodgerblue"),
             "CentralSerbia" = c(brewer.pal(n=9, name = "Greens")[c(7)]),
             "Hungary-Neo" = "blue",
             "LowerAustria" = "green",
             "SouthernGermany1" = c(brewer.pal(n=9, name = "Purples")[c(5)]), 
             "SouthernGermany2" =  brewer.pal(n=9, name = "Purples")[c(7)],
             "Balkan-Fisher" = c("black"),
             "NorthernEurope-Meso" = "orange3",
             "Caucasus" = "#F0E442",
             "WesternEurope-Meso" = c("#FF69B4"),
             "DanubeGorges-Meso" =  c(brewer.pal(n=9, name = "Reds")[c(5)])
              )


period <- list(modern = info[info$Period %in% "Modern",]$Region[c(4,1,2,3)],
               neo = info[info$Period %in% "Neo",]$Region[c(8,1,6,7,5,2,3,4)],
               fisher = info[info$Period %in% "Fisher",]$Region,
               meso = info[info$Period %in% "Meso",]$Region[c(3,2,1,4)]
               )  

ltype <- list( modern =c(4,1), neo = c(1, 1.6), 
               fisher = c(3,3.4), meso = c(2,2.4) )


mu <- 1.25e-8
gen <- 29

#pdf(args[2], width=10, height=8)
#pdf("pairPops_Ne_MSMC2.pdf", width=7, height=5)

plot(df.msmc$left_time_boundary/mu*gen, (1/df.msmc$lambda)/(2*mu), log="x",ylim=c(0,35000), xlim=c(1e+04,3e+06),
     type="n", cex.axis=1.3, ylab = "", xlab="")
mtext(side=1, line=3, "Years ago", col="black", font=1,cex=1.4)
mtext(side=2, line=3, "Ne", col="black", font=1, cex=1.4)
mtext("B", side=3, line=1, adj=0, cex=2, col="black", outer=F, font = 2 )

for (n in names(period)){
  for (p in period[[n]]){
    samples <- info[which (info$Region %in% p),]$Sample
    tmp <- df.msmc[which(df.msmc$Sample %in% samples),]
    lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu),  lty=ltype[[n]][1], lwd=ltype[[n]][2], type="s", col=cols[[p]])
  }
}

#par(fig=c(0, 1, 0, 1), oma=c(2, 0, 1, 2), mar=c(3, 0, 1, 1), new=TRUE)
par(fig=c(0, 1, 0, 1), oma=c(2, 0, 3, 3), mar=c(3, 2, 1, 0), new=TRUE)

#par(fig=c(0, 1, 0, 1), oma=c(0, 0, 4, 0), mar=c(0, 0, 1, 0), new=TRUE)
#par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
plot(0, type='n', bty='n', xaxt='n', yaxt='n')
#legend("topleft", cex=1, legend=unlist(period), 
#       lty=c(rep(4,4), rep(1,8), rep(2,3), 3), col=unlist(cols),
#       lwd=2, ncol=2)

legend("bottomright", cex=1.1, legend=unlist(period), 
       lty=c(rep(4,4), rep(1,8), 3, rep(2,4) ), col=unlist(cols),
       lwd=2, ncol=1, xpd = TRUE, horiz = F, inset = c(0, 0), bty='n')

dev.off()




