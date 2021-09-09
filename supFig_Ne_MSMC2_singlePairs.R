#
# Project: Ne MSMC analysis 
# Date: Mond 28 Aug 16:39:00 2021
# Author: Carlos S. Reyna-Blanco
#
if (!require("data.table")) install.packages("data.table"); library("data.table")
if (!require("RColorBrewer")) install.packages("RColorBrewer"); library("RColorBrewer")

#setwd("/home/reynac/Dropbox/Backup/MSMC2/withKK1/single/Ne/")
setwd("C://Users/creyn/Dropbox/Backup/MSMC2/withKK1/single/Ne/")


# Function to read MSMC output --------------------------------------------

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



# Single Ne data ----------------------------------------------------------

df.msmc <- createTable(suffix="final.txt")  
info <- read.table("samplelist", header = T)


period.single <- list(modern = levels(factor(info[info$Period %in% "Modern",]$Region)),
                      neo = levels(factor(info[info$Period %in% "Neo",]$Region)),
                      fisher = levels(factor(info[info$Period %in% "Fisher",]$Region)),
                      meso = levels(factor(info[info$Period %in% "Meso",]$Region))
)
cols.single <- list(
  "ZagrosRegion" = "brown",
  "EasternMarmara" = c(	"#CCCC00","#999900"),
  "NorthernGreece" = c("dodgerblue","dodgerblue2"),
  "CentralSerbia" = c(brewer.pal(n=9, name = "Greens")[c(6)],brewer.pal(n=9, name = "Greens")[c(7)]),
  "Hungary-Neo" = "blue",
  "LowerAustria" = c("green2", "green3"),
  "SouthernGermany" = c( brewer.pal(n=9, name = "Purples")[c(5)], brewer.pal(n=9, name = "Purples")[c(6)], 
                         brewer.pal(n=9, name = "Purples")[c(7)], brewer.pal(n=9, name = "Purples")[c(8)]) ,
  "Balkan-Fisher" = c(brewer.pal(n=9, name = "Greys")[c(8)],"black"),
  "NorthernEurope" = "orange3",
  "Caucasus" = "#F0E442",
  "WesternEurope-Meso"  = c("#FF69B4","pink"),
  "DanubeGorges-Meso" =  c(brewer.pal(n=9, name = "Reds")[c(5)],brewer.pal(n=9, name = "Reds")[c(6)])
)

n <- names(cols.single)
labels.single <- c(as.character(info[info$Region %in% n[1], ]$Sample),
                   as.character(info[info$Region %in% n[2], ]$Sample),
                   as.character(info[info$Region %in% n[3], ]$Sample),
                   as.character(info[info$Region %in% n[4], ]$Sample),
                   as.character(info[info$Region %in% n[5], ]$Sample),
                   as.character(info[info$Region %in% n[6], ]$Sample),
                   as.character(info[info$Region %in% n[7], ]$Sample),
                   as.character(info[info$Region %in% n[8], ]$Sample),
                   as.character(info[info$Region %in% n[9], ]$Sample),
                   as.character(info[info$Region %in% n[10], ]$Sample),
                   as.character(info[info$Region %in% n[11], ]$Sample),
                   as.character(info[info$Region %in% n[12], ]$Sample))

mu <- 1.25e-8
gen <- 29


# Plot Single Ne ----------------------------------------------------------

pdf("supFig_Ne_singlePairs_MSMC2.pdf",  width=11, height=8.5)
par(oma = c(1, 1, 2, 1),   mai=c(.7,.7,.5,.2) )# mar=c(3,1,3,1) )
layout.matrix <- matrix(c(1,2,3,4), nrow = 2, ncol = 2)
layout(mat = layout.matrix,
       heights = c(1,1), # Heights of the rows
       widths = c(2,0.6)) # Widths of the  columns

plot(df.msmc$left_time_boundary/mu*gen, (1/df.msmc$lambda)/(2*mu), log="x",ylim=c(0,35000), xlim=c(1e+04,3e+06),
     type="n",cex.axis=1.3, ylab = "", xlab="")
mtext(side=1, line=3, "Years ago", col="black", font=1,cex=1.4)
mtext(side=2, line=3, "Ne", col="black", font=1, cex=1.4)
mtext("A", side=3, line=1, adj=0, cex=2, col="black", outer=F, font = 2 )

for (n in period.single$neo){
  samples <- info[which (info$Region %in% n),]$Sample
  c <- 1
  for (s in samples){
    tmp <- df.msmc[which(df.msmc$Sample %in% s),]
    lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=2, lty=1, type="s", col=cols.single[[n]][c])
    c <- c+1
  }
}

for (m in period.single$meso){
  samples <- info[which (info$Region %in% m),]$Sample
  c <- 1
  for (s in samples){
    tmp <- df.msmc[which(df.msmc$Sample %in% s),]
    lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=2.5, lty=2, type="s", col=cols.single[[m]][c])
    c <- c+1
  }
}

for (f in period.single$fisher){
  samples <- info[which (info$Region %in% f),]$Sample
  c <- 1
  for (s in samples){
    tmp <- df.msmc[which(df.msmc$Sample %in% s),]
    lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=2.8, lty=3.4, type="s", col=cols.single[[f]][c])
    c <- c+1
  }
}

tmp <- df.msmc[grep("Mende", df.msmc$Sample),]
lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=0.8, lty=4,  type="s", col="gray")

tmp <- df.msmc[grep("French", df.msmc$Sample),]
lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=0.8, lty=4,  type="s", col="orange")

tmp <- df.msmc[grep("Han", df.msmc$Sample),]
lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=0.8, lty=4,  type="s", col="turquoise")

tmp <- df.msmc[grep("Karitiana", df.msmc$Sample),]
lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=0.8, lty=4,  type="s", col="chartreuse2")



# Pairs Ne data -----------------------------------------------------------
setwd("C://Users/creyn/Dropbox/Backup/MSMC2/withKK1/pair/Ne/")

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



# Plot Pair Ne ------------------------------------------------------------

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

# Legend panel A ----------------------------------------------------------

par(fig=c(0, 1, 0, 1), oma=c(0, 0, 3, 0), mar=c(1, 0, 1, 0), new=TRUE)
plot(0, type='n', bty='n', xaxt='n', yaxt='n')
legend("topright", cex=1.1, legend=c("Mende","French","Han","Karitiana",labels.single), 
       lty=c(rep(4,4), rep(1,14), rep(3,2), rep(2,6) ), col=c("gray","orange","turquoise","chartreuse2",unlist(cols.single)), 
       lwd=2, ncol=2, xpd = TRUE, horiz = F, inset = c(0, 0), bty='n')

# Legend for Panel B ------------------------------------------------------
par(fig=c(0, 1, 0, 1), oma=c(2, 0, 3, 3), mar=c(3, 2, 1, 0), new=TRUE)
plot(0, type='n', bty='n', xaxt='n', yaxt='n')
legend("bottomright", cex=1.1, legend=unlist(period), 
       lty=c(rep(4,4), rep(1,8), 3, rep(2,4) ), col=unlist(cols),
       lwd=2, ncol=1, xpd = TRUE, horiz = F, inset = c(0, 0), bty='n')

# END ----
dev.off()

