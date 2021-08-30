##################################################
## Project: MSMC analysis 
## Date: Tue April 22 16:00:44 2020
## Author: Carlos S. Reyna-Blanco
##################################################
if (!require("data.table")) install.packages("data.table"); library("data.table")
if (!require("ggplot2")) install.packages("ggplot2"); library("ggplot2")
if (!require("RColorBrewer")) install.packages("RColorBrewer"); library("RColorBrewer")

#setwd("/media/reynac/Elements/10X/MSMC/singleSampleCOutputMSMC/tmp/")
#setwd("/home/reynac/Dropbox/Backup/MSMC2/phased/v2/singleSamplesOutputMSMC/noPhasedmask/")
setwd("C://Users/creyn/Dropbox/Backup/MSMC2/phased/v2/singleSamplesOutputMSMC/noPhasedmask/")
#setwd("C://Users/creyn/Dropbox/Backup/MSMC2/noPhased/singleSamplesOutputMSMC/noPhasedmask/")



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

df.msmc <- createTable(suffix="final.txt")  
info <- read.table("samplelist2", header = T)

#region1 <- levels(factor(info$Region))[c(2,8,13,12)]
#region2 <- levels(factor(info$Region))[c(1,3,5,6)]
#region3 <- levels(factor(info$Region))[c(9,11)]
#region4 <- levels(factor(info$Region))[c(16)]
#region5 <- levels(factor(info$Region))[c(14)]
#cols <- c(brewer.pal(n = 9, name = "Blues"), brewer.pal(n=9, name = "Reds"), brewer.pal(n=9, name = "Greens"), brewer.pal(n=9, name = "Purples"),
#          brewer.pal(n=9, name = "Greys"),brewer.pal(n=11, name = "BrBG"),brewer.pal(n=11, name = "PiYG"),brewer.pal(n=8, name = "Set2"),
#          brewer.pal(n=12, name = "Paired"), brewer.pal(n=8, name = "Dark2") )
cols <- list(modern = c("orange", "pink", "black", "green", brewer.pal(n=9, name = "Greys")[5]),
             neo = c( brewer.pal(n=9, name = "Greens")[c(7)], "#CCCC00", brewer.pal(n=9, name = "Purples")[c(2)],  
                      brewer.pal(n=9, name = "Purples")[c(6)], "dodgerblue",  brewer.pal(n=9, name = "Purples")[c(9)],  "brown"),
             meso = c(brewer.pal(n=9, name = "Reds")[c(7)], "pink", "orange3"),
             fisher = "black" )
                
period <- list(modern = levels(factor(info[info$Period %in% "Modern",]$Region)),
               neo = levels(factor(info[info$Period %in% "Neo",]$Region)),
               meso = levels(factor(info[info$Period %in% "Meso",]$Region)),
               fisher = levels(factor(info[info$Period %in% "Fisher",]$Region)))


mu <- 1.25e-8
gen <- 29


#pdf("plotMSMC.pdf")
plot(df.msmc$left_time_boundary/mu*gen, (1/df.msmc$lambda)/(2*mu), log="x",ylim=c(0,35000),xlim=c(1e+04,3e+06),
     type="n", xlab="Years ago", ylab="Ne")

c <- 1
for (m in period$modern){
  s <- info[which (info$Region %in% m),]$Sample
  tmp <- df.msmc[which(df.msmc$Sample %in% s),]
  lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lty=1, type="s", col=cols$modern[c])
  c=c+1
}

c <- 1
for (n in period$neo){
  s <- info[which (info$Region %in% n),]$Sample
  tmp <- df.msmc[which(df.msmc$Sample %in% s),]
  lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=1.5, lty=2, type="s", col=cols$neo[c])
  c=c+1
}

c <- 1
for (m in period$meso){
  s <- info[which (info$Region %in% m),]$Sample
  tmp <- df.msmc[which(df.msmc$Sample %in% s),]
  lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=2.5, lty=3, type="s", col=cols$meso[c])
  c=c+1
}

c <- 1
for (f in period$fisher){
  s <- info[which (info$Region %in% f),]$Sample
  tmp <- df.msmc[which(df.msmc$Sample %in% s),]
  lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=2.5, lty=4, type="s", col=cols$fisher[c])
  c=c+1
}


legend("topleft", cex=0.6, legend=unlist(period), col=unlist(cols), lty=c(rep(1,5),rep(2,7), rep(3,3), 4), lwd=2.1, ncol=2)
#legend("topleft", cex=0.6, legend=unlist(period)[6:length(unlist(cols))], col=unlist(cols)[6:length(unlist(cols))], lty=c(rep(2,7), rep(3,3), 4), lwd=2.1)
######
######
######
cols2 <- list("Balkan-Neo" = c(brewer.pal(n=9, name = "Greens")[c(6)],brewer.pal(n=9, name = "Greens")[c(7)]),
              "EMarmara-Neo" = c(	"#CCCC00","#999900"),
              "Hungary-Neo" = brewer.pal(n=9, name = "Purples")[c(3)],
              "LAustria-Neo" = c( brewer.pal(n=9, name = "Purples")[c(8)], brewer.pal(n=9, name = "Purples")[c(9)]),
              "NGreece-Neo" = c("dodgerblue2","dodgerblue3"),
              "SGermany-Neo" = c( brewer.pal(n=9, name = "Purples")[c(4)], brewer.pal(n=9, name = "Purples")[c(5)], 
                                  brewer.pal(n=9, name = "Purples")[c(6)], brewer.pal(n=9, name = "Purples")[c(7)]) ,
              "Zagros-Neo" = "brown",
              "Balkan-Meso" =  c(brewer.pal(n=9, name = "Reds")[c(5)],brewer.pal(n=9, name = "Reds")[c(6)]),
              "NEurope-Meso" = "orange3",
              "WEurope-Meso" = c("#FF69B4","pink"),
              "Balkan-Fisher" = c("gray80","black"))


plot(df.msmc$left_time_boundary/mu*gen, (1/df.msmc$lambda)/(2*mu), log="x",ylim=c(0,35000), xlim=c(1e+04,3e+06),
     type="n", xlab="Years ago", ylab="Ne")

for (n in period$neo){
  samples <- info[which (info$Region %in% n),]$Sample
  c <- 1
  for (s in samples){
    tmp <- df.msmc[which(df.msmc$Sample %in% s),]
    lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=2, lty=1, type="s", col=cols2[[n]][c])
  c <- c+1
  }
}

for (m in period$meso){
  samples <- info[which (info$Region %in% m),]$Sample
  c <- 1
  for (s in samples){
    tmp <- df.msmc[which(df.msmc$Sample %in% s),]
    lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=2.5, lty=2, type="s", col=cols2[[m]][c])
    c <- c+1
  }
}

for (f in period$fisher){
  samples <- info[which (info$Region %in% f),]$Sample
  c <- 1
  for (s in samples){
    tmp <- df.msmc[which(df.msmc$Sample %in% s),]
    lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=2.8, lty=3, type="s", col=cols2[[f]][c])
    c <- c+1
  }
}
tmp <- df.msmc[grep("French", df.msmc$Sample),]
lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=0.8, lty=4,  type="s", col="orange")
tmp <- df.msmc[grep("Karitiana", df.msmc$Sample),]
lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=0.8, lty=4,  type="s", col="green")
tmp <- df.msmc[grep("Han", df.msmc$Sample),]
lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=0.8, lty=4,  type="s", col="pink")
tmp <- df.msmc[grep("Mende", df.msmc$Sample),]
lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu), lwd=0.8, lty=4,  type="s", col="gray")

n <- names(cols2)
labels <- c(info[info$Region %in% n[1], ]$Sample,
            info[info$Region %in% n[2], ]$Sample,
            info[info$Region %in% n[3], ]$Sample,
            info[info$Region %in% n[4], ]$Sample,
            info[info$Region %in% n[5], ]$Sample,
            info[info$Region %in% n[6], ]$Sample,
            info[info$Region %in% n[7], ]$Sample,
            info[info$Region %in% n[8], ]$Sample,
            info[info$Region %in% n[9], ]$Sample,
            info[info$Region %in% n[10], ]$Sample,
            info[info$Region %in% n[11], ]$Sample
            )

legend("topleft", cex=0.5, legend=c(labels,"French","Karitiana","Han"), 
       lty=c(rep(1,14), rep(2,5), rep(3,2), rep(4,3)), col=c(unlist(cols2),"gray","green","pink"), lwd=2, ncol=2)
dev.off()




