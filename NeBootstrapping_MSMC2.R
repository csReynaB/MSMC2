# Project: 10X - MSMC2, Ne bootstrapping analysis 
# Date: Mon 21 19:14:38 2020
# Author: Carlos S. Reyna-Blanco

if (!require("data.table")) install.packages("data.table"); library("data.table")
if (!require("RColorBrewer")) install.packages("RColorBrewer"); library("RColorBrewer")

setwd("C://Users/creyn/Dropbox/Backup/MSMC2/withKK1/pair/Ne/bootstrapping/")

createTable <- function(suffix){
  
  files <- list.files(pattern=suffix)
  samples <- lapply(files, read.table, header=T)
  msmc  <- data.table()
  count <- 1
  for (name in files){
    sample <- unlist(strsplit(name, '\\.'))[1]
    rep  <- unlist(strsplit(name, '\\.'))[4]
    msmc <- rbind(msmc, data.table(time_index=samples[[count]][,1], left_time_boundary=samples[[count]][ ,2], 
                                   right_time_boundary=samples[[count]][ ,3], lambda=samples[[count]][ ,4],
                                   Sample=sample, Rep=rep))
    count <- count + 1
  }
  return(msmc)
}



#########################
df.msmc <- createTable(suffix="final.txt")  
info <- read.table("samplelist", header = T)

period <- list(
  neo = info[info$Period %in% "Neo",]$Region,
  fisher = info[info$Period %in% "Fisher",]$Region,
  meso = info[info$Period %in% "Meso",]$Region
) 

l.clw <- list(
  "ZagrosRegion" = list("brown",c(1,1.6), 'A'),
  "NWAnatolia" = list(	"#CCCC00",c(1,1.6),'B'),
  "NorthernGreece" = list("dodgerblue",c(1,1.6),'C'),
  "CentralSerbia" = list(brewer.pal(n=9, name = "Greens")[c(7)], c(1,1.6), 'D') , 
  
  #"Hungary-Neo" = list("blue",c(1,1.6),'E'), 
  
  "LowerAustria" = list( "green",c(1,1.6),'E'),
  "SouthernGermany1" = list(brewer.pal(n=9, name = "Purples")[c(5)],c(1,1.6), 'F'), 
  "SouthernGermany2" =  list(brewer.pal(n=9, name = "Purples")[c(7)],c(1,1.6), 'G'),
  
  "Lepenski-Vir" = list("black",c(1,4), 'H'),
  
  "NorthernEurope-Meso" = list("orange3",c(2,2.5), 'I'),
  "Caucasus" = list("#F0E442", c(2,2.5), 'J'),
  "WesternEurope-Meso" = list("#FF69B4",c(2,2.5), 'K'),
  "DanubeGorges-Meso" =  list(brewer.pal(n=9, name = "Reds")[c(5)],c(2,2.5), 'L')
)



mu <- 1.25e-8
gen <- 29

pdf("supFig_bootstrappingAncientPopsNe_MSMC2.pdf", width=183*0.039370, height=183*0.039370)
par(oma = c(3, 1, 1, 1), mai=c(.6,.6,.2,.2) )
layout.matrix <- matrix(c(1:12), nrow = 4, ncol = 3)
layout(mat = layout.matrix,
       heights = c(1,1,1), # Heights of the rows
       widths = c(2,2,2)) # Widths of the  columns
#n <- "neo"
for (n in names(period)){
  #p <- "ZagrosRegion"
  for (p in period[[n]]){
    samples <- info[which (info$Region %in% p),]$Sample
    tmp <- df.msmc[which(df.msmc$Sample %in% samples),]
    
    plot(df.msmc$left_time_boundary/mu*gen, (1/df.msmc$lambda)/(2*mu), log="x",ylim=c(0,35000), xlim=c(1e+04,3e+06),
         type="n", xlab = "", ylab = "", cex.axis = 0.7)
    mtext(l.clw[[p]][[3]], side=3, line=1, adj=0, cex=0.7, col="black", outer=F, font = 2 )
    mtext(side=1, line=3, "Years ago", col="black", font=1,cex=0.7)
    mtext(side=2, line=3, paste("Ne - ", p, sep=" "), col="black", font=1, cex=0.7)
    text(5e+04, 33000, samples, cex=0.7, adj=0)
    lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu),  lty=1, lwd=2, type="s", col=brewer.pal(n=9, name = "Greys")[c(5)])

    tmp <- tmp[tmp$Rep %in% "final", ]
    lines(tmp$left_time_boundary/mu*gen, (1/tmp$lambda)/(2*mu),  lty=l.clw[[p]][[c(2,1)]], lwd=2, type="s", col=l.clw[[p]][[1]])
  }
}

dev.off()




# rCCR reading table ------------------------------------------------------

createTable1 <- function(suffix){
  
  files <- list.files(pattern=suffix)
  samples <- lapply(files, read.table, header=T)
  msmc  <- data.table()
  count <- 1
  for (name in files){
    pop1 <- strsplit(unlist(strsplit(name, '\\.'))[1], "_")[[c(1,1)]]
    pop2 <- strsplit(unlist(strsplit(name, '\\.'))[1], "_")[[c(1,2)]]
    rep  <- unlist(strsplit(name, '\\.'))[4]
    msmc <- rbind(msmc, data.table(time_index=samples[[count]][,1], 
                                   left_time_boundary=samples[[count]][ ,2], right_time_boundary=samples[[count]][ ,3], 
                                   lambda00=samples[[count]][ ,4], lambda01=samples[[count]][ ,5], lambda11=samples[[count]][ ,6], 
                                   Pop1=pop1, Pop2=pop2,Rep=rep))
    count <- count + 1
  }
  return(msmc)
}

setwd("C://Users/creyn/Dropbox/Backup/MSMC2/withKK1/multi/rCCR/bootstrapping") 

df.msmc <- createTable1(suffix="final.txt")  
info <- read.table("samplelist", header = T)

l.clw <- list(
  "ZagrosRegion" = list("brown",c(1,1.6), 'A'),
  #"NWAnatolia" = list(	"#CCCC00",c(1,1.6),'B'),
  "NorthernGreece" = list("dodgerblue",c(1,1.6),'B'),
  "CentralSerbia" = list(brewer.pal(n=9, name = "Greens")[c(7)], c(1,1.6), 'C') , 
  
  "Hungary-Neo" = list("blue",c(1,1.6),'D'), 
  
  "LowerAustria" = list( "green",c(1,1.6),'E'),
  "SouthernGermany1" = list(brewer.pal(n=9, name = "Purples")[c(5)],c(1,1.6), 'F'), 
  "SouthernGermany2" =  list(brewer.pal(n=9, name = "Purples")[c(7)],c(1,1.6), 'G'),
  
  "Lepenski-Vir" = list("black",c(1,4), 'H'),
  
  "NorthernEurope-Meso" = list("orange3",c(2,2.5), 'I'),
  "Caucasus" = list("#F0E442", c(2,2.5), 'J'),
  "WesternEurope-Meso" = list("#FF69B4",c(2,2.5), 'K'),
  "DanubeGorges-Meso" =  list(brewer.pal(n=9, name = "Reds")[c(5)],c(2,2.5), 'L')
)


pdf("supFig_bootstrapping_rCRR_MSMC2.pdf", width=183*0.039370 , height=183*0.039370)
par(oma = c(3, 1, 1, 1), mai=c(.6,.6,.2,.2) )
layout.matrix <- matrix(c(1:12), nrow = 4, ncol = 3)
layout(mat = layout.matrix,
       heights = c(1,1,1), # Heights of the rows
       widths = c(2,2,2)) # Widths of the  columns
  p <- "NWAnatolia"
    samples <- info[info$Region %in% p,]$Sample
    tmp <- df.msmc[df.msmc$Pop1 %in% samples,]
    for (i in unique(tmp$Pop2)[c(12,8,10,7,1,3,4,6,9,5,2,11)]){
      #i <- "NE1"
      pop <- as.character(info[info$Sample %in% i,]$Region)
      tmp2 <- tmp[tmp$Pop2 %in% i, ]

      plot(tmp2$left_time_boundary/mu*gen, 2 * tmp2$lambda01 / (tmp2$lambda00 + tmp2$lambda11), 
           xlim=c(1000,0.5e+05),ylim=c(0,1), log="x", type="n", xlab="", ylab="", cex.axis=0.8)
    
      mtext(l.clw[[pop]][[3]], side=3, line=1, adj=0, cex=0.6, col="black", outer=F, font = 2 )
      mtext(side=1, line=3, "Years ago", col="black", font=1,cex=0.6)
      mtext(side=2, line=3, "CCR from NW Anatolia", col="black", font=1, cex=0.6)
      abline(v=20000, lty=2,col="black")
      abline(h=0.5,lty=2,col="black")
      
      for (j in 1:20){
        #j <- 1
        tmp3 <- tmp2[tmp2$Rep %in% j,]
        lines(tmp3$left_time_boundary/mu*gen, 2 * tmp3$lambda01 / (tmp3$lambda00 + tmp3$lambda11),  
              lty=1, lwd=3, type="s", col=brewer.pal(n=9, name = "Greys")[c(5)])
      }

      tmp2 <- tmp2[tmp2$Rep %in% "final", ]
      lines(tmp2$left_time_boundary/mu*gen, 2 * tmp2$lambda01 / (tmp2$lambda00 + tmp2$lambda11),  
            lty=l.clw[[pop]][[c(2,1)]], lwd=2, type="s", col=l.clw[[pop]][[1]])
      text(1e+03, 0.8, i, cex=0.8,adj=0)
    }
    
dev.off()

