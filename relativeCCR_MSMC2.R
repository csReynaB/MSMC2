#
# Project: cRR - MSMC2 analysis 
# Date: Mon 30 Aug 18:03: 2021
# Author: Carlos S. Reyna-Blanco
#
if (!require("data.table")) install.packages("data.table"); library("data.table")
if (!require("RColorBrewer")) install.packages("RColorBrewer"); library("RColorBrewer")

#setwd("/home/reynac/Dropbox/Backup/MSMC2/phased/v2/Results/cRR/")
setwd("C://Users/creyn/Dropbox/Backup/MSMC2/withKK1/multi/rCCR/")


# Read combined MSMC2 output and estimate rCCR --------------------------------------------------
createTable <- function(suffix){
  
  files <- list.files(pattern=suffix)
  samples <- lapply(files, read.table, header=T)
  msmc  <- data.table()
  count <- 1
  for (name in files){
    pop1 <- strsplit(unlist(strsplit(name, '\\.'))[1], "_")[[c(1,1)]]
    pop2 <- strsplit(unlist(strsplit(name, '\\.'))[1], "_")[[c(1,2)]]
    model  <- unlist(strsplit(name, '\\.'))[3]
    msmc <- rbind(msmc, data.table(time_index=samples[[count]][,1], 
                                   left_time_boundary=samples[[count]][ ,2], right_time_boundary=samples[[count]][ ,3], 
                                   lambda00=samples[[count]][ ,4], lambda01=samples[[count]][ ,5], lambda11=samples[[count]][ ,6], 
                                   Pop1=pop1, Pop2=pop2, Model=model)) 
    count <- count + 1
  }
  return(msmc)
}


getCCRintersect <- function(df, val){
  
  xVec <- gen * ((df$left_time_boundary + df$right_time_boundary)/2) / mu
  yVec <- (2.0 * df$lambda01) / (df$lambda00 + df$lambda11)
  
  i <- 1
  while (yVec[i] < val){
    i <- i + 1
    if ( !(i > 1 && i <= length(yVec)) ){
      stop(paste("CCR intersection index out of bounds",i,sep = " "))
    } 
    intersectDistance <- (val - yVec[i - 1]) / (yVec[i] - yVec[i - 1])
  }
  return(xVec[i - 1] + intersectDistance * (xVec[i] - xVec[i - 1]))
}



# Reading and manipulating data -------------------------------------------
df.msmc <- createTable(suffix="combined.msmc2.final.txt")  #read files by providing the suffix
info <- read.table("samplelist", header = T) #file with samples to be plotted and extra information

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
      
             "Lepenski-Vir" = list("black",c(3,4), 'H'),
             
             "NorthernEurope-Meso" = list("orange3",c(2,2.4), 'I'),
             "Caucasus" = list("#F0E442", c(2,2.4), 'J'),
             "WesternEurope-Meso" = list("#FF69B4",c(2,2.4), 'K'),
             "DanubeGorges-Meso" =  list(brewer.pal(n=9, name = "Reds")[c(5)],c(2,2.4), 'L')
              )

mu <- 1.25e-8
gen <- 29

splitTime <- vector()
df.ccR <- data.frame()
pdf("2supFig_cRR_allPairwiseComparison_MSMC2.pdf", width=15.5, height=14.5)
par(oma = c(4, 1, 1, 1), mai=c(.6,.6,.2,.2) )
layout.matrix <- matrix(c(1:12,0,13,0,0), nrow = 4, ncol = 4)
layout(mat = layout.matrix,
      heights = c(1, 1, 1,1), # Heights of the rows
      widths = c(2,2,2,1)) # Widths of the  columns

#n <- "fisher"
for (n in names(period)){
  #p <- "Balkan-Fisher"
  for (p in period[[n]]){
    samples <- info[info$Region %in% p,]$Sample
    tmp <- df.msmc[df.msmc$Pop1 %in% samples,]
    plot(df.msmc$left_time_boundary/mu*gen, 2 * df.msmc$lambda01 / (df.msmc$lambda00 + df.msmc$lambda11), 
         xlim=c(1000,0.5e+05),ylim=c(0,1),log="x", type="n", xlab="", ylab="", cex.axis=1.4)
    
    mtext(l.clw[[p]][[3]], side=3, line=1, adj=0, cex=1, col="black", outer=F, font = 2 )
    mtext(side=1, line=3, "Years ago", col="black", font=1,cex=1.3)
    mtext(side=2, line=3, paste("CCR from",p,sep = " "), col="black", font=1, cex=1.3)
    
    abline(v=2e+04, lty=2,col="black")
    abline(h=0.5,lty=2,col="black")
    
    #i <- "AKT16-BAR25"
    for (i in unique(tmp$Pop2)){
      pop <- as.character(info[info$Sample %in% i,]$Region)
      s <- as.character(info[info$Region %in% pop,]$Sample)
      tmp2 <- tmp[tmp$Pop2 %in% i, ]
      splitTime <- getCCRintersect(tmp2, 0.5)
      df.ccR <- rbind(df.ccR, data.frame(Pop1=p, Pop2=pop, Ind1=samples, Ind2=s, DivergenceTime=splitTime))
      
      lines(tmp2$left_time_boundary/mu*gen, 2 * tmp2$lambda01 / (tmp2$lambda00 + tmp2$lambda11),  
            lty=l.clw[[pop]][[c(2,1)]], lwd=l.clw[[pop]][[c(2,2)]], type="s", col=l.clw[[pop]][[1]])
    }
    text(1e+03, 0.8, samples, cex=1.4,adj=0)
  }
}

par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
plot(0, type='n', bty='n', xaxt='n', yaxt='n')
legend("right", cex=1.5, legend=names(l.clw), 
       xpd = TRUE, horiz = F, inset = c(0, 0),
       lty=unlist(lapply(l.clw, '[[', c(2,1))), 
       col=unlist(lapply(l.clw, '[[', 1)), 
       lwd=2.3, ncol=1, bty="n")
dev.off()

save(df.ccR,file="ccR.RData")
write.table(df.ccR, file = "splitTimes_cRR.txt", sep = "\t", quote = FALSE, row.names = F)


