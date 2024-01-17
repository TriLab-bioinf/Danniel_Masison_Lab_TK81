library(tidyverse)
library(dplyr)
library(plyr)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

filelist <- c("proteinGroups","proteinGroupsIIItech1PA","proteinGroupsIIItech2PA","proteinGroupsII_Ssa1isLSsa2isMtech1","proteinGroupsII_Ssa1isLSsa2isMtech2")

# create namekey for rename column names
namekey <- c(
        Intensity.L.1="Intensity.L.rep1",
        Intensity.L.2="Intensity.L.rep2",
        Intensity.L.3="Intensity.L.rep3",
        Intensity.L.4="Intensity.L.rep4",
        Intensity.L.5="Intensity.L.rep5",
        Intensity.H.1="Intensity.H.rep1",
        Intensity.H.2="Intensity.H.rep2",
        Intensity.H.3="Intensity.H.rep3",
        Intensity.H.4="Intensity.H.rep4",
        Intensity.H.5="Intensity.H.rep5",
        Intensity.L.1one = "Intensity.L.rep1",
        Intensity.L.2two = "Intensity.L.rep2",
        Intensity.L.3three = "Intensity.L.rep3",
        Intensity.L.4four = "Intensity.L.rep4",
        Intensity.L.5five = "Intensity.L.rep5",
        Intensity.H.1one = "Intensity.H.rep1",
        Intensity.H.2two = "Intensity.H.rep2",
        Intensity.H.3three = "Intensity.H.rep3",
        Intensity.H.4four = "Intensity.H.rep4",
        Intensity.H.5five = "Intensity.H.rep5"
        )

genelist <- list()
for(i in 1:length(filelist)){
    data <- read.delim(paste0("../data/",filelist[i],".txt"),header=T)
    data <- data[which(data$Only.identified.by.site!="+" & data$Reverse!="+" & data$Potential.contaminant!="+"),]
    # rename column names to keep consistency
    data2 <- plyr::rename(
              data, 
              replace = namekey,
              warn_missing = FALSE
    )
    data3 <- data.frame(gene=data2$Gene.names,
                    L1=data2$Intensity.L.rep1,
                    L2=data2$Intensity.L.rep2,
                    L3=data2$Intensity.L.rep3,
                    L4=data2$Intensity.L.rep4,
                    L5=data2$Intensity.L.rep5,
                    H1=data2$Intensity.H.rep1,
                    H2=data2$Intensity.H.rep2,
                    H3=data2$Intensity.H.rep3,
                    H4=data2$Intensity.H.rep4,
                    H5=data2$Intensity.H.rep5)
    data3 <- data3[!data3$gene=="",]
    # count valid value numbers
    data3$valid_value1 <- apply(data3[,2:6], 1 , function(x) sum ( x > 0 ))
    data3$valid_value2 <- apply(data3[,7:11], 1 , function(x) sum ( x > 0 ))
    # filter proteins that have less than 3 valid values
    data3 <- data3[which(data3$valid_value1 >=3 & data3$valid_value2 >=3),1:(ncol(data3)-2)]
    data3[data3 == 0] <- NA
    data3$log2FC <- log2(apply(data3[,7:11],1,mean,na.rm = TRUE)/apply(data3[,2:6],1,mean,na.rm = TRUE))
    # t.test
    for(n in 1:dim(data3)[1]){
    data3[n,13] <- t.test(data3[n,7:11],data3[n,2:6])$p.value
    }
    names(data3)[names(data3) == "V13"] <- "pvalue"
    data3$class <- filelist[i]
    data3$L.median <- apply(data3[,2:6],1,median,na.rm = TRUE)
    data3$H.median <- apply(data3[,7:11],1,median,na.rm = TRUE)
    data3$L.mean <- apply(data3[,2:6],1,mean,na.rm = TRUE)
    data3$H.mean <- apply(data3[,7:11],1,mean,na.rm = TRUE)
    data3$L.sum <- apply(data3[,2:6],1,sum,na.rm = TRUE)
    data3$H.sum <- apply(data3[,7:11],1,sum,na.rm = TRUE)
    data3$rank.L[order(-data3$L.median)] <- 1:nrow(data3)
    data3$rank.H[order(-data3$H.median)] <- 1:nrow(data3)                           
    write.table(data3,paste0(filelist[i],"_intensity_ttest.txt"),sep="\t",quote=F,row.names=F) 
    genelist[[i]] <- as.list(data3[which(abs(data3$log2FC)>1 & data3$pvalue<0.05),]$gene)
    assign(filelist[i],data3)
}

res <- rbind(proteinGroups,
proteinGroupsIIItech1PA,
proteinGroupsIIItech2PA,
proteinGroupsII_Ssa1isLSsa2isMtech1,
proteinGroupsII_Ssa1isLSsa2isMtech2
)
write.table(res,"all_intensity_ttest.txt",sep="\t",quote=F,row.names=F) 

rank_ssa1 <- pivot_wider(res[,c(1,22,14)],  names_from = "class", values_from = "rank.H",values_fn = mean)
write.table(rank_ssa1,"ssa1_rank.txt",sep="\t",quote=F,row.names=F) 

intensity_ssa1 <- pivot_wider(res[,c(1,16,14)],  names_from = "class", values_from = "H.median",values_fn = mean)
write.table(intensity_ssa1,"ssa1_intensity.txt",sep="\t",quote=F,row.names=F) 

rank_ssa2 <- pivot_wider(res[,c(1,21,14)],  names_from = "class", values_from = "rank.L",values_fn = mean)
write.table(rank_ssa2,"ssa2_rank.txt",sep="\t",quote=F,row.names=F) 

intensity_ssa2 <- pivot_wider(res[,c(1,15,14)],  names_from = "class", values_from = "L.median",values_fn = mean)
write.table(intensity_ssa2,"ssa2_intensity.txt",sep="\t",quote=F,row.names=F) 

#samples <- c("proteinGroupsIIItech1PA","proteinGroupsIIItech2PA","proteinGroupsII_Ssa1isLSsa2isMtech1","proteinGroupsII_Ssa1isLSsa2isMtech2")

# ssa1
#for(i in 1:length(samples)){
#    samplename <- paste0("proteinGroupsI_vs_",samples[i],"_ssa1")
#    df <-full_join(res[which(res$class=="proteinGroups"),c(1,7:11,16,22)],res[which(res$class==samples[i]),c(1,7:11,16,22)],by="gene")
#    row1 = dim(df)[1]
#    for(n in 1:row1){
#        df$pvalue[n] <- ifelse(apply(df[n,2:6], 1 , function(x) sum ( x > 0 ))==0 | apply(df[n,9:13], 1 , function(x) sum ( x > 0 ))==0,
#                                     "NA",
#            wilcox.test(df[n,2:6][!is.na(df[n,2:6])],df[n,9:13][!is.na(df[n,9:13])])$p.value
#    )
#    }
#    write.table(df,paste0(samplename,"_ranksum.txt"),sep="\t",row.names=F,quote=F)
#}

# ssa2
#for(i in 1:length(samples)){
#    samplename <- paste0("proteinGroupsI_vs_",samples[i],"_ssa2")
#    df <-full_join(res[which(res$class=="proteinGroups"),c(1,2:6,15,21)],res[which(res$class==samples[i]),c(1,2:6,15,21)],by="gene")
#    row1 = dim(df)[1]
#    for(n in 1:row1){
#        df$pvalue[n] <- ifelse(apply(df[n,2:6], 1 , function(x) sum ( x > 0 ))==0 | apply(df[n,9:13], 1 , function(x) sum ( x > 0 ))==0,
#                                     "NA",
#            wilcox.test(df[n,2:6][!is.na(df[n,2:6])],df[n,9:13][!is.na(df[n,9:13])])$p.value
#    )
#    }
#    write.table(df,paste0(samplename,"_ranksum.txt"),sep="\t",row.names=F,quote=F)
#}




