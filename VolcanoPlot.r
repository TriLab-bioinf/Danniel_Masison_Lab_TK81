# load libraries
library(ggplot2)
library(ggrepel)
library(org.Sc.sgd.db)
library(clusterProfiler)
library(biomaRt)
library(EnhancedVolcano)
library(ggthemes)
library(enrichR)

filelist <- c("proteinGroups_perseus","proteinGroupsIIItech1PA_perseus","proteinGroupsIIItech2PA_perseus","proteinGroupsII_Ssa1isLSsa2isMtech1_perseus","proteinGroupsII_Ssa1isLSsa2isMtech2_perseus")

# volcanoplot
volcano <- function(data=data,prefix=filename,sig_genes=sig_genes){
    options(repr.plot.width = 5.5, repr.plot.height = 4, repr.plot.res = 300)
    p <-ggplot(data,aes(x = T.test.Difference, y = X.Log.T.test.p.value,color=color)) + 
      geom_point(size=1) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
      geom_vline(xintercept = c(log2(0.5), log2(2)),linetype = "dashed") +
      scale_color_manual(values = c("significant" = "red", "NS" = "grey")) +
      xlab("log2FoldChange") + ylab("-log10(p-value)")
    p2 <- p + geom_text_repel(data = sig_genes, aes(label = Gene.names),size = 3,max.overlaps = Inf) +
    theme_few()
    ggsave(paste0(prefix,"_volcano.pdf"),p2,height=4,width=5.5)
}

volcanoplot2 <- function(data=data,prefix=filename){
    options(repr.plot.width = 4, repr.plot.height = 4, repr.plot.res = 300)
    p <- EnhancedVolcano(data,
    lab = data$Gene.names,
    x = 'T.test.Difference',
    y = 'pvalue',
    pCutoff = 0.05,
    FCcutoff = 1)
    ggsave(paste0(prefix,"_volcano2.pdf"),p,height=8,width=7)
}

# GO enrichment
GO_enrich <- function(annotLookup=annotLookup,prefix=prefix){
    enrich_GO <- enrichGO(gene = annotLookup$ensembl_gene_id, 
                                   ont           = "ALL",
                                   OrgDb         = org.Sc.sgd.db,
                                   pAdjustMethod = "BH",
                                   keyType       = "ENSEMBL",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05)

    write.table(enrich_GO,paste0(prefix,"_GO_enrichment.txt"),sep="\t",quote=F)
    if(dim(enrich_GO)[1]>2){
    options(repr.plot.width = 7, repr.plot.height = 8, repr.plot.res = 300)
    p1 <- dotplot(enrich_GO,label_format = 100 ,showCategory=20, font.size=8)
    ggsave(paste0(prefix,"_GO_enrichment.pdf"),p1,height=8,width=7)}
}

# GO enrichment enrichr
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("YeastEnrichr") # Human genes   
}
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018","KEGG_2019")

GO_enrich2 <- function(genelist=genelist,prefix=prefix){
    enriched <- enrichr(genelist, dbs)
    res <- rbind(enriched[[1]],
                enriched[[2]],
                enriched[[3]],
                enriched[[4]])
    res2 <- res[which(res$Adjusted.P.value<0.05),]
    res3 <- res2[order(res2$Adjusted.P.value),]
    write.table(res3,paste0(prefix,"_enrichr.txt"),sep="\t",quote=F)
    options(repr.plot.width = 8, repr.plot.height = 6, repr.plot.res = 200)
    res3$GeneCount <- as.numeric(gsub("/.*$","",res3$Overlap))
    p<-ggplot(data=res3[1:20,], aes(x=reorder(Term,-Adjusted.P.value), y=GeneCount,fill=Adjusted.P.value)) +
      geom_bar(stat = "identity")
   
    # Horizontal bar plot
    p2 <- p + coord_flip() + scale_fill_gradient(low = "red", high = "yellow") +theme_few() + xlab('Enriched Terms')
    ggsave(paste0(prefix,"_top20_barplot.png"),p2,width=12,height=5)

}

# load data
for(i in 1:length(filelist)){
    data <- read.delim(paste0(filelist[i],".txt"),header=T)
    data <- data[2:dim(data)[1],]
    data$T.test.q.value <- as.numeric(data$T.test.q.value)
    data$T.test.Difference <- as.numeric(data$T.test.Difference)
    data$X.Log.T.test.p.value <- as.numeric(data$X.Log.T.test.p.value)
    data$pvalue <- 10^(-data$X.Log.T.test.p.value)
    data$color <- ifelse(abs(data$T.test.Difference)>1 & data$pvalue<0.05,"significant","NS")
    sig_genes <- data[which(data$pvalue <0.05 & (data$T.test.Difference < -1 | data$T.test.Difference>1)),]
    unique_symbols <- unique(unlist(strsplit(sig_genes$Gene.names, split = ";"), recursive = FALSE))
    mart <- useMart("ENSEMBL_MART_ENSEMBL")
    mart <- useDataset("scerevisiae_gene_ensembl", mart)
    annotLookup <- getBM(
      mart=mart,
      attributes=c("ensembl_gene_id", "gene_biotype", "external_gene_name"),
      filter="external_gene_name",
      values=unique_symbols,
      uniqueRows=TRUE)
    volcano(data=data,prefix=filelist[i],sig_genes=sig_genes)
    volcanoplot2(data=data,prefix=filelist[i])
    GO_enrich(annotLookup=annotLookup,prefix=filelist[i])
    GO_enrich2(genelist=unique_symbols,prefix=filelist[i])
}



