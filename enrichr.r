#https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html

# install.packages("enrichR")

library(enrichR)
library(xlsx)
library(r2excel)
library(ggplot2)
library(dplyr)
library(ggthemes)

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite("YeastEnrichr") # yeast genes   
}

# Then find the list of all available databases from Enrichr.
if (websiteLive) dbs <- listEnrichrDbs()

## if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

genelists <- list.files("/data/wangy80/TK81/perseus/genelist",pattern="*.list")

genelists

dir <- "/data/wangy80/TK81/perseus/genelist/"
files<- paste0(dir,genelists)

files

dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018","KEGG_2019")

file_name = "Functional_enrichment.xlsx"

for(i in 1:length(files)){
    samplename <- gsub("\\.list","",gsub("^.*list/","",files[i]))
    genelist <- read.delim(files[i],header=F)
    enriched <- enrichr(genelist$V1, dbs)
    res <- rbind(enriched[[1]],
                enriched[[2]],
                enriched[[3]],
                enriched[[4]])
    res2 <- res[which(res$Adjusted.P.value<0.05),]
    res3 <- res2[order(res2$Adjusted.P.value),]
    sheet_name = samplename
    # Check if an excel spreadsheet already exists, otherwise create one
    if(file.exists(file_name)){
      wb <- loadWorkbook(file = file_name)
    } else {
      wb <- createWorkbook(type="xlsx")
    }
  
    # Create new excel sheet, remove sheets if it already exits (if the sheet name is too long, the errors might occur)
    sheets <- getSheets(wb)
    if(is.element(sheet_name,names(sheets))){
      removeSheet(wb, sheetName=sheet_name)
    }
    sheet <- createSheet(wb, sheetName = sheet_name)
    xlsx.addTable(wb = wb, sheet = sheet, data = res3, startRow = 1, startCol = 1)
    # Write sorted table to Excel file as different worksheets. Need file name + Worksheet name !!!
    saveWorkbook(wb, file_name)
    assign(samplename,res3)
    # barplot
    options(repr.plot.width = 8, repr.plot.height = 6, repr.plot.res = 200)
    res3$GeneCount <- as.numeric(gsub("/.*$","",res3$Overlap))
    p<-ggplot(data=res3[1:20,], aes(x=reorder(Term,-Adjusted.P.value), y=GeneCount,fill=Adjusted.P.value)) +
      geom_bar(stat = "identity")
   
    # Horizontal bar plot
    p2 <- p + coord_flip() + scale_fill_gradient(low = "red", high = "yellow") +theme_few() + xlab('Enriched Terms') +ggtitle(samplename)
    ggsave(paste0(samplename,"_top20_barplot.png"),p2,width=12,height=5)
}




