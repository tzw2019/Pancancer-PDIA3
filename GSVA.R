options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
setwd("")
library(clusterProfiler)
library(limma)
library(ggplot2)
library(data.table)
library(ggpubr)
library(GSVA)
source("twoclasslimma.R")
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE)

rawAnno <- read.delim("merged_sample_quality_annotations.tsv",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T) 
rawAnno$simple_barcode <- substr(rawAnno$aliquot_barcode,1,15)
samAnno <- rawAnno[!duplicated(rawAnno$simple_barcode),c("cancer type", "simple_barcode")]
samAnno <- samAnno[which(samAnno$`cancer type` != ""),]

write.table(samAnno,"easy_input_sample_annotation.txt",sep = "\t",row.names = F,col.names = T,quote = F)

samAnno <- read.table("easy_input_sample_annotation.txt", sep = "\t",header = T, check.names = F)

expr <- fread("EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",sep = "\t",stringsAsFactors = F,check.names = F,header = T)
expr[1:3,1:3]


expr <- as.data.frame(expr); rownames(expr) <- expr[,1]; expr <- expr[,-1]
gene <- sapply(strsplit(rownames(expr),"|",fixed = T), "[",1)
expr$gene <- gene
expr <- expr[!duplicated(expr$gene),]
rownames(expr) <- expr$gene; expr <- expr[,-ncol(expr)]
expr[expr < 0] <- 0 
colnames(expr) <- substr(colnames(expr),1,15)
gc()

es <- log2(expr[rownames(expr) == "",] + 1) 

msigdb.hallmark <- read.gmt("h.all.v7.2.symbols.gmt") 

pct <- 0.3 

gseaTab <- NULL

tumors <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC",
            "ESCA","GBM","HNSC","KICH","KIRC",
            "KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD",
            "PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM") 

for (i in tumors) {
  message("--",i,"...")
  sam <- samAnno[which(samAnno$`cancer type` == i),"simple_barcode"]
  comsam <- intersect(colnames(es), sam) 
  
  tumsam <- comsam[substr(comsam,14,14) == "0"] 
  # tumsam <- comsam 
  es_subset <- as.data.frame(t(es[,tumsam]))
  es_subset <- es_subset[order(es_subset$CDKN2A,decreasing = T),,drop = F]
  
  high.es <- rownames(es_subset[1:(nrow(es_subset) * pct),,drop = F])
  low.es <- rownames(es_subset[nrow(es_subset):(nrow(es_subset) - nrow(es_subset) * pct + 1),,drop = F]) 
 
  subt <- data.frame(condition = rep(c("high","low"),c(length(high.es),length(low.es))),
                     row.names = c(high.es, low.es),
                     stringsAsFactors = F)
  gset <- log2(na.omit(expr[,rownames(subt)]) + 1)
  twoclasslimma(subtype  = subt, 
                featmat  = gset, 
                treatVar = "high",
                ctrlVar  = "low", 
                prefix   = paste0("TCGA_enrichment_score_",i),
                overwt   = T, 
                sort.p   = F, 
                verbose  = TRUE, 
                res.path = ".") 

  res <- read.table(paste0("TCGA_enrichment_score_",i,"_limma_test_result.high_vs_low.txt"),sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

  res <- res[order(res$log2fc, decreasing = T),]
  glist <- res$log2fc; names(glist) <- rownames(res)

  set.seed(20211114)
  gsea <- GSEA(geneList = glist,
               pvalueCutoff = 1, 
               seed = TRUE,
               TERM2GENE = msigdb.hallmark)
  gc()
  gsea.df <- as.data.frame(gsea) 
  write.table(gsea.df,file = "output_gsea between high and low group of enrichment score in pancancer.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
  
  gseaTab <- rbind.data.frame(gseaTab,
                              data.frame(term = gsea.df$ID,
                                         NES = gsea.df$NES,
                                         FDR = gsea.df$p.adjust,
                                         tumor = i,
                                         stringsAsFactors = F),
                              stringsAsFactors = F)
}

write.table(gseaTab, "output_TCGA_pancan_gsea_regarding_es_group.txt",sep = "\t",row.names = F,col.names = T,quote = F)



darkblue <- "#303B7F"
darkred <- "#D51113"
yellow <- "#EECA1F"

tmp <- gseaTab
tmp$term <- gsub("HALLMARK_","",tmp$term)
my_palette <- colorRampPalette(c(darkblue,yellow,darkred), alpha=TRUE)(n=128)
ggplot(tmp, aes(x=tumor,y=term)) +
  geom_point(aes(size=-log10(FDR),color=NES)) +
  scale_color_gradientn('NES', 
                        colors=my_palette) + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, size = 12, hjust = 0.3, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,1,1), "lines"))
ggsave("GSEA regarding CDKN2A group in pancancer.pdf", width = 12,height = 15)