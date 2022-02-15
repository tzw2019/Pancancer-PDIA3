
options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
library(xlsx)
library(tidyverse)
library(GEOquery)
library(plyr)
library(circlize)
library(ComplexHeatmap)


setwd("")
options(java.parameters = "-Xmx8000m")
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE) 

GPL96 <- getGEO("GPL96", destdir = getwd())
GPL96 <- Table(GPL96)[, c("ID", "Gene Symbol")]

pan_n <- 33
pancancertype <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", 
                   "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV",
                   "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM",
                   "UCEC", "UCS", "UVM")
diflist <- list()

for (i in 1:pan_n){
  print(paste0("Load ", pancancertype[i], " differentially expressed gene data"))
  difdata <- openxlsx::read.xlsx("", sheet = i, colNames = TRUE) %>% 
    merge(GPL96, ., by.x = "Gene Symbol", by.y = "id")
  difdata$log2fc <- as.numeric(difdata$log2fc)
  difdata <- difdata[order(difdata$log2fc, decreasing = T), ]
  diflist[[i]] <- difdata
}

names(diflist) <- pancancertype

if (!file.exists("./grpfile")){
  dir.create("./grpfile")
}
for (i in 1:pan_n){
  print(paste0(pancancertype[i], " processing"))
  tmpdata <- diflist[[i]]
  tmp_updata <- tmpdata[tmpdata$log2fc > 0, ]
  tmp_downdata <- tmpdata[tmpdata$log2fc< 0, ]
  if (nrow(tmp_updata) > 500){
    top500up <- tmp_updata[, "ID"][1:500]
  } else {top500up <- tmp_updata[, "ID"]}
  
  if (nrow(tmp_downdata) > 500){
    top500down <- rev(tmp_downdata[, "ID"])[1:500]
  } else {top500down <- rev(tmp_downdata[, "ID"])}
  write.table(top500up, paste0("grpfile/", pancancertype[i], "_up500.grp"),
              row.names = F, sep = "\t", quote = F, col.names = F)
  write.table(top500down, paste0("grpfile/", pancancertype[i], "_down500.grp"),
              row.names = F, sep = "\t", quote = F, col.names = F)
}

camp_outputlist <- list()
for (i in 1:pan_n){
  camp_output_tmp <- openxlsx::read.xlsx("00Connectivity_Map_Output.xlsx", sheet = i, colNames = T)
  camp_outputlist[[i]] <- camp_output_tmp
}
names(camp_outputlist) <- pancancertype
camp_output <- plyr::ldply(camp_outputlist, data.frame)
colnames(camp_output)[1] <- "Cancer"
str(camp_output)

oncoplotdata <- reshape2::dcast(camp_output[, c("Cancer", "cmap.name", "enrichment")], Cancer ~ cmap.name)
rownames(oncoplotdata) <- oncoplotdata$Cancer
oncoplotdata <- oncoplotdata[, -1]
mean(apply(oncoplotdata, 1, function(x) sum(!is.na(x)))) #average of 74 compounds per tumor type
oncoplotdata[, "zardaverine"]
oncoplotdata2 <- oncoplotdata[, apply(oncoplotdata, 2, function(x) sum(!is.na(x)) > 3)]
# compounds ordered by number of significantly enriched
col_sum <- apply(oncoplotdata2, 2, function(x)sum(!is.na(x)))
oncoplotdata2 <- oncoplotdata2[, order(col_sum, decreasing = F)]
calculate_sum <- function(x){
  sumdata <- data.frame(positive = sum(x > 0,na.rm = T), negative = sum(x < 0, na.rm = T))
  return(sumdata)
}
campsum <- apply(oncoplotdata2, 2, calculate_sum)
campsum <- plyr::ldply(campsum, data.frame)
rownames(campsum) <- campsum$.id
campsum <- campsum[, -1]
column_ha = HeatmapAnnotation(cancernumber = anno_barplot(campsum, axis = TRUE,
                                                          axis_param = list(side = "right"),
                                                          bar_width = 1, border = F,
                                                          gp = gpar(fill = c("#E41A1C",
                                                                             "#377EB8")),
                                                          height = unit(2.5, "cm")),
                              show_annotation_name = FALSE)
col_fun = colorRamp2(c(-1, 0, 1), c("#377EB8", "white", "#E41A1C"))

pdf("CMap_heatmap PDIA3.pdf", width = 15, height = 8)
ht_list = Heatmap(as.matrix(oncoplotdata2), col = col_fun, name = "Enrichment",
                  heatmap_width = unit(1, "npc"),
                  heatmap_height = unit(0.8, "npc"),
                  column_names_side = "top", show_row_dend = FALSE, show_column_dend = FALSE,
                  show_column_names = TRUE, show_row_names = TRUE,
                  bottom_annotation = column_ha,na_col = "white", rect_gp = gpar(col = "grey"),
                  cluster_rows =  F, cluster_columns = F, 
                  heatmap_legend_param = list(grid_width = unit(2, "cm"),
                                              grid_height = unit(1.5, "cm")))
draw(ht_list, column_title = "specific inhibitors", column_title_gp = gpar(fontsize = 25))
decorate_annotation("cancernumber", {
  grid.text("Number of cancer \ntype with p < 0.05", 
            unit(1, "npc") + unit(8, "mm"), just = "left",
            gp = gpar(fontsize = 5, col="black"))})
dev.off()