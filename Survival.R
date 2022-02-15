library(data.table)
library(survival)
library(pheatmap)
library(survminer)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

setwd("")
panexpr <- fread("tcga_RSEM_gene_tpm",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
annopb <- read.delim("gencode.v23.annotation.gene.probemap",row.names = 1,check.names = F,stringsAsFactors = F,header = T,sep = "\t")

panexpr <- as.data.frame(panexpr)
rownames(panexpr) <- panexpr$sample; panexpr <- panexpr[,-1]
comgene <- intersect(rownames(annopb), rownames(panexpr))
panexpr <- panexpr[comgene,]; annopb <- annopb[comgene,]
panexpr$genename <- annopb$gene; panexpr <- panexpr[!duplicated(panexpr$genename),]
rownames(panexpr) <- panexpr$genename; panexpr <- panexpr[,-ncol(panexpr)]
panexpr[1:3,1:3]

gene_group <- read.table("easy_input_gene.txt",sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)

panexpr <- panexpr[gene_group$Symbol,]
gc() 

panexpr <- 2^panexpr - 0.001 
panexpr[panexpr < 0] <- 0
panexpr <- log2(panexpr + 1) 

panexpr <- panexpr[,apply(panexpr, 1, sd) > 0] 

pansurv <- read.table("pancancerPFIData.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
rownames(pansurv) <- paste0(rownames(pansurv),"-01")

tumors <- unique(pansurv$type) 
comsam <- intersect(colnames(panexpr),rownames(pansurv)) 
pansurv <- pansurv[comsam,]
panexpr <- panexpr[,comsam]
tumors <- unique(pansurv$type)

tumors <- unique(pansurv$type)

survland.coxhr <- matrix(NA,nrow = nrow(gene_group),ncol = length(tumors),dimnames = list(gene_group$Symbol,tumors)) 
survland.coxp <- matrix(NA,nrow = nrow(gene_group),ncol = length(tumors),dimnames = list(gene_group$Symbol,tumors)) 
survland.coxplot <- matrix(0,nrow = nrow(gene_group),ncol = length(tumors),dimnames = list(gene_group$Symbol,tumors)) 

survland.logrankhr <- matrix(NA,nrow = nrow(gene_group),ncol = length(tumors),dimnames = list(gene_group$Symbol,tumors)) 
survland.logrankp <- matrix(NA,nrow = nrow(gene_group),ncol = length(tumors),dimnames = list(gene_group$Symbol,tumors)) 
survland.logrankplot <- matrix(0,nrow = nrow(gene_group),ncol = length(tumors),dimnames = list(gene_group$Symbol,tumors)) 


for(t in tumors) {
  for (g in gene_group$Symbol) { 
    sam <- rownames(pansurv[which(pansurv$type == t),]) 
    expr <- as.numeric(panexpr[g,sam]) 
    
    expr.surv <- data.frame(futime = pansurv[sam,"PFI.time"], 
                            fustat = pansurv[sam,"PFI"],
                            expr = expr, 
                            stringsAsFactors = F)
    
    cox <- coxph(Surv(futime,fustat) ~ expr, data = expr.surv) 
    coxSummary <- summary(cox)
    hr <- as.numeric(coxSummary$coefficients[,"exp(coef)"])[1]
    pvalue <- as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1] 
    survland.coxhr[g,t] <- hr
    survland.coxp[g,t] <- pvalue

    if(pvalue < 0.05) { 
      survland.coxplot[g,t] <- ifelse(hr > 1, 1, -1)
    }

    res.cut<- surv_cutpoint(expr.surv, time="futime", event="fustat",variables = c("expr"))
    expr.surv$group = ifelse(expr > res.cut$cutpoint$cutpoint,"high","low")
    expr.surv$group <- factor(expr.surv$group, levels = c("low", "high"))
    
    data.survdiff <- survdiff(Surv(futime,fustat) ~ group, data = expr.surv)
    pvalue <- 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    hr <- (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
    survland.logrankhr[g,t] <- hr
    survland.logrankp[g,t] <- pvalue
    
    if(pvalue < 0.05) {
      survland.logrankplot[g,t] <- ifelse(hr > 1, 1, -1)
    }
  }
}


write.table(survland.coxplot, file = "cox_genes associated with the OS.txt", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(survland.coxhr,file = "cox HR in survival landscape.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(survland.coxp,file = "cox pvalue in survival landscape.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

write.table(survland.logrankplot, file = "logrank_genes associated with the OS.txt", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(survland.logrankhr,file = "logrank HR in survival landscape.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(survland.logrankp,file = "logrank pvalue in survival landscape.txt",sep = "\t",row.names = T,col.names = NA,quote = F)