
library(ggplot2)
library(ggpubr)
library(patchwork)
library(showtext)
library(EPIC)
library(IOBR)
showtext.auto(enable = TRUE)
font.add('arial', 'arial.ttf') 
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 
setwd("")
file <- dir()
data <- list()
for (i in (1:length(file))) {
  data[[i]] <- read.csv(file[i],header = T)
}
name <- c("UVM (n=80)", "UCS (n=57)", "UCEC (n=545)", "THYM (n=120)", 
          "THCA (n=509)", "TGCT (n=150)", "STAD (n=415)", "SKCM (n=471)", 
          "SKCM-Primary (n=103)", "SKCM-Metastasis (n=368)", "SARC (n=260)", 
          "READ (n=166)", "PRAD (n=498)", "PCPG (n=181)", "PAAD (n=179)", 
          "OV (n=303)", "MESO (n=87)", "LUSC (n=501)", "LUAD (n=515)", 
          "LIHC (n=371)", "LGG (n=516)", "KIRP (n=290)", "KIRC (n=533)", 
          "KICH (n=66)", "HNSC (n=522)", "HNSC-HPV +  (n=98)", "HNSC-HPV- (n=422)", 
          "GBM (n=153)", "ESCA (n=185)", "DLBC (n=48)", "COAD (n=458)", 
          "CHOL (n=36)", "CESC (n=306)", "BRCA (n=1100)", "BRCA-LumB (n=219)", 
          "BRCA-LumA (n=568)", "BRCA-Her2 (n=82)", "BRCA-Basal (n=191)", 
          "BLCA (n=408)", "ACC (n=79)")

plot_tme <- function(x){
  x$pvalue = ifelse(x$adj.p >= 0.05, "p≥0.05", "p<0.05")
  x$cancer = factor(x$cancer, levels = name) 
  
  ggplot(x, aes(infiltrates, cancer, 
                shape = pvalue, 
                color = rho)) +
    geom_point(size = 3) + 
    scale_shape_manual(values = c(15, 7)) + 
    scale_color_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c") + 
    #scale_x_discrete(position = "top") + 
    theme_bw() + 
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90,hjust = 0,vjust = 0),
          axis.text.x.bottom = element_text(family = "arial"))
}

dd <- lapply(data, plot_tme)

breakpoint <- 13

up <- dd[[1]] + 
  scale_x_discrete(position = "top") + 
  guides(color=FALSE) + guides(shape=FALSE) 

for (i in 2:(breakpoint)) {
  up <- up + dd[[i]] + 
    scale_x_discrete(position = "top") + 
    theme(axis.text.y = element_blank()) + 
    guides(color = FALSE) + guides(shape = FALSE) 
}

width_up <- NULL
for (i in 1:breakpoint) {
  width_up <- c(width_up,length(unique(data[[i]]$infiltrates)))
}

up + plot_layout(guides = 'collect',widths = width_up)

down <- dd[[breakpoint + 1]] + scale_x_discrete(position = "bottom") + 
  guides(color = FALSE) + guides(shape = FALSE) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0))

for (i in (breakpoint + 2):(length(file)-1)) {
  down <- down + dd[[i]] + scale_x_discrete(position = "bottom") + 
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0)) + 
    guides(color=FALSE) + guides(shape=FALSE) 
}

down <- down + dd[[length(file)]] + scale_x_discrete(position = "bottom") + 
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0))

down <- down & scale_fill_continuous(limits = c(-1, 1))

width_down <- NULL
for (i in (breakpoint + 1):length(file)) {
  width_down <- c(width_down,length(unique(data[[i]]$infiltrates)))
}

down + plot_layout(guides = 'collect',widths = width_down)

(up + plot_layout(guides = 'collect',widths = width_up))-(down + plot_layout(guides = 'collect',widths = width_down)) + plot_layout(nrow = 2)
ggsave("TIMER.pdf", 
       height=15, width=15)