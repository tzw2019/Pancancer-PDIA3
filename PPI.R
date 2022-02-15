inputFile=""
setwd("")
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
idata <- read.table(file ="comppi--interactors_of_CDKN2A.txt", sep = "\t", header = T, quote = "")
idata[1:4,]
ilinks <- NULL
for(i in 1:nrow(idata)){
  for(j in strsplit(idata$Major.Loc.With.Loc.Score[i], split = "\\|")[[1]]){
    tmp <- strsplit(j, split = "\\(")[[1]]
    ilinks <- rbind(ilinks, c(idata$Interactor[i], 
                              tmp[1],
                              gsub(pattern = "\\)", "", tmp[2])))
  }
}
colnames(ilinks) <- c("Interactor", "Loc","LocScore")
ilinks <- as.data.frame(ilinks)
head(ilinks)
write.table(ilinks,file="PDIA3-PPI.xls",sep="\t",row.names=T,quote=F)



library(ggplot2)
library(reshape2)
library(igraph)
library(magrittr)
library(clusterProfiler)
sourceGene <- ""
sourceName <- ""
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, lheight = 30, add=TRUE, inches=FALSE)
         })
}

add.vertex.shape("fcircle", clip=igraph.shape.noclip,
                 plot=mycircle, parameters=list(vertex.frame.color=1,
                                                vertex.frame.width=1))
inputFile=""
links=read.table(inputFile,sep="\t",header=T,check.names=F)
head(links)
net <- graph_from_data_frame(links, directed = F)

pal_loc <- setNames(RColorBrewer::brewer.pal(8, name = "Set3")[-2], 
                    sort(unique(links$Loc)))

# see plot.igraph for more settings 
vertexes <- get.vertex.attribute(net, "name")
net %<>% 
  set.vertex.attribute("label", 
                       value = ifelse(vertexes == sourceGene, 
                                      sourceName, vertexes)) %<>%
  set.edge.attribute("color", 
                     value = pal_loc[links$Loc]) %<>%
  set.edge.attribute("weight", 
                     value = links$score)

set.seed(777)
lay <- layout_with_dh(net) #weight.node.dist = unique(links[,c("target","score")])$score) # weight here may be invalid!

# plot in SVG device, then you may use Adobe Illustrator 
# to fine tune the network map. 
# Such as adjust the text position to avoid overlapping.

pdf(file = "interaction.pdf", width = 6, height = 6)
plot(net, layout = lay, 
     vertex.shape = 'fcircle',
     vertex.size = ifelse(vertexes == sourceGene, 15, 3),
     vertex.color = "hotpink",
     vertex.frame.color = ggplot2::alpha("white", alpha = 0.5),
     vertex.frame.width = ifelse(vertexes == sourceGene, 30, 0.1),
     vertex.label.color = ifelse(vertexes == sourceGene, "black", "grey30"),
     vertex.label.family = "sans",
     vertex.label.cex = ifelse(vertexes == sourceGene, 0.9, 0.6),
     edge.curved = curve_multiple(net, start = 0.06)) 
# fine tune "start" to get better visual effect. To shrink multiple edges using smaller value.
legend(-1, -1, legend= names(pal_loc),
       col= pal_loc, lty=1, cex=0.6,
       box.lty=1, box.lwd=1, box.col="black")
dev.off()