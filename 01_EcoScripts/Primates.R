# Data loading

# Abundance data
setwd("~/MicroDynamics/02_Data")
abs_bf <- read.table(file = "16S_ASV_filtered_Rel_Abun.txt", header = TRUE)
rnames <- abs_bf[,1]
abs_bf <- abs_bf[,-1]
rownames(abs_bf)<-rnames

# Taxonomic data
library(readr)
taxdata <- read_tsv(file = "16S_taxonomy.tsv", col_names = T)
taxdata<-as.data.frame(taxdata)
library(tidyverse)
taxdata<-taxdata %>%
  separate(Taxon, into = paste0("TL", 1:7), sep = ";")
r2names<-taxdata[,1]
taxdata<-taxdata[,-1]
rownames(taxdata)<-r2names

# Collapse the table

# Unique taxa identification
unq<-unique(taxdata[,"TL6"]); unq
# Sum of abundances of the same taxa
abs_genus <- matrix(, nrow = length(abs_bf), ncol = 0)
for(i in 1:length(unq)){
  a<-which(taxdata[,"TL6"] %in% unq[i]); a
  b<-as.matrix(colSums(abs_bf[i,]))
  abs_genus<-cbind(abs_genus, b)
}
# Taxa names assignment
colnames(abs_genus)<-unq

# Clean the data bases
abs_genus<-abs_genus[,-which(colnames(abs_genus) == " g__")]
abs_genus<-abs_genus[,-which(is.na(colnames(abs_genus)))]
rownames(abs_genus) <- gsub("^X", "", rownames(abs_genus))

# Load meta data
metadata<-read_table("Metadata_Filtered.txt")
metadata$Sample_Type

# Separate data
types <- unique(metadata$Sample_Type)

# Ba Aka humans
BaAka <- metadata$SampleID[which(metadata$Sample_Type == "BaAka-Human")]
T_BaAka <- abs_genus[which(rownames(abs_genus) %in% BaAka),]
# Bantu humans
Bantu <- metadata$SampleID[which(metadata$Sample_Type == "Bantu-Human")]
T_Bantu <- abs_genus[which(rownames(abs_genus) %in% Bantu),]
# USA humans
USA <- metadata$SampleID[which(metadata$Sample_Type == "USA-Human")]
T_USA <- abs_genus[which(rownames(abs_genus) %in% USA),]

# Network inference

library(minet)
library(igraph)

# Ba Aka humans
mim <- build.mim(T_BaAka, estimator="spearman")
aracne_mat <- aracne(mim)
BaAka_Net<-graph.adjacency(aracne_mat)
BaAka_Net<-as.undirected(BaAka_Net)
# Bantu humans
mim <- build.mim(T_Bantu, estimator="spearman")
aracne_mat <- aracne(mim)
Bantu_Net<-graph.adjacency(aracne_mat)
Bantu_Net<-as.undirected(Bantu_Net)
# USA humans
mim <- build.mim(T_USA, estimator="spearman")
aracne_mat <- aracne(mim)
USA_Net<-graph.adjacency(aracne_mat)
USA_Net<-as.undirected(USA_Net)

# Color vector
phylunq<-unique(taxdata$TL2)
phylunq<-phylunq[-c(which(phylunq == " p__"), which(is.na(phylunq)))]
colors<-rainbow(length(phylunq))

# Multilayer plot
humans_ml<-list(BaAka_Net, Bantu_Net, USA_Net)
humans_ml<-v_colored_ml(humans_ml, taxdata, "TL2", "TL6", colors)
humans_abs<-list(T_BaAka, T_Bantu, T_USA)
humans_abs<-abs_mat(humans_abs, humans_ml, 20)
phylacols<-node_color_mat(humans_ml, "phylo")

library(muxViz)
lay <- layoutMultiplex(humans_ml, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(humans_ml, layer.layout=lay,
                 layer.colors=c("green3", "skyblue", "red3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto", node.size.scale=humans_abs,
                 node.colors=phylacols,
                 edge.colors="black",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Rad-abundance analysis
library(RADanalysis)
sample_classes <- c(rep(1, nrow(T_BaAka)),rep(2, nrow(T_Bantu)),
                    rep(3, nrow(T_USA)))
line_cols <- c("green3", "skyblue", "red3")

severity_mat<-rbind(T_BaAka, T_Bantu, T_USA)
# Normalized abundances
n_severity_mat<-severity_mat/rowSums(severity_mat)
n_severity_mat<-as.matrix(n_severity_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_severity_mat)){
  sorted_abs[[i]]<-sort(n_severity_mat[i,], decreasing = T)
}
n_severity_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                       ncol=length(sorted_abs[[1]]), byrow=TRUE)

# Plot the axis
plot(1e10,xlim = c(1,15),ylim = c(0,0.3),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3),
               las = 0)

# Rank-abundance colors per disease severity
a <- representative_RAD(norm_rad = n_severity_mat, sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_severity_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_severity_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Ba Aka","Bantu", "USA"),
       col = line_cols, lwd = 3)

# MSD analysis

d <- dist(x = n_severity_mat,method = "manhattan")
mds <- cmdscale(d = d,k = 5,eig = TRUE)
plot(mds$points,xlab = "First coordinate",ylab = "Second coordinate",pch = 19,
     cex =1,col = line_cols[sample_classes],
     main = "MDS plot with representative points \n of ba aka, bantu and americans")
a <- representative_point(input = mds$points,ids = which(sample_classes == 1),
                          col = scales::alpha(line_cols[1],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 2),
                          col = scales::alpha(line_cols[2],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 3),
                          col = scales::alpha(line_cols[3],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
legend("topleft",bty = "n",legend = c("Ba Aka","Bantu", "USA"),
       col = line_cols,pch = 19)


# Separate data

# Savage chimpanzees
Chimps <- metadata$SampleID[which(metadata$Sample_Type == "Chimps")]
T_Chimps <- abs_genus[which(rownames(abs_genus) %in% Chimps),]

mim <- build.mim(T_Chimps, estimator="spearman")
aracne_mat <- aracne(mim)
Chimps_Net<-graph.adjacency(aracne_mat)
Chimps_Net<-as.undirected(Chimps_Net)

# Captive chimpanzees
Captive <- metadata$SampleID[which(metadata$Sample_Type == "Captive")]
T_Captive <- abs_genus[which(rownames(abs_genus) %in% Captive),]

mim <- build.mim(T_Captive, estimator="spearman")
aracne_mat <- aracne(mim)
Captive_Net<-graph.adjacency(aracne_mat)
Captive_Net<-as.undirected(Captive_Net)

# Multilayer plot
chimp_ml<-list(Chimps_Net, Captive_Net)
chimp_ml<-v_colored_ml(chimp_ml, taxdata, "TL2", "TL6", colors)
chimp_abs<-list(T_Chimps, T_Captive)
chimp_abs<-abs_mat(chimp_abs, chimp_ml, 2)
phylacols<-node_color_mat(chimp_ml, "phylo")

lay <- layoutMultiplex(chimp_ml, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(chimp_ml, layer.layout=lay,
                 layer.colors=c("green3", "red3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto", node.size.scale=chimp_abs,
                 node.colors=phylacols,
                 edge.colors="black",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Rad-abundance analysis
library(RADanalysis)
sample_classes <- c(rep(1, nrow(T_Chimps)),rep(2, nrow(T_Captive)))
line_cols <- c("green3", "red3")

severity_mat<-rbind(T_Chimps, T_Captive)
# Normalized abundances
n_severity_mat<-severity_mat/rowSums(severity_mat)
n_severity_mat<-as.matrix(n_severity_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_severity_mat)){
  sorted_abs[[i]]<-sort(n_severity_mat[i,], decreasing = T)
}
n_severity_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                       ncol=length(sorted_abs[[1]]), byrow=TRUE)

# Plot the axis
plot(1e10,xlim = c(1,15),ylim = c(0,0.2),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2),
               las = 0)

# Rank-abundance colors per disease severity
a <- representative_RAD(norm_rad = n_severity_mat, sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_severity_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Wild","Captivity"),
       col = line_cols, lwd = 3)

# MSD analysis

d <- dist(x = n_severity_mat,method = "manhattan")
mds <- cmdscale(d = d,k = 5,eig = TRUE)
plot(mds$points,xlab = "First coordinate",ylab = "Second coordinate",pch = 19,
     cex =1,col = line_cols[sample_classes],
     main = "MDS plot with representative points of wild \n chimpanzees vs chimpanzees in captivity")
a <- representative_point(input = mds$points,ids = which(sample_classes == 1),
                          col = scales::alpha(line_cols[1],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 2),
                          col = scales::alpha(line_cols[2],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
legend("topleft",bty = "n",legend = c("Wild","Captivity"),
       col = line_cols,pch = 19)

