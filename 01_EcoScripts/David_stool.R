library(seqtime)
library(igraph)
library(phyloseq)
library(mlBioNets)

# David stool A - Travel

data("david_stool_lineages")
data("david_stoolA_otus")
# Sort abundances
sorted <- sort(apply(david_stoolA_otus,1,sum),
               decreasing = TRUE, index.return = TRUE)
# Plot the most abundant taxa
tsplot(david_stoolA_otus[sorted$ix[1:10],])
# Order the data
O_table <- david_stoolA_otus
T_table <- david_stool_lineages
# Table collapse
T_Collapsed<-T_collapse(F, T_table = T_table, O_table = O_table,
                        names_level = "V7")
# Clean data
T_Collapsed<-T_Collapsed[, -c(which(colnames(T_Collapsed) == "none"))]
dim(T_Collapsed)
# Separate data
T1_mat <- T_Collapsed[1:70,]
T2_mat <- T_Collapsed[71:122,]
T3_mat <- T_Collapsed[123:dim(T_Collapsed)[1],]

# Network's inference
library(minet)

mim <- build.mim(T1_mat,estimator="spearman")
T_matA <- aracne(mim)
T1_aracne<-graph.adjacency(T_matA)
T1_aracne<-as.undirected(T1_aracne)
plot_network(T1_aracne)

mim <- build.mim(T2_mat,estimator="spearman")
T_matB <- aracne(mim)
T2_aracne<-graph.adjacency(T_matB)
T2_aracne<-as.undirected(T2_aracne)
plot_network(T2_aracne)

mim <- build.mim(T3_mat,estimator="spearman")
T_matC <- aracne(mim)
T3_aracne<-graph.adjacency(T_matC)
T3_aracne<-as.undirected(T3_aracne)
plot_network(T3_aracne)

# Degree distribution
hist(degree(T1_aracne), freq = FALSE, main = "Before travel", ylab = "Density",
     col = "khaki")
lines(dx1, lwd = 2, col = "red")
rug(jitter(degree(T1_aracne)))

hist(degree(T2_aracne), freq = FALSE, main = "During travel", ylab = "Density",
     col = "khaki")
lines(dx2, lwd = 2, col = "blue")
rug(jitter(degree(T2_aracne)))

hist(degree(T3_aracne), freq = FALSE, main = "After travel", ylab = "Density",
     col = "khaki")
lines(dx3, lwd = 2, col = "green")
rug(jitter(degree(T2_aracne)))

par(mfrow=c(1,1))

dx1 <- density(degree(T1_aracne), bw = 0.1)
dx2 <- density(degree(T2_aracne), bw = 0.1)
dx3 <- density(degree(T3_aracne), bw = 0.1)

plot(dx1, lwd = 2, main = "Degree dendities", xlab = "",
     col = "red", xlim = c(-4, 6), ylim = c(0, 0.5))
rug(jitter(degree(T1_aracne)), col = "red")

lines(dx2, lwd = 2, col = "blue")
rug(jitter(degree(T2_aracne)), col = "blue")

lines(dx3, lwd = 2, col = "green")
rug(jitter(degree(T3_aracne)), col = "green")

Grupos<-c("Before travel", "During travel", "After travel")
legend("topleft", Grupos, col = c("Red", "Blue", "Green"), lty = 1) 

unq<-unique(T_table[,"V3"])
colors<-rainbow(length(unq))
T1_aracne<-v_colored(T1_aracne, T_table, g_tax = "V3",
                     p_tax = "V7", g_colors = colors)
T2_aracne<-v_colored(T2_aracne, T_table, g_tax = "V3",
                     p_tax = "V7", g_colors = colors)
T3_aracne<-v_colored(T3_aracne, T_table, g_tax = "V3",
                     p_tax = "V7", g_colors = colors)

T1_aracne<-g_abundance(layer_mat = T1_mat, g = T1_aracne)
T2_aracne<-g_abundance(layer_mat = T2_mat, g = T2_aracne)
T3_aracne<-g_abundance(layer_mat = T3_mat, g = T3_aracne)

library(muxViz)
g.list<-list(T1_aracne, T2_aracne, T3_aracne)
g.list<-ctr_ml(g.list, "degree")
matctr<-node_color_mat(g.list, "centrality")
matsize<-abs_mat(list(T1_mat, T2_mat, T3_mat), g.list, 10)

lay <- layoutMultiplex(g.list, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(g.list, layer.layout=lay,
                 layer.colors=c("green3", "red3", "blue3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL,
                 layer.labels.cex=1.5, node.size.values="auto",
                 node.size.scale=matsize,
                 node.colors=matctr,
                 edge.colors="black",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Networks' connectivity analysis

degree_df <- data.frame(Species = vertex.attributes(g.list[[1]])$name,
                        Before = degree(g.list[[1]]),
                        During = degree(g.list[[2]]),
                        After = degree(g.list[[3]]))

degree_df<-degree_df[-which(rowSums(degree_df[,c(2,3,4)])/
                              mean(rowSums(degree_df[,c(2,3,4)]))<0.9),]

# Plot of degree between layers

library(gridExtra)
library(ggplot2)
p1<-ggplot(data=degree_df, aes(x=Species, y=Before, fill=NULL)) +
  geom_bar(stat="identity")
p2<-ggplot(data=degree_df, aes(x=Species, y=During, fill=NULL)) +
  geom_bar(stat="identity")
p3<-ggplot(data=degree_df, aes(x=Species, y=After, fill=NULL)) +
  geom_bar(stat="identity")

grid.arrange(p1 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.position = "none") + ylab("Before travel"),
             p2 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.position = "none") + ylab("During travel"),
             p3 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.position = "none") + ylab("After travel"),
             ncol=3)

# Rank-abundance analysis

library(RADanalysis)
sample_classes <- c(rep(1, nrow(T1_mat)),rep(2, nrow(T2_mat)),
                    rep(3, nrow(T3_mat)))
line_cols <- c("green3", "red3", "blue3")

stool_mat<-rbind(T1_mat, T2_mat, T3_mat)

# Normalized abundances
n_stool_mat<-stool_mat/rowSums(stool_mat)
n_stool_mat<-as.matrix(n_stool_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_stool_mat)){
  sorted_abs[[i]]<-sort(n_stool_mat[i,], decreasing = T)
}
n_stool_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                       ncol=length(sorted_abs[[1]]), byrow=TRUE)

# Plot the axis
plot(1e10,xlim = c(1,15),ylim = c(0,0.4),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4),
               las = 0)

# Rank-abundance colors per disease severity
a <- representative_RAD(norm_rad = n_stool_mat, sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_stool_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_stool_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Before travel","During travel", "After travel"),
       col = line_cols, lwd = 3)

# MSD analysis

d <- dist(x = n_stool_mat,method = "manhattan")
mds <- cmdscale(d = d,k = 5,eig = TRUE)
plot(mds$points,xlab = "First coordinate",ylab = "Second coordinate",pch = 19,
     cex =1,col = line_cols[sample_classes],
     main = "MDS plot with representative points \n of each group and error bars")
a <- representative_point(input = mds$points,ids = which(sample_classes == 1),
                          col = scales::alpha(line_cols[1],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 2),
                          col = scales::alpha(line_cols[2],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 3),
                          col = scales::alpha(line_cols[3],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
legend("topright",bty = "n",legend = c("Severe","Moderate","Mild"),
       col = line_cols,pch = 19)

# Phyloseq data
T_otus <- otu_table(t(T_Collapsed), taxa_are_rows = T)

sampledata <- sample_data(data.frame(
  Samples = c(rep("basal", nrow(T1_mat)), rep("pert", nrow(T2_mat)),
              rep("rec", nrow(T3_mat))),
  Ab_mean = colSums(otu_table(T_otus))/length(otu_table(T_otus)),
  row.names = colnames(otu_table(T_otus)),
  stringsAsFactors = F
))

davidA_physeq <- merge_phyloseq(T_otus, sample_data(sampledata))

# Dysbiosis analysis

library(dysbiosisR)
# Bray-Curtis distance matrix
dist.mat <- phyloseq::distance(davidA_physeq, "bray")
# Get reference samples
ref.samples <- sample_names(subset_samples(davidA_physeq, 
                                           Samples == "basal"))
# Community level variation analysis
dysbiosis_1 <- dysbiosisMedianCLV(davidA_physeq,
                                  dist_mat = dist.mat,
                                  reference_samples = ref.samples)

# We sample the data set identifying as dysbiotic the data under the 90th percentile
dysbiosis_thres <- quantile(subset(dysbiosis_1, Samples == "pert")$score, 0.9)
normobiosis_thres <- quantile(subset(dysbiosis_1, Samples == "pert")$score, 0.1)

library(dplyr)
dysbiosis_1 <- dysbiosis_1 |> 
  mutate(isDysbiostic = ifelse(score >= dysbiosis_thres, TRUE, FALSE))

# Dysbiosis plot measures according to CLV method
p1 <- plotDysbiosis(df=dysbiosis_1,
                    xvar="Samples",
                    yvar="score",
                    colors=c(basal="green4", pert="red3",
                             post="orange3", rec="blue4"),
                    show_points = F) +
  labs(x="", y="Dysbiosis Score") +
  theme_bw(base_size = 14)
p1

# Dysbiosis plot measures according to euclidean method
dysbiosis_2 <- euclideanDistCentroids(davidA_physeq,
                                      dist_mat = dist.mat,
                                      use_squared = TRUE,
                                      group_col = "Samples",
                                      control_label = "basal",
                                      case_label = "pert")

p2 <- plotDysbiosis(df=dysbiosis_2,
                    xvar="Samples",
                    yvar="CentroidDist_score",
                    colors=c(basal="green4", pert="red3", rec="blue4"),
                    show_points = FALSE) +
  labs(x="", y="Dysbiosis Score (Centroid)") +
  theme_bw(base_size = 14)
p2


# Log-fold change analysis

log_fc<-function(g.list, layer_names, control_layer, test_layer){
  centralities<-list()
  for(i in 1:length(g.list)){
    centralities[[i]]<-degree(g.list[[i]])
  }
  d_mat<-matrix(unlist(centralities), length(centralities[[1]]),
                length(centralities))
  df_degree<-as.data.frame(d_mat)
  rownames(df_degree)<-vertex.attributes(g.list[[1]])$name
  colnames(df_degree)<-layer_names
  control_n<-which(colnames(df_degree)==control_layer)
  test_n<-which(colnames(df_degree)==test_layer)
  log_fc<--log((df_degree[,control_n]+1)/(df_degree[,test_n]+1), 2)
  df_logFC<-data.frame(
    node_names = vertex.attributes(g.list[[1]])$name,
    log_fc = log_fc
  )
  not_fc <- which(df_logFC$log_fc == 0)
  df_logFC<-df_logFC[-not_fc,]
  return(df_logFC)
}

lfc1<-log_fc(g.list, c("Basal", "Perturbated", "Recovered"), "Basal", "Perturbated")
lfc2<-log_fc(g.list, c("Basal", "Perturbated", "Recovered"), "Basal", "Recovered")

# Plot log-fold change

library(ggpubr)
library(stringr)

ggbarplot(lfc1, x = "node_names", y = "log_fc",
          fill = "node_names",
          color = "blue",
          sort.val = "desc",
          sort.by.groups = FALSE,
          x.text.angle = 90,
          ylab = "Log-fold change",
          rotate = TRUE,
          ggtheme = theme_minimal()) +
  scale_x_discrete(labels = function(x) str_wrap(x, 1)) +
  guides(fill = F)

ggbarplot(lfc2, x = "node_names", y = "log_fc",
          fill = "node_names",
          color = "blue",
          sort.val = "desc",
          sort.by.groups = FALSE,
          x.text.angle = 90,
          ylab = "Log-fold change",
          rotate = TRUE,
          ggtheme = theme_minimal()) +
  scale_x_discrete(labels = function(x) str_wrap(x, 1)) +
  guides(fill = F)






# David stool B - Salmonella infection

data("david_stoolB_otus")
# Sort abundances
sorted <- sort(apply(david_stoolB_otus,1,sum),
               decreasing = TRUE, index.return = TRUE)
# Plot the most abundant taxa
tsplot(david_stoolB_otus[sorted$ix[1:10],])
# Order the data
O_table <- david_stoolB_otus
# Table collapse
T_Collapsed<-T_collapse(T_table = T_table, O_table = O_table,
                        names_level = "V7", is_phyloseq = F)
# Clean data
T_Collapsed<-T_Collapsed[, -c(which(colnames(T_Collapsed) == "none"))]
dim(T_Collapsed)
# Separate data
T1_mat <- T_Collapsed[1:150,]
T2_mat <- T_Collapsed[151:159,]
T3_mat <- T_Collapsed[160:dim(T_Collapsed)[1],]

# Network's inference
library(minet)

mim <- build.mim(T1_mat,estimator="spearman")
T_matA <- aracne(mim)
T1_aracne<-graph.adjacency(T_matA)
T1_aracne<-as.undirected(T1_aracne)
plot_network(T1_aracne)

mim <- build.mim(T2_mat,estimator="spearman")
T_matB <- aracne(mim)
T2_aracne<-graph.adjacency(T_matB)
T2_aracne<-as.undirected(T2_aracne)
plot_network(T2_aracne)

mim <- build.mim(T3_mat,estimator="spearman")
T_matC <- aracne(mim)
T3_aracne<-graph.adjacency(T_matC)
T3_aracne<-as.undirected(T3_aracne)
plot_network(T3_aracne)

# Degree distribution
hist(degree(T1_aracne), freq = FALSE, main = "Before infection", ylab = "Density",
     col = "khaki")
lines(dx1, lwd = 2, col = "red")
rug(jitter(degree(T1_aracne)))

hist(degree(T2_aracne), freq = FALSE, main = "During infection", ylab = "Density",
     col = "khaki")
lines(dx2, lwd = 2, col = "blue")
rug(jitter(degree(T2_aracne)))

hist(degree(T3_aracne), freq = FALSE, main = "After infection", ylab = "Density",
     col = "khaki")
lines(dx3, lwd = 2, col = "green")
rug(jitter(degree(T2_aracne)))

par(mfrow=c(1,1))

dx1 <- density(degree(T1_aracne), bw = 0.1)
dx2 <- density(degree(T2_aracne), bw = 0.1)
dx3 <- density(degree(T3_aracne), bw = 0.1)

plot(dx1, lwd = 2, main = "Degree dendities", xlab = "",
     col = "red", xlim = c(-4, 6), ylim = c(0, 0.5))
rug(jitter(degree(T1_aracne)), col = "red")

lines(dx2, lwd = 2, col = "blue")
rug(jitter(degree(T2_aracne)), col = "blue")

lines(dx3, lwd = 2, col = "green")
rug(jitter(degree(T3_aracne)), col = "green")

Grupos<-c("Before infection", "During infection", "After infection")
legend("topleft", Grupos, col = c("Red", "Blue", "Green"), lty = 1) 

unq<-unique(T_table[,"V3"])
colors<-rainbow(length(unq))
T1_aracne<-v_colored(T1_aracne, T_table, g_tax = "V3",
                     p_tax = "V7", g_colors = colors)
T2_aracne<-v_colored(T2_aracne, T_table, g_tax = "V3",
                     p_tax = "V7", g_colors = colors)
T3_aracne<-v_colored(T3_aracne, T_table, g_tax = "V3",
                     p_tax = "V7", g_colors = colors)

T1_aracne<-g_abundance(layer_mat = T1_mat, g = T1_aracne)
T2_aracne<-g_abundance(layer_mat = T2_mat, g = T2_aracne)
T3_aracne<-g_abundance(layer_mat = T3_mat, g = T3_aracne)

library(muxViz)
g.list<-list(T1_aracne, T2_aracne, T3_aracne)
g.list<-ctr_ml(g.list, "degree")
matctr<-node_color_mat(g.list, "centrality")
matsize<-abs_mat(list(T1_mat, T2_mat, T3_mat), g.list, 20)

lay <- layoutMultiplex(g.list, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(g.list, layer.layout=lay,
                 layer.colors=c("green3", "red3", "blue3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL,
                 layer.labels.cex=1.5, node.size.values="auto",
                 node.size.scale=matsize,
                 node.colors=matctr,
                 edge.colors="black",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Networks' connectivity analysis

degree_df <- data.frame(Species = vertex.attributes(g.list[[1]])$name,
                        Before = degree(g.list[[1]]),
                        During = degree(g.list[[2]]),
                        After = degree(g.list[[3]]))

degree_df<-degree_df[-which(rowSums(degree_df[,c(2,3,4)])/
                              mean(rowSums(degree_df[,c(2,3,4)]))<0.9),]

# Plot of degree between layers

library(gridExtra)
library(ggplot2)
p1<-ggplot(data=degree_df, aes(x=Species, y=Before, fill=NULL)) +
  geom_bar(stat="identity")
p2<-ggplot(data=degree_df, aes(x=Species, y=During, fill=NULL)) +
  geom_bar(stat="identity")
p3<-ggplot(data=degree_df, aes(x=Species, y=After, fill=NULL)) +
  geom_bar(stat="identity")

grid.arrange(p1 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.position = "none") + ylab("Before infection"),
             p2 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.position = "none") + ylab("During infection"),
             p3 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.position = "none") + ylab("After infection"),
             ncol=3)

# Rank-abundance analysis

library(RADanalysis)
sample_classes <- c(rep(1, nrow(T1_mat)),rep(2, nrow(T2_mat)),
                    rep(3, nrow(T3_mat)))
line_cols <- c("green3", "red3", "blue3")

stool_mat<-rbind(T1_mat, T2_mat, T3_mat)

# Normalized abundances
n_stool_mat<-stool_mat/rowSums(stool_mat)
n_stool_mat<-as.matrix(n_stool_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_stool_mat)){
  sorted_abs[[i]]<-sort(n_stool_mat[i,], decreasing = T)
}
n_stool_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                    ncol=length(sorted_abs[[1]]), byrow=TRUE)

# Plot the axis
plot(1e10,xlim = c(1,15),ylim = c(0,0.7),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7),
               las = 0)

# Rank-abundance colors per disease severity
a <- representative_RAD(norm_rad = n_stool_mat, sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_stool_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_stool_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Before travel","During travel", "After travel"),
       col = line_cols, lwd = 3)

# MSD analysis

d <- dist(x = n_stool_mat,method = "manhattan")
mds <- cmdscale(d = d,k = 5,eig = TRUE)
plot(mds$points,xlab = "First coordinate",ylab = "Second coordinate",pch = 19,
     cex =1,col = line_cols[sample_classes],
     main = "MDS plot with representative points \n of each group and error bars")
a <- representative_point(input = mds$points,ids = which(sample_classes == 1),
                          col = scales::alpha(line_cols[1],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 2),
                          col = scales::alpha(line_cols[2],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
a <- representative_point(input = mds$points,ids = which(sample_classes == 3),
                          col = scales::alpha(line_cols[3],0.5),
                          plot = TRUE,standard_error_mean = TRUE,pch = 19, cex = 4)
legend("topright",bty = "n",legend = c("Severe","Moderate","Mild"),
       col = line_cols,pch = 19)



# Phyloseq data
T_otus <- otu_table(t(T_Collapsed), taxa_are_rows = T)

sampledata <- sample_data(data.frame(
  Samples = c(rep("basal", nrow(T1_mat)), rep("pert", nrow(T2_mat)),
              rep("rec", nrow(T3_mat))),
  Ab_mean = colSums(otu_table(T_otus))/length(otu_table(T_otus)),
  row.names = colnames(otu_table(T_otus)),
  stringsAsFactors = F
))

davidA_physeq <- merge_phyloseq(T_otus, sample_data(sampledata))

# Dysbiosis analysis

library(dysbiosisR)
# Bray-Curtis distance matrix
dist.mat <- phyloseq::distance(davidA_physeq, "bray")
# Get reference samples
ref.samples <- sample_names(subset_samples(davidA_physeq, 
                                           Samples == "basal"))
# Community level variation analysis
dysbiosis_1 <- dysbiosisMedianCLV(davidA_physeq,
                                  dist_mat = dist.mat,
                                  reference_samples = ref.samples)

# We sample the data set identifying as dysbiotic the data under the 90th percentile
dysbiosis_thres <- quantile(subset(dysbiosis_1, Samples == "pert")$score, 0.9)
normobiosis_thres <- quantile(subset(dysbiosis_1, Samples == "pert")$score, 0.1)

library(dplyr)
dysbiosis_1 <- dysbiosis_1 |> 
  mutate(isDysbiostic = ifelse(score >= dysbiosis_thres, TRUE, FALSE))

# Dysbiosis plot measures according to CLV method
p1 <- plotDysbiosis(df=dysbiosis_1,
                    xvar="Samples",
                    yvar="score",
                    colors=c(basal="green4", pert="red3",
                             post="orange3", rec="blue4"),
                    show_points = F) +
  labs(x="", y="Dysbiosis Score") +
  theme_bw(base_size = 14)
p1

# Dysbiosis plot measures according to euclidean method
dysbiosis_2 <- euclideanDistCentroids(davidA_physeq,
                                      dist_mat = dist.mat,
                                      use_squared = TRUE,
                                      group_col = "Samples",
                                      control_label = "basal",
                                      case_label = "pert")

p2 <- plotDysbiosis(df=dysbiosis_2,
                    xvar="Samples",
                    yvar="CentroidDist_score",
                    colors=c(basal="green4", pert="red3", rec="blue4"),
                    show_points = FALSE) +
  labs(x="", y="Dysbiosis Score (Centroid)") +
  theme_bw(base_size = 14)
p2


# Log-fold change analysis

lfc1<-log_fc(g.list, c("Basal", "Perturbated", "Recovered"), "Basal", "Perturbated")
lfc2<-log_fc(g.list, c("Basal", "Perturbated", "Recovered"), "Basal", "Recovered")

# Plot log-fold change

library(ggpubr)
library(stringr)

ggbarplot(lfc1, x = "node_names", y = "log_fc",
          fill = "node_names",
          color = "blue",
          sort.val = "desc",
          sort.by.groups = FALSE,
          x.text.angle = 90,
          ylab = "Log-fold change",
          rotate = TRUE,
          ggtheme = theme_minimal()) +
  scale_x_discrete(labels = function(x) str_wrap(x, 1)) +
  guides(fill = F)

ggbarplot(lfc2, x = "node_names", y = "log_fc",
          fill = "node_names",
          color = "blue",
          sort.val = "desc",
          sort.by.groups = FALSE,
          x.text.angle = 90,
          ylab = "Log-fold change",
          rotate = TRUE,
          ggtheme = theme_minimal()) +
  scale_x_discrete(labels = function(x) str_wrap(x, 1)) +
  guides(fill = F)
