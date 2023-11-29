setwd("~/MicroDynamics/02_Data")
# Data loading
library(readxl)
otu_bacter<-read.csv("3_Fungi.csv")
meta_data<-read.csv("clin_attr_long.csv")

# Severe patient's data
severe<-meta_data$sample.number[which(meta_data$disease.severity == "severe")]

severe_data<-list()
for(i in 1:length(severe)){
  severe_data[[i]]<-subset(otu_bacter, grepl(severe[i], PatientID))
}
severe_data<-do.call(rbind, severe_data)
severe_data<-severe_data[,-1]

# Moderate patient's data
moderate<-meta_data$sample.number[which(meta_data$disease.severity == "moderate")]

moderate_data<-list()
for(i in 1:length(moderate)){
  moderate_data[[i]]<-subset(otu_bacter, grepl(moderate[i], PatientID))
}
moderate_data<-do.call(rbind, moderate_data)
moderate_data<-moderate_data[,-1]

# Mild patient's data
mild<-meta_data$sample.number[which(meta_data$disease.severity == "mild")]

mild_data<-list()
for(i in 1:length(mild)){
  mild_data[[i]]<-subset(otu_bacter, grepl(mild[i], PatientID))
}
mild_data<-do.call(rbind, mild_data)
mild_data<-mild_data[,-1]

# Networks' inference

library(minet)
library(igraph)

# Severe's network
mim <- build.mim(severe_data, estimator="spearman")
aracne_mat <- aracne(mim)
severeNet<-graph.adjacency(aracne_mat)
severeNet<-as.undirected(severeNet)

# Moderate's network
mim <- build.mim(moderate_data, estimator="spearman")
aracne_mat <- aracne(mim)
moderateNet<-graph.adjacency(aracne_mat)
moderateNet<-as.undirected(moderateNet)

# Mild's network
mim <- build.mim(mild_data, estimator="spearman")
aracne_mat <- aracne(mim)
mildNet<-graph.adjacency(aracne_mat)
mildNet<-as.undirected(mildNet)

# muxViz object
library(muxViz)
library(mlBioNets)
fungi_ml<-list(severeNet, moderateNet, mildNet)
fungi_ml<-ctr_ml(fungi_ml, "degree")
colsmat<-node_color_mat(fungi_ml, "centrality")
abs_data<-list(severe_data, moderate_data, mild_data)
F_abs_mat<-abs_mat(abs_data, fungi_ml, 20)

# Degree 3D plot
lay <- layoutMultiplex(fungi_ml, layout="fr", ggplot.format=F, box=T)
plot_multiplex3D(fungi_ml, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=c("Severe", "Moderate", "Mild"), layer.labels.cex=1.5,
                 node.size.values="auto", node.size.scale=F_abs_mat,
                 node.colors=colsmat,
                 edge.colors="black",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Networks' connectivity analysis

degree_df <- data.frame(Species = vertex.attributes(fungi_ml[[1]])$name,
                        Severe = degree(fungi_ml[[1]]),
                        Moderate = degree(fungi_ml[[2]]),
                        Mild = degree(fungi_ml[[3]]))

degree_df<-degree_df[-which(rowSums(degree_df[,c(2,3,4)])/
                              mean(rowSums(degree_df[,c(2,3,4)]))<0.9),]

# Plot of degree between layers

library(gridExtra)
library(ggplot2)
p1<-ggplot(data=degree_df, aes(x=Species, y=Severe, fill=NULL)) +
  geom_bar(stat="identity")
p2<-ggplot(data=degree_df, aes(x=Species, y=Moderate, fill=NULL)) +
  geom_bar(stat="identity")
p3<-ggplot(data=degree_df, aes(x=Species, y=Mild, fill=NULL)) +
  geom_bar(stat="identity")

grid.arrange(p1 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.position = "none") + ylab("Severe"),
             p2 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.position = "none") + ylab("Moderate"),
             p3 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.position = "none") + ylab("Mild"),
             ncol=3)

# Diversity analysis

# Shannon diversity
Div_Shannon <- function(abundancias_ab){
  abs_rel <- abundancias_ab/sum(abundancias_ab)
  Shannon <- -sum(abs_rel*log(abs_rel))
  return(Shannon)
}

# Simpson dominance
Dom_Simpson <- function(abundancias_ab){
  abs_rel <- abundancias_ab/sum(abundancias_ab)
  Simpson <- sum(abs_rel^2)
  return(Simpson)
}

# Pielou evenness
Eq_Pielou <- function(abundancias_ab){
  abs_rel <- abundancias_ab/sum(abundancias_ab)
  Shannon <- -sum(abs_rel*log(abs_rel))
  Pielou <- Shannon/log(length(abundancias_ab))
  return(Pielou)
}

# Abundance-table diversity function
ab_table_div<-function(ab_table, diversity_type){
  require(vegan)
  if(diversity_type == "shannon"){
    div_table<-c()
    for(i in 1:ncol(ab_table)){
      div_table[i]<-diversity(ab_table[, i])
    }
  } else
    if(diversity_type == "simpson"){
      div_table<-c()
      for(i in 1:ncol(ab_table)){
        div_table[i]<-Dom_Simpson(ab_table[, i])
      }
    } else
      if(diversity_type == "pielou"){
        div_table<-c()
        for(i in 1:ncol(ab_table)){
          S <- specnumber(ab_table[, i])
          div_table[i] <- diversity(ab_table[, i])/log(S)
        }
      } else
        if(diversity_type == "ginisimpson"){
          div_table<-c()
          for(i in 1:ncol(ab_table)){
            div_table[i]<-1-Dom_Simpson(ab_table[, i])
          }
        }
  return(div_table)
}

# Dominance between patients
dom_severe<-ab_table_div(t(severe_data), "simpson")
dom_severe<-data.frame("Severity" = rep("Severe", length(dom_severe)),
                       "Data" = dom_severe)
dom_moderate<-ab_table_div(t(moderate_data), "simpson")
dom_moderate<-data.frame("Severity" = rep("Moderate", length(dom_moderate)),
                         "Data" = dom_moderate)
dom_mild<-ab_table_div(t(mild_data), "simpson")
dom_mild<-data.frame("Severity" = rep("Mild", length(dom_mild)),
                     "Data" = dom_mild)
dom_severity<-rbind(dom_severe, dom_moderate, dom_mild)

# Violin plot
library(ggplot2)
ggplot(dom_severity, aes(x=dom_severity$Severity, y=dom_severity$Data)) +
  geom_violin(fill="darkslategray3", color="black") +
  geom_boxplot(width=0.15, notch=TRUE,
               fill=c("green3", "darkorange1", "red3"),
               color="black") + ggtitle("Simpson dominance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + ylab("")


# Rank-abundance analysis

library(RADanalysis)
sample_classes <- c(rep(1, nrow(dom_severe)),rep(2, nrow(dom_moderate)),
                    rep(3, nrow(dom_mild)))
line_cols <- c("red3","darkorange1","green3")

severity_mat<-rbind(severe_data, moderate_data, mild_data)

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
plot(1e10,xlim = c(1,15),ylim = c(0,0.8),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),
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
       legend = c("Severe","Moderate", "Mild"),
       col = line_cols, lwd = 3)


# MSD analysis

d <- dist(x = n_severity_mat,method = "manhattan")
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

lfc1<-log_fc(fungi_ml, c("Severe", "Moderate", "Mild"), "Mild", "Moderate")
lfc2<-log_fc(fungi_ml, c("Severe", "Moderate", "Mild"), "Mild", "Severe")

# Plot log-fold change

library(ggpubr)

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

# Dominance violin plot
ggbetweenstats(
  data  = dom_severity,
  x     = Severity,
  y     = Data
)
