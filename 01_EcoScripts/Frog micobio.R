# Data of frog's micobiota under different temperatures
setwd("C:/Users/LaV_V/Downloads")
# Data loading
library(readxl)
meta_hongos<-read.csv("Metadata_Hongos.csv")
hongos<-read_excel("Hongos.xlsx")

# Identify frog´s life stages
life_stage<-unique(meta_hongos$Life_stage)

# Separation of data by life stages type
Tadpole<-meta_hongos[which(meta_hongos$Life_stage == life_stage[1]),]
Metamorphic<-meta_hongos[which(meta_hongos$Life_stage == life_stage[2]),]
Sub_adult<-meta_hongos[which(meta_hongos$Life_stage == life_stage[3]),]

# Tadpole

# Treatment identification
Tad_T1<-Tadpole$SampleID[which(Tadpole$Treatment=="Treatment1")]
Tad_T2<-Tadpole$SampleID[which(Tadpole$Treatment=="Treatment2")]
Tad_Ctr<-Tadpole$SampleID[which(Tadpole$Treatment=="Control")]
# Separation of data by treatment types
HTad_T1<-as.data.frame(hongos[,which(colnames(hongos) %in% Tad_T1)])
HTad_T2<-as.data.frame(hongos[,which(colnames(hongos) %in% Tad_T2)])
HTad_Ctr<-as.data.frame(hongos[,which(colnames(hongos) %in% Tad_Ctr)])

# Metamorphic

# Treatment identification
Met_T1<-Metamorphic$SampleID[which(Metamorphic$Treatment=="Treatment1")]
Met_T2<-Metamorphic$SampleID[which(Metamorphic$Treatment=="Treatment2")]
Met_Ctr<-Metamorphic$SampleID[which(Metamorphic$Treatment=="Control")]
# Separation of data by treatment types
HMet_T1<-as.data.frame(hongos[,which(colnames(hongos) %in% Met_T1)])
HMet_T2<-as.data.frame(hongos[,which(colnames(hongos) %in% Met_T2)])
HMet_Ctr<-as.data.frame(hongos[,which(colnames(hongos) %in% Met_Ctr)])

# Sub-adult

# Treatment identification
Adl_T1<-Sub_adult$SampleID[which(Sub_adult$Treatment=="Treatment1")]
Adl_T2<-Sub_adult$SampleID[which(Sub_adult$Treatment=="Treatment2")]
Adl_Ctr<-Sub_adult$SampleID[which(Sub_adult$Treatment=="Control")]
# Separation of data by treatment types
HAdl_T1<-as.data.frame(hongos[,which(colnames(hongos) %in% Adl_T1)])
HAdl_T2<-as.data.frame(hongos[,which(colnames(hongos) %in% Adl_T2)])
HAdl_Ctr<-as.data.frame(hongos[,which(colnames(hongos) %in% Adl_Ctr)])

# Taxa table
t1<-which(colnames(hongos)=="Kingdom")
t2<-which(colnames(hongos)=="Species")
tax_fungi<-as.data.frame(hongos[,t1:t2])

# Collapse function
T_collapse<-function(T_table, O_table, names_level){
  # Unique taxa identification
  unq<-unique(T_table[,names_level]); unq
  # Sum of abundances of the same taxa
  mat <- matrix(, nrow = length(O_table), ncol = 0)
  for(i in 1:length(unq)){
    a<-which(T_table[,names_level] %in% unq[i]); a
    b<-as.matrix(colSums(O_table[a,]))
    mat<-cbind(mat, b)
  }
  # Taxa names assignment
  colnames(mat)<-unq
  return(mat)
}

# Abundance tables' collapse

# Tadpole under treatment 1
FTad_T1<-T_collapse(tax_fungi, HTad_T1, "Genus")
FTad_T1<-FTad_T1[, -which(colnames(FTad_T1) == "unidentified")]
FTad_T1<-FTad_T1[, -which(is.na(colnames(FTad_T1)))]

# Tadpole under treatment 2
FTad_T2<-T_collapse(tax_fungi, HTad_T2, "Genus")
FTad_T2<-FTad_T2[, -which(colnames(FTad_T2) == "unidentified")]
FTad_T2<-FTad_T2[, -which(is.na(colnames(FTad_T2)))]

# Tadpole control
FTad_Ctr<-T_collapse(tax_fungi, HTad_Ctr, "Genus")
FTad_Ctr<-FTad_Ctr[, -which(colnames(FTad_Ctr) == "unidentified")]
FTad_Ctr<-FTad_Ctr[, -which(is.na(colnames(FTad_Ctr)))]

# Metamorphic under treatment 1
FMet_T1<-T_collapse(tax_fungi, HMet_T1, "Genus")
FMet_T1<-FMet_T1[, -which(colnames(FMet_T1) == "unidentified")]
FMet_T1<-FMet_T1[, -which(is.na(colnames(FMet_T1)))]

# Metamorphic under treatment 2
FMet_T2<-T_collapse(tax_fungi, HMet_T2, "Genus")
FMet_T2<-FMet_T2[, -which(colnames(FMet_T2) == "unidentified")]
FMet_T2<-FMet_T2[, -which(is.na(colnames(FMet_T2)))]

# Metamorphic control
FMet_Ctr<-T_collapse(tax_fungi, HMet_Ctr, "Genus")
FMet_Ctr<-FMet_Ctr[, -which(colnames(FMet_Ctr) == "unidentified")]
FMet_Ctr<-FMet_Ctr[, -which(is.na(colnames(FMet_Ctr)))]

# Sub-adult under treatment 1
FAdl_T1<-T_collapse(tax_fungi, HAdl_T1, "Genus")
FAdl_T1<-FAdl_T1[, -which(colnames(FAdl_T1) == "unidentified")]
FAdl_T1<-FAdl_T1[, -which(is.na(colnames(FAdl_T1)))]

# Sub-adult under treatment 2
FAdl_T2<-T_collapse(tax_fungi, HAdl_T2, "Genus")
FAdl_T2<-FAdl_T2[, -which(colnames(FAdl_T2) == "unidentified")]
FAdl_T2<-FAdl_T2[, -which(is.na(colnames(FAdl_T2)))]

# Sub-adult control
FAdl_Ctr<-T_collapse(tax_fungi, HAdl_Ctr, "Genus")
FAdl_Ctr<-FAdl_Ctr[, -which(colnames(FAdl_Ctr) == "unidentified")]
FAdl_Ctr<-FAdl_Ctr[, -which(is.na(colnames(FAdl_Ctr)))]


# Networks' inference

library(minet)
library(igraph)

# Tadpole under treatment 1
mim <- build.mim(FTad_T1,estimator="spearman")
aracne_mat <- aracne(mim)
FTad_T1Net<-graph.adjacency(aracne_mat)
FTad_T1Net<-as.undirected(FTad_T1Net)

# Tadpole under treatment 2
mim <- build.mim(FTad_T2,estimator="spearman")
aracne_mat <- aracne(mim)
FTad_T2Net<-graph.adjacency(aracne_mat)
FTad_T2Net<-as.undirected(FTad_T2Net)

# Tadpole control
mim <- build.mim(FTad_Ctr,estimator="spearman")
aracne_mat <- aracne(mim)
FTad_CtrNet<-graph.adjacency(aracne_mat)
FTad_CtrNet<-as.undirected(FTad_CtrNet)

# Metamorphic under treatment 1
mim <- build.mim(FMet_T1,estimator="spearman")
aracne_mat <- aracne(mim)
FMet_T1Net<-graph.adjacency(aracne_mat)
FMet_T1Net<-as.undirected(FMet_T1Net)

# Metamorphic under treatment 2
mim <- build.mim(FMet_T2,estimator="spearman")
aracne_mat <- aracne(mim)
FMet_T2Net<-graph.adjacency(aracne_mat)
FMet_T2Net<-as.undirected(FMet_T2Net)

# Metamorphic control
mim <- build.mim(FMet_Ctr,estimator="spearman")
aracne_mat <- aracne(mim)
FMet_CtrNet<-graph.adjacency(aracne_mat)
FMet_CtrNet<-as.undirected(FMet_CtrNet)

# Sub-adult under treatment 1
mim <- build.mim(FAdl_T1,estimator="spearman")
aracne_mat <- aracne(mim)
FAdl_T1Net<-graph.adjacency(aracne_mat)
FAdl_T1Net<-as.undirected(FAdl_T1Net)

# Sub-adult under treatment 2
mim <- build.mim(FAdl_T2,estimator="spearman")
aracne_mat <- aracne(mim)
FAdl_T2Net<-graph.adjacency(aracne_mat)
FAdl_T2Net<-as.undirected(FAdl_T2Net)

# Sub-adult control
mim <- build.mim(FAdl_Ctr,estimator="spearman")
aracne_mat <- aracne(mim)
FAdl_CtrNet<-graph.adjacency(aracne_mat)
FAdl_CtrNet<-as.undirected(FAdl_CtrNet)

# Node-color assignment function

v_colored<-function(g, T_table, g_tax, p_tax, g_colors){
  require(igraph)
  #Identificación de elementos únicos del taxón de mayor jerarquía
  unq<-unique(T_table[,g_tax])
  #Asignación de colores asociados a los elementos de g_tax en el objeto igraph
  for(i in 1:length(unq)){
    IDs<-which(unq[i] == T_table[,g_tax])
    t_names<-unique(T_table[p_tax][IDs,])
    vertex<-which(vertex.attributes(g)$name %in% t_names)
    V(g)[vertex]$color<-g_colors[i]
  }
  return(g)
}

# Colors
unq<-unique(tax_fungi[,"Class"])
unq<-unq[-c(which(unq == "unidentified"), which(is.na(unq)))]
library(viridis)
colors <- sample(viridis_pal()(length(unq)))

# Tadpole under treatment 1
FTad_T1Net<-v_colored(FTad_T1Net, tax_fungi, g_tax = "Phylum",
                         p_tax = "Genus", g_colors = colors)
plot(FTad_T1Net, vertex.label.color="black",
     vertex.color = vertex.attributes(FTad_T1Net)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Tadpole under treatment 1")

# Tadpole under treatment 2
FTad_T2Net<-v_colored(FTad_T2Net, tax_fungi, g_tax = "Phylum",
                      p_tax = "Genus", g_colors = colors)
plot(FTad_T2Net, vertex.label.color="black",
     vertex.color = vertex.attributes(FTad_T2Net)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Tadpole under treatment 2")

# Tadpole control
FTad_CtrNet<-v_colored(FTad_CtrNet, tax_fungi, g_tax = "Phylum",
                      p_tax = "Genus", g_colors = colors)
plot(FTad_CtrNet, vertex.label.color="black",
     vertex.color = vertex.attributes(FTad_CtrNet)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Tadpole control")

# Metamorph under treatment 1
FMet_T1Net<-v_colored(FMet_T1Net, tax_fungi, g_tax = "Phylum",
                      p_tax = "Genus", g_colors = colors)
plot(FMet_T1Net, vertex.label.color="black",
     vertex.color = vertex.attributes(FMet_T1Net)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Metamorph under treatment 1")

# Metamorph under treatment 2
FMet_T2Net<-v_colored(FMet_T2Net, tax_fungi, g_tax = "Phylum",
                      p_tax = "Genus", g_colors = colors)
plot(FMet_T2Net, vertex.label.color="black",
     vertex.color = vertex.attributes(FMet_T2Net)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Metamorph under treatment 2")

# Metamorph control
FMet_CtrNet<-v_colored(FMet_CtrNet, tax_fungi, g_tax = "Phylum",
                      p_tax = "Genus", g_colors = colors)
plot(FMet_CtrNet, vertex.label.color="black",
     vertex.color = vertex.attributes(FMet_CtrNet)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Metamorph control")

# Sub-adult under treatment 1
FAdl_T1Net<-v_colored(FAdl_T1Net, tax_fungi, g_tax = "Phylum",
                      p_tax = "Genus", g_colors = colors)
plot(FAdl_T1Net, vertex.label.color="black",
     vertex.color = vertex.attributes(FAdl_T1Net)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Sub-adult under treatment 1")

# Sub-adult under treatment 2
FAdl_T2Net<-v_colored(FAdl_T2Net, tax_fungi, g_tax = "Phylum",
                      p_tax = "Genus", g_colors = colors)
plot(FAdl_T2Net, vertex.label.color="black",
     vertex.color = vertex.attributes(FAdl_T2Net)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Sub-adult under treatment 2")

# Sub-adult control
FAdl_CtrNet<-v_colored(FAdl_CtrNet, tax_fungi, g_tax = "Class",
                      p_tax = "Genus", g_colors = colors)
plot(FAdl_CtrNet, vertex.label.color="black",
     vertex.color = vertex.attributes(FAdl_CtrNet)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Sub-adult control")

# Node's Phyla function

TaxGroup<-function(g, T_table, g_tax, p_tax){
  require(igraph)
  unq<-unique(T_table[,g_tax])
  for(i in 1:length(unq)){
    IDs<-which(unq[i] == T_table[,g_tax])
    t_names<-unique(T_table[p_tax][IDs,])
    vertex<-which(vertex.attributes(g)$name %in% t_names)
    V(g)[vertex]$Taxon<-unq[i]
  }
  return(g)
}

FTad_CtrNet<-TaxGroup(FTad_CtrNet, tax_fungi, "Phylum", "Genus")
FTad_T1Net<-TaxGroup(FTad_T1Net, tax_fungi, "Phylum", "Genus")
FTad_T2Net<-TaxGroup(FTad_T2Net, tax_fungi, "Phylum", "Genus")
FMet_CtrNet<-TaxGroup(FMet_CtrNet, tax_fungi, "Phylum", "Genus")
FMet_T1Net<-TaxGroup(FMet_T1Net, tax_fungi, "Phylum", "Genus")
FMet_T2Net<-TaxGroup(FMet_T2Net, tax_fungi, "Phylum", "Genus")
FAdl_CtrNet<-TaxGroup(FAdl_CtrNet, tax_fungi, "Phylum", "Genus")
FAdl_T1Net<-TaxGroup(FAdl_T1Net, tax_fungi, "Phylum", "Genus")
FAdl_T2Net<-TaxGroup(FAdl_T2Net, tax_fungi, "Phylum", "Genus")

# Multilayer networks

library(muxViz)

ml_FTad<-list(FTad_T2Net, FTad_T1Net, FTad_CtrNet) # Tadpole
ml_FMet<-list(FMet_T2Net, FMet_T1Net, FMet_CtrNet) # Metamorphic
ml_FAdl<-list(FAdl_T2Net, FAdl_T1Net, FAdl_CtrNet) # Adult

# Matrix vertex-color function for multilayer network
mat_colors<-function(colors, g.list){
  colmat<-matrix(colors, nrow = length(colors),
                 ncol = length(g.list))
  return(colmat)
}
  
# 3D plot for tadpoles
lay <- layoutMultiplex(ml_FTad, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_FTad, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=c("Treatment 2", "Treatment 1", "Control"),
                 layer.labels.cex=1.5, node.size.values=10,
                 node.size.scale=0.6,
                 node.colors=mat_colors(V(ml_FTad[[1]])$color, ml_FTad),
                 edge.colors="black",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# 3D plot for metamorphics
lay <- layoutMultiplex(ml_FMet, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_FMet, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=c("Treatment 2", "Treatment 1", "Control"), layer.labels.cex=1.5,
                 node.size.values=10, node.size.scale=0.6,
                 node.colors=mat_colors(V(FMet_T1Net)$color, ml_FMet),
                 edge.colors="black",
                 node.colors.aggr=mat_colors(V(FMet_T1Net)$color, ml_FMet),
                 show.aggregate=F)

# 3D plot for the adults
lay <- layoutMultiplex(ml_FAdl, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_FAdl, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=c("Treatment 2", "Treatment 1", "Control"), layer.labels.cex=1.5,
                 node.size.values=10, node.size.scale=0.6,
                 node.colors=mat_colors(V(FAdl_T1Net)$color, ml_FAdl),
                 edge.colors="black",
                 node.colors.aggr=mat_colors(V(FAdl_T1Net)$color, ml_FAdl),
                 show.aggregate=F)


# Centrality analysis

# Centrality function
ctr<-function(g.list, ctr_type){
  require(igraph)
  if(ctr_type == "degree"){
    for(i in 1:length(g.list)){
      ctr_max<-which.max(degree(g.list[[i]]))
      ctr<-degree(g.list[[i]])/degree(g.list[[i]])[ctr_max]
      Redpal<-heat.colors(5, alpha=1)
      hl1<-which(ctr > 0 & ctr <= 0.2)
      V(g.list[[i]])[hl1]$hl<-Redpal[5]
      hl2<-which(ctr > 0.2 & ctr <= 0.4)
      V(g.list[[i]])[hl2]$hl<-Redpal[4]
      hl3<-which(ctr > 0.4 & ctr <= 0.6)
      V(g.list[[i]])[hl3]$hl<-Redpal[3]
      hl4<-which(ctr > 0.6 & ctr <= 0.8)
      V(g.list[[i]])[hl4]$hl<-Redpal[2]
      hl5<-which(ctr > 0.8 & ctr <= 1)
      V(g.list[[i]])[hl5]$hl<-Redpal[1]
      vertex.attributes(g.list[[i]])$hl
    }
  }
  if(ctr_type == "betweenness"){
    for(i in 1:length(g.list)){
      ctr_max<-which.max(betweenness(g.list[[i]]))
      ctr<-betweenness(g.list[[i]])/betweenness(g.list[[i]])[ctr_max]
      Redpal<-heat.colors(5, alpha=1)
      hl1<-which(ctr > 0 & ctr <= 0.2)
      V(g.list[[i]])[hl1]$hl<-Redpal[5]
      hl2<-which(ctr > 0.2 & ctr <= 0.4)
      V(g.list[[i]])[hl2]$hl<-Redpal[4]
      hl3<-which(ctr > 0.4 & ctr <= 0.6)
      V(g.list[[i]])[hl3]$hl<-Redpal[3]
      hl4<-which(ctr > 0.6 & ctr <= 0.8)
      V(g.list[[i]])[hl4]$hl<-Redpal[2]
      hl5<-which(ctr > 0.8 & ctr <= 1)
      V(g.list[[i]])[hl5]$hl<-Redpal[1]
      vertex.attributes(g.list[[i]])$hl
    }
  }
  if(ctr_type == "closeness"){
    for(i in 1:length(g.list)){
      ctr_max<-which.max(closeness(g.list[[i]]))
      ctr<-closeness(g.list[[i]])/closeness(g.list[[i]])[ctr_max]
      Redpal<-heat.colors(5, alpha=1)
      hl1<-which(ctr > 0 & ctr <= 0.2)
      V(g.list[[i]])[hl1]$hl<-Redpal[5]
      hl2<-which(ctr > 0.2 & ctr <= 0.4)
      V(g.list[[i]])[hl2]$hl<-Redpal[4]
      hl3<-which(ctr > 0.4 & ctr <= 0.6)
      V(g.list[[i]])[hl3]$hl<-Redpal[3]
      hl4<-which(ctr > 0.6 & ctr <= 0.8)
      V(g.list[[i]])[hl4]$hl<-Redpal[2]
      hl5<-which(ctr > 0.8 & ctr <= 1)
      V(g.list[[i]])[hl5]$hl<-Redpal[1]
      vertex.attributes(g.list[[i]])$hl
    }
  }
  return(g.list)
}

# Degree centrality
ml_FTad<-ctr(ml_FTad, "degree")
ml_FMet<-ctr(ml_FMet, "degree")
ml_FAdl<-ctr(ml_FAdl, "degree")

colmat<-matrix(c(V(ml_FTad[[1]])$hl,V(ml_FTad[[2]])$hl,V(ml_FTad[[3]])$hl),
               nrow = length(V(ml_FTad[[1]])), ncol = length(ml_FTad))

# Degree 3D plot for tadpoles
lay <- layoutMultiplex(ml_FTad, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_FTad, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=NULL, layer.labels.cex=1.5,
                 node.size.values="auto", node.size.scale=abs_mat(a, ml_FTad, 20),
                 node.colors=mat_colors(V(ml_FTad[[1]])$hl, ml_FTad),
                 edge.colors="black",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Degree 3D plot for metamorphic
lay <- layoutMultiplex(ml_FMet, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_FMet, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=c("Treatment 2", "Treatment 1", "Control"), layer.labels.cex=1.5,
                 node.size.values=10, node.size.scale=0.6,
                 node.colors=mat_colors(V(ml_FMet[[1]])$hl, ml_FMet),
                 edge.colors="black",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Degree 3D plot for adults
lay <- layoutMultiplex(ml_FAdl, layout="kk", ggplot.format=F, box=T)
plot_multiplex3D(ml_FAdl, layer.layout=lay,
                 layer.colors=c("red3", "orange", "green3"),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=c("Treatment 2", "Treatment 1", "Control"), layer.labels.cex=1.5,
                 node.size.values=10, node.size.scale=0.6,
                 node.colors=mat_colors(V(ml_FAdl[[1]])$hl, ml_FAdl),
                 edge.colors="black",
                 node.colors.aggr=NULL,
                 show.aggregate=F)


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

# Diversity between treatments

# Simpson dominance for tadpoles
T1_s<-ab_table_div(t(FTad_T1), "simpson")
T2_s<-ab_table_div(t(FTad_T2), "simpson")
T3_s<-ab_table_div(t(FTad_Ctr), "simpson")

T1df_s<-data.frame("Period" = rep("Treatment 1", length(T1_s)),
                   "Data" = T1_s)
T2df_s<-data.frame("Period" = rep("Treatment 2", length(T2_s)),
                   "Data" = T2_s)
T3df_s<-data.frame("Period" = rep("Control", length(T3_s)),
                   "Data" = T3_s)

Tdf_s<-rbind(T1df_s, T2df_s, T3df_s)

# Violin plot
library(ggplot2)
ggplot(Tdf_s, aes(x=Tdf_s$Period, y=Tdf_s$Data)) +
  geom_violin(fill="darkslategray3", color="black") +
  geom_boxplot(width=0.15, notch=TRUE,
               fill=c("green3", "red3",
                      "darkorange1"),
               color="black") + ggtitle("Simpson dominance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + ylab("")


# Simpson dominance for metamorphic
M1_s<-ab_table_div(t(FMet_T1), "simpson")
M2_s<-ab_table_div(t(FMet_T2), "simpson")
M3_s<-ab_table_div(t(FMet_Ctr), "simpson")

M1df_s<-data.frame("Period" = rep("Treatment 1", length(M1_s)),
                   "Data" = M1_s)
M2df_s<-data.frame("Period" = rep("Treatment 2", length(M2_s)),
                   "Data" = M2_s)
M3df_s<-data.frame("Period" = rep("Control", length(M3_s)),
                   "Data" = M3_s)

Mdf_s<-rbind(M1df_s, M2df_s, M3df_s)

# Violin plot
library(ggplot2)
ggplot(Mdf_s, aes(x=Mdf_s$Period, y=Mdf_s$Data)) +
  geom_violin(fill="darkslategray3", color="black") +
  geom_boxplot(width=0.15, notch=TRUE,
               fill=c("green3", "red3",
                      "darkorange1"),
               color="black") + ggtitle("Simpson dominance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + ylab("")

# Simpson dominance for adults
A1_s<-ab_table_div(t(FAdl_T1), "simpson")
A2_s<-ab_table_div(t(FAdl_T2), "simpson")
A3_s<-ab_table_div(t(FAdl_Ctr), "simpson")

A1df_s<-data.frame("Period" = rep("Treatment 1", length(A1_s)),
                   "Data" = A1_s)
A2df_s<-data.frame("Period" = rep("Treatment 2", length(A2_s)),
                   "Data" = A2_s)
A3df_s<-data.frame("Period" = rep("Control", length(A3_s)),
                   "Data" = A3_s)

Adf_s<-rbind(A1df_s, A2df_s, A3df_s)

# Violin plot
library(ggplot2)
ggplot(Adf_s, aes(x=Adf_s$Period, y=Adf_s$Data)) +
  geom_violin(fill="darkslategray3", color="black") +
  geom_boxplot(width=0.15, notch=TRUE,
               fill=c("green3", "red3",
                      "darkorange1"),
               color="black") + ggtitle("Simpson dominance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("") + ylab("")



# Rank-abundance analysis


# Tadpoles

library(RADanalysis)
sample_classes <- c(rep(1, length(T1_s)),rep(2, length(T2_s)),
                    rep(3, length(T3_s)))
line_cols <- c("darkorange1","red3","green3")

FTad_mat<-rbind(FTad_T1, FTad_T2, FTad_Ctr)

# Normalized abundances
n_FTad_mat<-FTad_mat/rowSums(FTad_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_FTad_mat)){
  sorted_abs[[i]]<-sort(n_FTad_mat[i,], decreasing = T)
}
n_FTad_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                     ncol=length(sorted_abs[[1]]), byrow=TRUE)


# Plot the axis
plot(1e10,xlim = c(1,30),ylim = c(0,0.6),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15,20,25,30))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4,0.5,0.6),las = 0)

# Rank-abundance colors per treatment
a <- representative_RAD(norm_rad = n_FTad_mat,sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_FTad_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_FTad_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Treatment 1","Treatment 2", "Control"),
       col = line_cols, lwd = 3)


# Metamorphic

sample_classes <- c(rep(1, length(M1_s)),rep(2, length(M2_s)),
                    rep(3, length(M3_s)))
line_cols <- c("darkorange1","red3","green3")

FMet_mat<-rbind(FMet_T1, FMet_T2, FMet_Ctr)

# Normalized abundances
n_FMet_mat<-FMet_mat/rowSums(FMet_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_FMet_mat)){
  sorted_abs[[i]]<-sort(n_FMet_mat[i,], decreasing = T)
}
n_FMet_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                   ncol=length(sorted_abs[[1]]), byrow=TRUE)


# Plot the axis
plot(1e10,xlim = c(1,30),ylim = c(0,0.6),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15,20,25,30))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4,0.5,0.6),las = 0)

# Rank-abundance colors per treatment
a <- representative_RAD(norm_rad = n_FMet_mat,sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_FMet_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_FMet_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Treatment 1","Treatment 2", "Control"),
       col = line_cols, lwd = 3)


# Adults

sample_classes <- c(rep(1, length(A1_s)),rep(2, length(A2_s)),
                    rep(3, length(A3_s)))
line_cols <- c("darkorange1","red3","green3")

FAdl_mat<-rbind(FAdl_T1, FAdl_T2, FAdl_Ctr)

# Normalized abundances
n_FAdl_mat<-FAdl_mat/rowSums(FAdl_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_FAdl_mat)){
  sorted_abs[[i]]<-sort(n_FAdl_mat[i,], decreasing = T)
}
n_FAdl_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                   ncol=length(sorted_abs[[1]]), byrow=TRUE)


# Plot the axis
plot(1e10,xlim = c(1,30),ylim = c(0,0.6),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15,20,25,30))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4,0.5,0.6),las = 0)

# Rank-abundance colors per treatment
a <- representative_RAD(norm_rad = n_FAdl_mat,sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_FAdl_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_FAdl_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Treatment 1","Treatment 2", "Control"),
       col = line_cols, lwd = 3)






# MSD analysis


# Tadpoles
d <- dist(x = n_FTad_mat,method = "manhattan")
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
legend("bottomright",bty = "n",legend = c("Treatment 1","Treatment 2","Control"),
       col = line_cols,pch = 19)


# Metamorphic
d <- dist(x = n_FMet_mat,method = "manhattan")
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
legend("bottomright",bty = "n",legend = c("Treatment 1","Treatment 2","Control"),
       col = line_cols,pch = 19)

# Adults
d <- dist(x = n_FAdl_mat,method = "manhattan")
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
legend("bottomright",bty = "n",legend = c("Treatment 1","Treatment 2","Control"),
       col = line_cols,pch = 19)






# Dysbiosis analysis

# Load required packages
library(dysbiosisR)
library(ggplot2)
library(microbiome)
library(dplyr)
library(phyloseq)

# Transform the data into a phyloseq object
tax_fungi<-as.matrix(hongos[,t1:t2])
otus_fungi<-as.matrix(hongos[,-(t1:t2)])
otus_fungi<-otus_fungi[,-1]
# OTUs IDs as row names and sample IDs as column names
ID_otus <- otus_fungi[,1]
otus_fungi<-otus_fungi[,-1]
samples_fungi<-colnames(otus_fungi)
otus_fungi<-matrix(as.numeric(otus_fungi),
                   ncol = length(samples_fungi),
                   nrow = length(ID_otus))
colnames(otus_fungi)<-samples_fungi
rownames(otus_fungi)<-ID_otus
# OTUs IDs as row names in the taxa table
rownames(tax_fungi)<-ID_otus
# Phyloseq objects
otus_fungi<-otu_table(otus_fungi, taxa_are_rows = TRUE)
tax_fungi<-tax_table(tax_fungi)
physeq = phyloseq(otus_fungi, tax_fungi)
# Sample metadata
rownames(meta_hongos)<-meta_hongos[,1]
meta_hongos<-meta_hongos[,-1]
sample_data<-sample_data(meta_hongos)
# Merge of phyloseq objects
physeq2<-merge_phyloseq(physeq, sample_data)

# Dysbiosis analysis

# Bray-Curtis distance matrix
dist.mat <- phyloseq::distance(physeq2, "bray")
# Get reference samples
ref.samples <- sample_names(subset_samples(physeq2, 
                                           Treatment == "Control"))
# Community level variation analysis
dysbiosis_1 <- dysbiosisMedianCLV(physeq2,
                                  dist_mat = dist.mat,
                                  reference_samples = ref.samples)
# We sample the data set identifying as dysbiotic the data under the 90th percentile
dysbiosis_thres <- quantile(subset(dysbiosis_1, Treatment == "Treatment2")$score, 0.9)
normobiosis_thres <- quantile(subset(dysbiosis_1, Treatment == "Treatment2")$score, 0.1)

dysbiosis_1 <- dysbiosis_1 |> 
  mutate(isDysbiostic = ifelse(score >= dysbiosis_thres, TRUE, FALSE))


# Dysbiosis plot measures according to CLV method
p1 <- plotDysbiosis(df=dysbiosis_1,
                    xvar="Treatment",
                    yvar="score",
                    colors=c(Treatment1="orange", Treatment2="red",
                             Control="green"),
                    show_points = FALSE) +
  labs(x="", y="Dysbiosis Score") +
  theme_bw(base_size = 14)
p1


# Dysbiosis plot measures according to euclidean method
dysbiosis_2 <- euclideanDistCentroids(physeq2,
                                      dist_mat = dist.mat,
                                      use_squared = TRUE,
                                      group_col = "Treatment",
                                      control_label = "Control",
                                      case_label = "Treatment2")

p2 <- plotDysbiosis(df=dysbiosis_2,
                    xvar="Treatment",
                    yvar="CentroidDist_score",
                    colors=c(Treatment1="orange", Treatment2="red",
                             Control="green"),
                    show_points = FALSE) +
  labs(x="", y="Dysbiosis Score (Centroid)") +
  theme_bw(base_size = 14)
p2


# Dysbiosis plot measures according to the Combined alpha-beta diversity based score
dysbiosis_3 <- combinedShannonJSD(physeq2,
                                  reference_samples = ref.samples)


p3 <- plotDysbiosis(df=dysbiosis_3,
                    xvar="Treatment",
                    yvar="ShannonJSDScore",
                    colors=c(Treatment1="orange", Treatment2="red",
                             Control="green"),
                    show_points = FALSE)+
  labs(x="", y="Shannon-JSD\nDysbiosis Score") +
  theme_bw(base_size = 14)
p3


cloud.results <- cloudStatistic(physeq2,
                                dist_mat = dist.mat,
                                reference_samples = ref.samples,
                                ndim=-1,
                                k_num=5)

p4 <- plotDysbiosis(df=cloud.results,
                    xvar="Treatment",
                    yvar="log2Stats",
                    colors=c(Treatment1="red", Treatment2="orange",
                             Control="green"),
                    show_points = FALSE) +
  labs(x="", y="Dysbiosis CLOUD Score")
p4

# Degree analysis

degree_df <- data.frame(Species = vertex.attributes(FTad_CtrNet)$name,
                 color = vertex.attributes(FTad_CtrNet)$color,
                 Phylum = vertex.attributes(FTad_CtrNet)$Taxon,
                 Ctr_degree = degree(FTad_CtrNet),
                 T1_degree = degree(FTad_T1Net),
                 T2_degree = degree(FTad_T2Net))

degree_df<-degree_df[-which(rowSums(degree_df[,c(4,5,6)])/
                              mean(rowSums(degree_df[,c(4,5,6)]))<0.9),]

library(gridExtra)
p1<-ggplot(data=degree_df, aes(x=Species, y=Ctr_degree, fill=Phylum)) +
  geom_bar(stat="identity")
p2<-ggplot(data=degree_df, aes(x=Species, y=T1_degree, fill=Phylum)) +
  geom_bar(stat="identity")
p3<-ggplot(data=degree_df, aes(x=Species, y=T2_degree, fill=Phylum)) +
  geom_bar(stat="identity")

grid.arrange(p1 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.position = "none") + ylab("Control"),
             p2 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.position = "none") + ylab("Treatment 1"),
             p3 + coord_flip() +
               theme(axis.text.y = element_text(size = 0.01),
                     legend.text = element_text(size = 6),
                     legend.position="bottom", legend.box = "horizontal") +
               ylab("Treatment 2"),
             ncol=3)


