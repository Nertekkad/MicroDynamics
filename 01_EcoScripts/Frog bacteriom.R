# Data of frog's microbiota under different temperatures

setwd("C:/Users/LaV_V/Downloads")
# Data loading
library(readxl)
meta_bacterias<-read.csv("Metadata_Bacterias.csv")
bacterias<-read_excel("Bacterias.xlsx")

# Identify frog´s life stages
life_stage<-unique(meta_bacterias$Life_stage)

# Separation of data by life stages type
Tadpole<-meta_bacterias[which(meta_bacterias$Life_stage == life_stage[1]),]
Metamorphic<-meta_bacterias[which(meta_bacterias$Life_stage == life_stage[2]),]
Sub_adult<-meta_bacterias[which(meta_bacterias$Life_stage == life_stage[3]),]

# Tadpole

# Treatment identification
Tad_T1<-Tadpole$sample.id[which(Tadpole$Treatment=="Treatment1")]
Tad_T2<-Tadpole$sample.id[which(Tadpole$Treatment=="Treatment2")]
Tad_Ctr<-Tadpole$sample.id[which(Tadpole$Treatment=="Control")]
# Separation of data by treatment types
BaTad_T1<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Tad_T1)])
BaTad_T2<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Tad_T2)])
BaTad_Ctr<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Tad_Ctr)])

# Metamorphic

# Treatment identification
Met_T1<-Metamorphic$sample.id[which(Metamorphic$Treatment=="Treatment1")]
Met_T2<-Metamorphic$sample.id[which(Metamorphic$Treatment=="Treatment2")]
Met_Ctr<-Metamorphic$sample.id[which(Metamorphic$Treatment=="Control")]
# Separation of data by treatment types
BaMet_T1<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Met_T1)])
BaMet_T2<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Met_T2)])
BaMet_Ctr<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Met_Ctr)])

# Sub-adult

# Treatment identification
Adl_T1<-Sub_adult$sample.id[which(Sub_adult$Treatment=="Treatment1")]
Adl_T2<-Sub_adult$sample.id[which(Sub_adult$Treatment=="Treatment2")]
Adl_Ctr<-Sub_adult$sample.id[which(Sub_adult$Treatment=="Control")]
# Separation of data by treatment types
BaAdl_T1<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Adl_T1)])
BaAdl_T2<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Adl_T2)])
BaAdl_Ctr<-as.data.frame(bacterias[,which(colnames(bacterias) %in% Adl_Ctr)])

# Taxa table
t1<-which(colnames(bacterias)=="Kingdom")
t2<-which(colnames(bacterias)=="Species")
tax_bacter<-as.data.frame(bacterias[,t1:t2])

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
BTad_T1<-T_collapse(tax_bacter, BaTad_T1, "Genus")
BTad_T1<-BTad_T1[, -which(colnames(BTad_T1) == "uncultured")]
BTad_T1<-BTad_T1[, -which(is.na(colnames(BTad_T1)))]

# Tadpole under treatment 2
BTad_T2<-T_collapse(tax_bacter, BaTad_T2, "Genus")
BTad_T2<-BTad_T2[, -which(colnames(BTad_T2) == "uncultured")]
BTad_T2<-BTad_T2[, -which(is.na(colnames(BTad_T2)))]

# Tadpole control
BTad_Ctr<-T_collapse(tax_bacter, BaTad_Ctr, "Genus")
BTad_Ctr<-BTad_Ctr[, -which(colnames(BTad_Ctr) == "uncultured")]
BTad_Ctr<-BTad_Ctr[, -which(is.na(colnames(BTad_Ctr)))]

# Metamorphic under treatment 1
BMet_T1<-T_collapse(tax_bacter, BaMet_T1, "Genus")
BMet_T1<-BMet_T1[, -which(colnames(BMet_T1) == "uncultured")]
BMet_T1<-BMet_T1[, -which(is.na(colnames(BMet_T1)))]

# Metamorphic under treatment 2
BMet_T2<-T_collapse(tax_bacter, BaMet_T2, "Genus")
BMet_T2<-BMet_T2[, -which(colnames(BMet_T2) == "uncultured")]
BMet_T2<-BMet_T2[, -which(is.na(colnames(BMet_T2)))]

# Metamorphic control
BMet_Ctr<-T_collapse(tax_bacter, BaMet_Ctr, "Genus")
BMet_Ctr<-BMet_Ctr[, -which(colnames(BMet_Ctr) == "uncultured")]
BMet_Ctr<-BMet_Ctr[, -which(is.na(colnames(BMet_Ctr)))]

# Sub-adult under treatment 1
BAdl_T1<-T_collapse(tax_bacter, BaAdl_T1, "Genus")
BAdl_T1<-BAdl_T1[, -which(colnames(BAdl_T1) == "uncultured")]
BAdl_T1<-BAdl_T1[, -which(is.na(colnames(BAdl_T1)))]

# Sub-adult under treatment 2
BAdl_T2<-T_collapse(tax_bacter, BaAdl_T2, "Genus")
BAdl_T2<-BAdl_T2[, -which(colnames(BAdl_T2) == "uncultured")]
BAdl_T2<-BAdl_T2[, -which(is.na(colnames(BAdl_T2)))]

# Sub-adult control
BAdl_Ctr<-T_collapse(tax_bacter, BaAdl_Ctr, "Genus")
BAdl_Ctr<-BAdl_Ctr[, -which(colnames(BAdl_Ctr) == "uncultured")]
BAdl_Ctr<-BAdl_Ctr[, -which(is.na(colnames(BAdl_Ctr)))]

# Networks' inference

library(minet)
library(igraph)

# Tadpole under treatment 1
mim <- build.mim(BTad_T1,estimator="spearman")
aracne_mat <- aracne(mim)
BTad_T1Net<-graph.adjacency(aracne_mat)
BTad_T1Net<-as.undirected(BTad_T1Net)

# Tadpole under treatment 2
mim <- build.mim(BTad_T2,estimator="spearman")
aracne_mat <- aracne(mim)
BTad_T2Net<-graph.adjacency(aracne_mat)
BTad_T2Net<-as.undirected(BTad_T2Net)

# Tadpole control
mim <- build.mim(BTad_Ctr,estimator="spearman")
aracne_mat <- aracne(mim)
BTad_CtrNet<-graph.adjacency(aracne_mat)
BTad_CtrNet<-as.undirected(BTad_CtrNet)

# Metamorphic under treatment 1
mim <- build.mim(BMet_T1,estimator="spearman")
aracne_mat <- aracne(mim)
BMet_T1Net<-graph.adjacency(aracne_mat)
BMet_T1Net<-as.undirected(BMet_T1Net)

# Metamorphic under treatment 2
mim <- build.mim(BMet_T2,estimator="spearman")
aracne_mat <- aracne(mim)
BMet_T2Net<-graph.adjacency(aracne_mat)
BMet_T2Net<-as.undirected(BMet_T2Net)

# Metamorphic control
mim <- build.mim(BMet_Ctr,estimator="spearman")
aracne_mat <- aracne(mim)
BMet_CtrNet<-graph.adjacency(aracne_mat)
BMet_CtrNet<-as.undirected(BMet_CtrNet)

# Sub-adult under treatment 1
mim <- build.mim(BAdl_T1,estimator="spearman")
aracne_mat <- aracne(mim)
BAdl_T1Net<-graph.adjacency(aracne_mat)
BAdl_T1Net<-as.undirected(BAdl_T1Net)

# Sub-adult under treatment 2
mim <- build.mim(BAdl_T2,estimator="spearman")
aracne_mat <- aracne(mim)
BAdl_T2Net<-graph.adjacency(aracne_mat)
BAdl_T2Net<-as.undirected(BAdl_T2Net)

# Sub-adult control
mim <- build.mim(BAdl_Ctr,estimator="spearman")
aracne_mat <- aracne(mim)
BAdl_CtrNet<-graph.adjacency(aracne_mat)
BAdl_CtrNet<-as.undirected(BAdl_CtrNet)

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
unq<-unique(tax_bacter[,"Phylum"])
unq<-unq[-c(which(unq == "uncultured"), which(is.na(unq)))]
library(viridis)
colors <- sample(viridis_pal()(length(unq)))

# Tadpole under treatment 1
BTad_T1Net<-v_colored(BTad_T1Net, tax_bacter, g_tax = "Phylum",
                      p_tax = "Genus", g_colors = colors)
plot(BTad_T1Net, vertex.label.color="black",
     vertex.color = vertex.attributes(BTad_T1Net)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Tadpole under treatment 1")

# Tadpole under treatment 2
BTad_T2Net<-v_colored(BTad_T2Net, tax_bacter, g_tax = "Phylum",
                      p_tax = "Genus", g_colors = colors)
plot(BTad_T2Net, vertex.label.color="black",
     vertex.color = vertex.attributes(BTad_T2Net)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Tadpole under treatment 2")

# Tadpole control
BTad_CtrNet<-v_colored(BTad_CtrNet, tax_bacter, g_tax = "Phylum",
                       p_tax = "Genus", g_colors = colors)
plot(BTad_CtrNet, vertex.label.color="black",
     vertex.color = vertex.attributes(BTad_CtrNet)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Tadpole control")

# Metamorphic under treatment 1
BMet_T1Net<-v_colored(BMet_T1Net, tax_bacter, g_tax = "Phylum",
                      p_tax = "Genus", g_colors = colors)
plot(BMet_T1Net, vertex.label.color="black",
     vertex.color = vertex.attributes(BMet_T1Net)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Metamorphic under treatment 1")

# Metamorphic under treatment 2
BMet_T2Net<-v_colored(BMet_T2Net, tax_bacter, g_tax = "Phylum",
                      p_tax = "Genus", g_colors = colors)
plot(BMet_T2Net, vertex.label.color="black",
     vertex.color = vertex.attributes(BMet_T2Net)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Metamorphic under treatment 2")

# Metamorphic control
BMet_CtrNet<-v_colored(BMet_CtrNet, tax_bacter, g_tax = "Phylum",
                       p_tax = "Genus", g_colors = colors)
plot(BMet_CtrNet, vertex.label.color="black",
     vertex.color = vertex.attributes(BMet_CtrNet)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Metamorphic control")

# Sub-adult under treatment 1
BAdl_T1Net<-v_colored(BAdl_T1Net, tax_bacter, g_tax = "Phylum",
                      p_tax = "Genus", g_colors = colors)
plot(BAdl_T1Net, vertex.label.color="black",
     vertex.color = vertex.attributes(BAdl_T1Net)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Sub-adult under treatment 1")

# Sub-adult under treatment 2
BAdl_T2Net<-v_colored(BAdl_T2Net, tax_bacter, g_tax = "Phylum",
                      p_tax = "Genus", g_colors = colors)
plot(BAdl_T2Net, vertex.label.color="black",
     vertex.color = vertex.attributes(BAdl_T2Net)$color, vertex.label.cex=.5,
     vertex.label.dist=1,layout=layout_with_kk, vertex.size = 5,
     main = "Sub-adult under treatment 2")

# Sub-adult control
BAdl_CtrNet<-v_colored(BAdl_CtrNet, tax_bacter, g_tax = "Phylum",
                       p_tax = "Genus", g_colors = colors)
plot(BAdl_CtrNet, vertex.label.color="black",
     vertex.color = vertex.attributes(BAdl_CtrNet)$color, vertex.label.cex=.5,
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

BTad_CtrNet<-TaxGroup(BTad_CtrNet, tax_bacter, "Phylum", "Genus")
BTad_T1Net<-TaxGroup(BTad_T1Net, tax_bacter, "Phylum", "Genus")
BTad_T2Net<-TaxGroup(BTad_T2Net, tax_bacter, "Phylum", "Genus")
BMet_CtrNet<-TaxGroup(BMet_CtrNet, tax_bacter, "Phylum", "Genus")
BMet_T1Net<-TaxGroup(BMet_T1Net, tax_bacter, "Phylum", "Genus")
BMet_T2Net<-TaxGroup(BMet_T2Net, tax_bacter, "Phylum", "Genus")
BAdl_CtrNet<-TaxGroup(BAdl_CtrNet, tax_bacter, "Phylum", "Genus")
BAdl_T1Net<-TaxGroup(BAdl_T1Net, tax_bacter, "Phylum", "Genus")
BAdl_T2Net<-TaxGroup(BAdl_T2Net, tax_bacter, "Phylum", "Genus")


# Multilayer networks

library(muxViz)

ml_BTad<-list(BTad_T2Net, BTad_T1Net, BTad_CtrNet) # Tadpole
ml_BMet<-list(BMet_T2Net, BMet_T1Net, BMet_CtrNet) # Metamorphic
ml_BAdl<-list(BAdl_T2Net, BAdl_T1Net, BAdl_CtrNet) # Adult

# Node-color matrix function
node_color_mat<-function(g.list, type){
  if(type == "phylo"){
    # Colors for each node
    taxacolors<-list()
    for(i in 1:length(g.list)){
      taxacolors[[i]]<-V(g.list[[i]])$color
    }
    # Input color vector
    taxacolors<-unlist(taxacolors)
    # Color matrix at cetain taxonomic level
    taxamat<-matrix(taxacolors, nrow = length(V(g.list[[1]])),
                    ncol = length(g.list))
    return(taxamat)
  } else
    if(type == "centrality"){
      # Colors for each node
      centrality<-list()
      for(i in 1:length(g.list)){
        centrality[[i]]<-V(g.list[[i]])$hl
      }
      # Input color vector
      centrality<-unlist(centrality)
      # Color matrix at cetain taxonomic level
      heatmat<-matrix(centrality, nrow = length(V(g.list[[1]])),
                      ncol = length(g.list))
      return(heatmat)
    }
}

# Abundances matrix function
abs_mat<-function(abs.list, g.list, n){
  # Colors for each node
  abundances<-list()
  for(i in 1:length(abs.list)){
    abundances[[i]]<-log(colSums(abs.list[[i]])+2)/n
  }
  # Input color vector
  abundances<-as.vector(unlist(abundances))
  # Color matrix at cetain taxonomic level
  absmat<-matrix(abundances, nrow = length(V(g.list[[1]])),
                 ncol = length(g.list))
  return(absmat)
}

# Abundance tables list
abs_BTad<-list(BTad_T2, BTad_T1, BTad_Ctr)

# 3D plot for tadpoles
lay <- layoutMultiplex(ml_BTad, layout="fr", ggplot.format=F, box=T)
plot_multiplex3D(ml_BTad, layer.layout=lay,
                 layer.colors=rainbow(length(ml_BTad)),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=c("Treatment 2", "Treatment 1", "Control"), layer.labels.cex=1.5,
                 node.size.values="auto",
                 node.size.scale=abs_mat(abs_BTad, ml_BTad, 2),
                 node.colors=node_color_mat(ml_BTad, "phylo"),
                 edge.colors="black",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# 3D plot for metamorphics
lay <- layoutMultiplex(ml_BMet, layout="fr", ggplot.format=F, box=T)
plot_multiplex3D(ml_BMet, layer.layout=lay,
                 layer.colors=rainbow(length(ml_BMet)),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=c("Treatment 1", "Treatment 2", "Control"), layer.labels.cex=1.5,
                 node.size.values=10, node.size.scale=0.6,
                 node.colors=mat_colors(V(BMet_T1Net)$color, ml_BMet),
                 edge.colors="black",
                 node.colors.aggr=mat_colors(V(BMet_T1Net)$color, ml_BMet),
                 show.aggregate=F)

# 3D plot for the adults
lay <- layoutMultiplex(ml_BAdl, layout="fr", ggplot.format=F, box=T)
plot_multiplex3D(ml_BAdl, layer.layout=lay,
                 layer.colors=rainbow(length(ml_BAdl)),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=c("Treatment 1", "Treatment 2", "Control"), layer.labels.cex=1.5,
                 node.size.values=10, node.size.scale=0.6,
                 node.colors=mat_colors(V(BAdl_T1Net)$color, ml_BAdl),
                 edge.colors="black",
                 node.colors.aggr=mat_colors(V(BAdl_T1Net)$color, ml_BAdl),
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
ml_BTad<-ctr(ml_BTad, "degree")
ml_BMet<-ctr(ml_BMet, "degree")
ml_BAdl<-ctr(ml_BAdl, "degree")

# Degree 3D plot for tadpoles
lay <- layoutMultiplex(ml_BTad, layout="fr", ggplot.format=F, box=T)
plot_multiplex3D(ml_BTad, layer.layout=lay,
                 layer.colors=rainbow(length(ml_BTad)),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=c("Treatment 1", "Treatment 2", "Control"), layer.labels.cex=1.5,
                 node.size.values="auto", node.size.scale=abs_mat(),
                 node.colors=node_color_mat(ml_BTad, "phylo"),
                 edge.colors="black",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Degree 3D plot for metamorphic
lay <- layoutMultiplex(ml_BMet, layout="fr", ggplot.format=F, box=T)
plot_multiplex3D(ml_BMet, layer.layout=lay,
                 layer.colors=rainbow(length(ml_BMet)),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=c("Treatment 1", "Treatment 2", "Control"), layer.labels.cex=1.5,
                 node.size.values=10, node.size.scale=0.6,
                 node.colors=mat_colors(V(ml_BMet[[1]])$hl, ml_BMet),
                 edge.colors="black",
                 node.colors.aggr=NULL,
                 show.aggregate=F)

# Degree 3D plot for adults
lay <- layoutMultiplex(ml_BAdl, layout="fr", ggplot.format=F, box=T)
plot_multiplex3D(ml_BAdl, layer.layout=lay,
                 layer.colors=rainbow(length(ml_BAdl)),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=c("Treatment 1", "Treatment 2", "Control"), layer.labels.cex=1.5,
                 node.size.values=10, node.size.scale=0.6,
                 node.colors=mat_colors(V(ml_BAdl[[1]])$hl, ml_BAdl),
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
T1_s<-ab_table_div(t(BTad_T1), "simpson")
T2_s<-ab_table_div(t(BTad_T2), "simpson")
T3_s<-ab_table_div(t(BTad_Ctr), "simpson")

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
M1_s<-ab_table_div(t(BMet_T1), "simpson")
M2_s<-ab_table_div(t(BMet_T2), "simpson")
M3_s<-ab_table_div(t(BMet_Ctr), "simpson")

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
A1_s<-ab_table_div(t(BAdl_T1), "simpson")
A2_s<-ab_table_div(t(BAdl_T2), "simpson")
A3_s<-ab_table_div(t(BAdl_Ctr), "simpson")

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

# Tadpole
library(RADanalysis)
sample_classes <- c(rep(1, length(T1_s)),rep(2, length(T2_s)),
                    rep(3, length(T3_s)))
line_cols <- c("darkorange1","red3","green3")

BTad_mat<-rbind(BTad_T1, BTad_T2, BTad_Ctr)

# Normalized abundances
n_BTad_mat<-BTad_mat/rowSums(BTad_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_BTad_mat)){
  sorted_abs[[i]]<-sort(n_BTad_mat[i,], decreasing = T)
}
n_BTad_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                   ncol=length(sorted_abs[[1]]), byrow=TRUE)

# Plot the axis
plot(1e10,xlim = c(1,30),ylim = c(0,0.6),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15,20,25,30))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4,0.5,0.6),las = 0)

# Rank-abundance colors per treatment
a <- representative_RAD(norm_rad = n_BTad_mat,sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_BTad_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_BTad_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Treatment 1","Treatment 2", "Control"),
       col = line_cols, lwd = 3)


# Metamorphic
sample_classes <- c(rep(1, length(M1_s)),rep(2, length(M2_s)),
                    rep(3, length(M3_s)))
line_cols <- c("darkorange1","red3","green3")

BMet_mat<-rbind(BMet_T1, BMet_T2, BMet_Ctr)

# Normalized abundances
n_BMet_mat<-BMet_mat/rowSums(BMet_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_BMet_mat)){
  sorted_abs[[i]]<-sort(n_BMet_mat[i,], decreasing = T)
}
n_BMet_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                   ncol=length(sorted_abs[[1]]), byrow=TRUE)

# Plot the axis
plot(1e10,xlim = c(1,30),ylim = c(0,0.6),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15,20,25,30))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4,0.5,0.6),las = 0)

# Rank-abundance colors per treatment
a <- representative_RAD(norm_rad = n_BMet_mat,sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_BMet_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_BMet_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Treatment 1","Treatment 2", "Control"),
       col = line_cols, lwd = 3)

# Adults
sample_classes <- c(rep(1, length(A1_s)),rep(2, length(A2_s)),
                    rep(3, length(A3_s)))
line_cols <- c("darkorange1","red3","green3")

BAdl_mat<-rbind(BAdl_T1, BAdl_T2, BAdl_Ctr)

# Normalized abundances
n_BAdl_mat<-BAdl_mat/rowSums(BAdl_mat)

# Sort the abundances in a decreasing order
sorted_abs<-list()
for(i in 1:nrow(n_BAdl_mat)){
  sorted_abs[[i]]<-sort(n_BAdl_mat[i,], decreasing = T)
}
n_BAdl_mat<-matrix(unlist(sorted_abs), nrow=length(sorted_abs),
                   ncol=length(sorted_abs[[1]]), byrow=TRUE)

# Plot the axis
plot(1e10,xlim = c(1,30),ylim = c(0,0.6),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15,20,25,30))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4,0.5,0.6),las = 0)

# Rank-abundance colors per treatment
a <- representative_RAD(norm_rad = n_BAdl_mat,sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),border = NA)
a <- representative_RAD(norm_rad = n_BAdl_mat,sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),border = NA)
a <- representative_RAD(norm_rad = n_BAdl_mat,sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),border = NA)
legend("topright",bty = "n",
       legend = c("Treatment 1","Treatment 2", "Control"),
       col = line_cols, lwd = 3)



# MSD analysis

# Tadpole
d <- dist(x = n_BTad_mat,method = "manhattan")
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
d <- dist(x = n_BMet_mat,method = "manhattan")
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


# Networks' connectivity analysis

degree_df <- data.frame(Species = vertex.attributes(BTad_CtrNet)$name,
                        color = vertex.attributes(BTad_CtrNet)$color,
                        Phylum = vertex.attributes(BTad_CtrNet)$Taxon,
                        Ctr_degree = degree(BTad_CtrNet),
                        T1_degree = degree(BTad_T1Net),
                        T2_degree = degree(BTad_T2Net))

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
