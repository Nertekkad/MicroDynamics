library(seqtime)
library(igraph)
library(phyloseq)
library(mlBioNets)
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
T_Collapsed<-T_collapse(T_table = T_table, O_table = O_table,
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
abs_mat()

lay <- layoutMultiplex(g.list, layout="fr", ggplot.format=F, box=T)
plot_multiplex3D(g.list, layer.layout=lay,
                 layer.colors=rainbow(length(g.list)),
                 layer.shift.x=0.5, layer.space=2,
                 layer.labels=c("Before travel", "During travel", "After travel"),
                 layer.labels.cex=1.5, node.size.values=node.ab.matrix,
                 node.size.scale=0.6,
                 node.colors=vertex.attributes(T1_aracne)$color,
                 edge.colors="white",
                 node.colors.aggr=vertex.attributes(T1_aracne)$color,
                 show.aggregate=T)