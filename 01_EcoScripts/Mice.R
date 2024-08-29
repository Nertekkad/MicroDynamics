#### Data procesing ####

setwd("~/Mice")
library(readr)
# ASV abundance table
asv_table<-as.data.frame(read_tsv("counts.tsv"))
ASVs<-as.vector(asv_table[,1])
asv_table<-asv_table[,-1]
rownames(asv_table)<-ASVs

# ASV metadata table
asv_meta<-as.data.frame(read_tsv("metadata.tsv"))

# Perturbations' table
asv_pert<-as.data.frame(read_tsv("perturbations.tsv"))

# Taxonomic table
asv_taxa<-as.data.frame(read_tsv("rdp_species.tsv"))

# Table collapse at genus level
library(mlBioNets)
asv_table2<-T_collapse(is_phyloseq = F, T_table = asv_taxa, O_table = asv_table,
                       names_level = "Genus")

# Subject 2
S2<-asv_meta[which(asv_meta$subject == 2),]
S2<-S2[order(S2$time), ]

# Basal

a<-S2$time[1]
b<-asv_pert[which(asv_pert$subject == 2),]$start[1]
b<-S2$time[which(S2$time==b)-1]

basal_S2<-S2[which(S2$time == a):which(S2$time == b),]$sampleID
basal_S2<-asv_table2[basal_S2,]

# Fat-diet

a<-asv_pert[which(asv_pert$subject == 2),]$start[1]
b<-asv_pert[which(asv_pert$subject == 2),]$end[1]

fat_S2<-S2[which(S2$time == a):which(S2$time == b),]$sampleID
fat_S2<-asv_table2[fat_S2,]

# Recovered 1

a<-asv_pert[which(asv_pert$subject == 2),]$end[1]
a<-S2$time[which(S2$time==a)+1]
b<-asv_pert[which(asv_pert$subject == 2),]$start[2]
b<-S2$time[which(S2$time==b)-1]

R1_S2<-S2[which(S2$time == a):which(S2$time == b),]$sampleID
R1_S2<-asv_table2[R1_S2,]

# Vancomycin

a<-asv_pert[which(asv_pert$subject == 2),]$start[2]
b<-asv_pert[which(asv_pert$subject == 2),]$end[2]

van_S2<-S2[which(S2$time == a):which(S2$time == b),]$sampleID
van_S2<-asv_table2[van_S2,]

# Recovered 2

a<-asv_pert[which(asv_pert$subject == 2),]$end[2]
a<-S2$time[which(S2$time==a)+1]
b<-asv_pert[which(asv_pert$subject == 2),]$start[3]
b<-S2$time[which(S2$time==b)-1]

R2_S2<-S2[which(S2$time == a):which(S2$time == b),]$sampleID
R2_S2<-asv_table2[R2_S2,]

# Gentamicin

a<-asv_pert[which(asv_pert$subject == 2),]$start[3]
b<-asv_pert[which(asv_pert$subject == 2),]$end[3]

gen_S2<-S2[which(S2$time == a):which(S2$time == b),]$sampleID
gen_S2<-asv_table2[gen_S2,]

# Recovered 3

a<-asv_pert[which(asv_pert$subject == 2),]$end[3]
a<-S2$time[which(S2$time==a)+1]
b<-S2$time[which.max(S2$time)]

R3_S2<-S2[which(S2$time == a):which(S2$time == b),]$sampleID
R3_S2<-asv_table2[R3_S2,]

#### Topological network analysis ####

ml_S2<-list(basal_S2, fat_S2, R1_S2, van_S2, R2_S2, gen_S2, R3_S2)

Nbasal_S2<-net_inference(basal_S2, "aracne")
Nfat_S2<-net_inference(fat_S2, "aracne")
NR1_S2<-net_inference(R1_S2, "aracne")
Nvan_S2<-net_inference(van_S2, "aracne")
NR2_S2<-net_inference(R2_S2, "aracne")
Ngen_S2<-net_inference(gen_S2, "aracne")
NR3_S2<-net_inference(R3_S2, "aracne")

fat_ml<-list(Nbasal_S2, Nfat_S2, NR1_S2)
van_ml<-list(NR1_S2, Nvan_S2, NR2_S2)
gen_ml<-list(NR2_S2, Ngen_S2, NR3_S2)

treatments <- c(1, 2, 3)
mice_ml_properties<-rbind(ml_properties(fat_ml, treatments),
                          ml_properties(van_ml, treatments),
                          ml_properties(gen_ml, treatments))

Stage = c(rep("Fat diet", length(fat_ml)),
          rep("Vancomycin",length(van_ml)),
          rep("Gentamycin", length(gen_ml)))

mice_ml_properties <- cbind(mice_ml_properties, Stage)

require(ggplot2)

p1<-ggplot(mice_ml_properties, aes(x = Treatments, y = Mean_degree, color = Stage)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  labs(title = "",
       x = "Period",
       y = "Mean degree") +
  scale_color_manual(values = c("darkgreen", "purple4", "brown4")) +
  theme_minimal()

p2<-ggplot(mice_ml_properties, aes(x = Treatments, y = sd_degree, color = Stage)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  labs(title = "",
       x = "Period",
       y = "SD of degree") +
  scale_color_manual(values = c("darkgreen", "purple4", "brown4")) +
  theme_minimal()

p3<-ggplot(mice_ml_properties, aes(x = Treatments, y = Clusterization, color = Stage)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  labs(title = "",
       x = "Period",
       y = "Transitivity") +
  scale_color_manual(values = c("darkgreen", "purple4", "brown4")) +
  theme_minimal()

p4<-ggplot(mice_ml_properties, aes(x = Treatments, y = Edge_density, color = Stage)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  labs(title = "",
       x = "Period",
       y = "Edge density") +
  scale_color_manual(values = c("darkgreen", "purple4", "brown4")) +
  theme_minimal()

p5<-ggplot(mice_ml_properties, aes(x = Treatments, y = Connected_nodes, color = Stage)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  labs(title = "",
       x = "Period",
       y = "Proportion of linked nodes") +
  scale_color_manual(values = c("darkgreen", "purple4", "brown4")) +
  theme_minimal()

p6<-ggplot(mice_ml_properties, aes(x = Treatments, y = Modularity, color = Stage)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +  
  labs(title = "",
       x = "Period",
       y = "Modularity") +
  scale_color_manual(values = c("darkgreen", "purple4", "brown4")) +
  theme_minimal()

library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol = 3)



##### Phyloseq object #####

library(phyloseq)

asv_taxa2<-asv_taxa[,-c(1,2)]
rownames(asv_taxa2) <- ASVs
asv_taxa2
asv_table2<-otu_table(asv_table, taxa_are_rows = T)

sample_IDs<-c(rownames(basal_S2), rownames(fat_S2), rownames(R1_S2), rownames(van_S2),
              rownames(R2_S2), rownames(gen_S2), rownames(R3_S2))
treatments<-c(rep("Basal", length(rownames(basal_S2))),
              rep("Fat diet", length(rownames(fat_S2))),
              rep("Recovered 1", length(rownames(R1_S2))),
              rep("Vancomycin", length(rownames(van_S2))),
              rep("Recovered 2", length(rownames(R2_S2))),
              rep("Gentamycin", length(rownames(gen_S2))),
              rep("Recovered 3", length(rownames(R3_S2))))
meta_samples<-data.frame(
  "Samples"=sample_IDs,
  "Treatments"=treatments
)
row.names(meta_samples)<-sample_IDs

asv_table2<-asv_table2[,which(colnames(asv_table2) %in% sample_IDs)]

meta_samples2<-sample_data(meta_samples)

physeq_mice<-merge_phyloseq(asv_table2, asv_taxa2)
physeq_mice2<-phyloseq(physeq_mice, meta_samples2)

#### Richness comparison ####

library(MicrobiotaProcess)
library(UpSetR)

upsetda2 <- get_upset(obj=physeq_mice2, factorNames="Treatments")
upset(upsetda2, sets=unique(as.vector(sample_data(physeq_mice2)$Treatments)), 
      sets.bar.color = "#56B4E9",
      order.by = "freq", 
      empty.intersections = "on")

#### Dysbiosis ####

library(dysbiosisR)
library(ggplot2)
library(microbiome)
library(dplyr)

# Bray-Curtis distance matrix
dist.mat <- phyloseq::distance(physeq_mice2, "bray")
# Get reference samples
ref.samples <- sample_names(subset_samples(physeq_mice2, 
                                           Treatments == "Basal"))


# Community level variation analysis
dysbiosis_1 <- dysbiosisMedianCLV(physeq_mice2,
                                  dist_mat = dist.mat,
                                  reference_samples = ref.samples)
# We sample the data set identifying as dysbiotic the data under the 90th percentile
dysbiosis_thres <- quantile(subset(dysbiosis_1, Treatments == "Fat diet")$score, 0.9)
normobiosis_thres <- quantile(subset(dysbiosis_1, Treatments == "Fat diet")$score, 0.1)

dysbiosis_1 <- dysbiosis_1 |> 
  mutate(isDysbiostic = ifelse(score >= dysbiosis_thres, TRUE, FALSE))

# Dysbiosis plot measures according to CLV method
p1 <- plotDysbiosis(df=dysbiosis_1,
                    xvar="Treatments",
                    yvar="score",
                    show_points = FALSE) +
  labs(x="", y="Dysbiosis Score") +
  theme_bw(base_size = 14)
p1


# Dysbiosis plot measures according to euclidean method
dysbiosis_2 <- euclideanDistCentroids(physeq_mice2,
                                      dist_mat = dist.mat,
                                      use_squared = TRUE,
                                      group_col = "Treatments",
                                      control_label = "Basal",
                                      case_label = "Fat diet")

p2 <- plotDysbiosis(df=dysbiosis_2,
                    xvar="Treatments",
                    yvar="CentroidDist_score",
                    show_points = FALSE) +
  labs(x="", y="Dysbiosis Score (Centroid)") +
  theme_bw(base_size = 14)
p2


# Dysbiosis plot measures according to the Combined alpha-beta diversity based score
dysbiosis_3 <- combinedShannonJSD(physeq_mice2,
                                  reference_samples = ref.samples)


p3 <- plotDysbiosis(df=dysbiosis_3,
                    xvar="Treatments",
                    yvar="ShannonJSDScore",
                    show_points = FALSE)+
  labs(x="", y="Shannon-JSD\nDysbiosis Score") +
  theme_bw(base_size = 14)
p3


# Test for outliers’ detection that accounts for the wide range of
# microbiome phenotypes observed in a typical set of healthy individuals
# and for intra-individual temporal variation.
cloud.results <- cloudStatistic(physeq_mice2,
                                dist_mat = dist.mat,
                                reference_samples = ref.samples,
                                ndim=-1,
                                k_num=5)

p4 <- plotDysbiosis(df=cloud.results,
                    xvar="Treatments",
                    yvar="log2Stats",
                    show_points = FALSE) +
  labs(x="", y="Dysbiosis CLOUD Score") +
  theme_bw(base_size = 14)
p4


#### Early warnings metrics ####

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

# Diversity function for abundance tables
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
          S <- length(ab_table[, i])
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

library(earlywarnings)
library(EWS)
library(EWSmethods)
library(codyn)
library(vegan)

mice_abundance<-rbind(basal_S2, fat_S2, R1_S2,
                      van_S2, R2_S2, gen_S2, R3_S2)

mice_data <- data.frame(time = seq(1:length(ab_table_div(t(mice_abundance), "shannon"))),
                       abundance = ab_table_div(t(mice_abundance), "shannon"))

ews_metrics <- c("SD","ar1","skew")

roll_ews <- uniEWS(data = mice_data, metrics =  ews_metrics, method = "rolling", winsize = 50)

plot(roll_ews,  y_lab = "Abundances")

exp_ews <- uniEWS(data = mice_data, metrics =  ews_metrics, method = "expanding",
                  burn_in = 10, threshold = 2,  tail.direction = "one.tailed")
plot(exp_ews, y_lab = "Abundances")
