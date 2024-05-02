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
asv_table2<-t(asv_table2)

# Subject 2
S2<-asv_meta[which(asv_meta$subject == 2),]
S2<-S2[order(S2$time), ]

# Basal

a<-S2$time[1]
b<-asv_pert[which(asv_pert$subject == 2),]$start[1]
b<-S2$time[which(S2$time==b)-1]

basal_S2<-S2[which(S2$time == a):which(S2$time == b),]$sampleID
basal_S2<-asv_table2[, basal_S2]

# Fat-diet

a<-asv_pert[which(asv_pert$subject == 2),]$start[1]
b<-asv_pert[which(asv_pert$subject == 2),]$end[1]

fat_S2<-S2[which(S2$time == a):which(S2$time == b),]$sampleID
fat_S2<-asv_table2[, fat_S2]

# Recovered 1

a<-asv_pert[which(asv_pert$subject == 2),]$end[1]
a<-S2$time[which(S2$time==a)+1]
b<-asv_pert[which(asv_pert$subject == 2),]$start[2]
b<-S2$time[which(S2$time==b)-1]

R1_S2<-S2[which(S2$time == a):which(S2$time == b),]$sampleID
R1_S2<-asv_table2[, R1_S2]

# Vancomycin

a<-asv_pert[which(asv_pert$subject == 2),]$start[2]
b<-asv_pert[which(asv_pert$subject == 2),]$end[2]

van_S2<-S2[which(S2$time == a):which(S2$time == b),]$sampleID
van_S2<-asv_table2[, van_S2]

# Recovered 2

a<-asv_pert[which(asv_pert$subject == 2),]$end[2]
a<-S2$time[which(S2$time==a)+1]
b<-asv_pert[which(asv_pert$subject == 2),]$start[3]
b<-S2$time[which(S2$time==b)-1]

R2_S2<-S2[which(S2$time == a):which(S2$time == b),]$sampleID
R2_S2<-asv_table2[, R2_S2]

# Gentamicin

a<-asv_pert[which(asv_pert$subject == 2),]$start[3]
b<-asv_pert[which(asv_pert$subject == 2),]$end[3]

gen_S2<-S2[which(S2$time == a):which(S2$time == b),]$sampleID
gen_S2<-asv_table2[, gen_S2]

# Recovered 3

a<-asv_pert[which(asv_pert$subject == 2),]$end[3]
a<-S2$time[which(S2$time==a)+1]
b<-S2$time[which.max(S2$time)]

R3_S2<-S2[which(S2$time == a):which(S2$time == b),]$sampleID
R3_S2<-asv_table2[, R3_S2]


