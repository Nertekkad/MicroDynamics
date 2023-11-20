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

# Clean the data base
abs_genus<-abs_genus[,-which(colnames(abs_genus) == " g__")]
abs_genus<-abs_genus[,-which(is.na(colnames(abs_genus)))]
rownames(abs_genus) <- gsub("^X", "", rownames(abs_genus))

# Load meta data
metadata<-read_table("Metadata_Filtered.txt")
View(metadata)
metadata$Sample_Type

# Separate data
types <- unique(metadata$Sample_Type)

# BaAka humans
BaAka <- metadata$SampleID[which(metadata$Sample_Type == types[1]]
which(rownames(abs_genus) %in% BaAka)
T_BaAka <- abs_genus[which(rownames(abs_genus) %in% BaAka),]
# Bantu humans
Bantu <- metadata$SampleID[which(metadata$Sample_Type == "Bantu-Human"]
which(rownames(abs_genus) %in% Bantu)
T_Bantu <- abs_genus[which(rownames(abs_genus) %in% Bantu),]
# USA humans
BaAka <- metadata$SampleID[which(metadata$Sample_Type == "USA-Human"]
which(rownames(abs_genus) %in% BaAka)
T_BaAka <- abs_genus[which(rownames(abs_genus) %in% BaAka),]


abs_genus[grepl("^BaAka", rownames(abs_genus)), ]

metadata<-read_table("Metadata_Filtered.txt")
metadata$Sample_Type




