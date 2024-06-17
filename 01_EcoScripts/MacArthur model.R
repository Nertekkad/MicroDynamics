library(miaSim)
# Number of species ans resources
n_species <- 10
n_resources <- 5
# Random efficiency matrix for consumer resource model
matE <- randomE(
  n_species = n_species, n_resources = n_resources,
  mean_consumption = 2, mean_production = 2, maintenance = 0.4
)
# Initial concentrations of resources
resources_c <- rep(100, n_resources)
# MacArthur consumer-resource model with a perturbation
tse <- simulateConsumerResource(
  n_species = n_species,
  n_resources = n_resources, names_species = letters[seq_len(n_species)],
  names_resources = paste0("res", LETTERS[seq_len(n_resources)]), E = matE,
  x0 = rep(0.001, n_species), resources = resources_c,
  growth_rates = runif(n_species),
  error_variance = 0.001,
  t_end = 10000,
  t_step = 1,
  t_external_events = 5000,
  t_external_durations = 1000,
  sigma_external = 0.05,
  stochastic = T,
  norm = T
)
# Extract abundance and resource tables
abundances<-as.data.frame(assay(tse))
resources<-as.data.frame(metadata(tse)$resources)
# Time units
times<-as.vector(resources$time)
resources<-resources[,-ncol(resources)]
# Plot populations
matplot(x = times, y = t(abundances), type = "l",
        xlab = "Time", ylab = "Population",
        main = "McArthur consumer-resource model", lwd = 2, lty = 1)
# Plot resources
matplot(x = times, y = resources, type = "l",
        xlab = "Time", ylab = "Resources",
        main = "McArthur consumer-resource model", lwd = 2, lty = 1)


# Before simulate the system, we shall find at which time step it gets stability

# Function for estimate the Shannon, Pielou and Simpson diversity
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

# Function for estimate the standard deviation of the diversity of n-iterations
# of the model
models_sd<-function(ab_tables){
  div_list<-list()
  for(i in 1:length(ab_tables)){
    div_list[[i]]<-ab_table_div(ab_tables[[i]], "shannon")
  }
  div_mat<-matrix(unlist(div_list), ncol = length(ab_tables))
  mat_sd<-c()
  for(i in 1:nrow(div_mat)){
    mat_sd[i]<-sd(div_mat[i,])
  }
  return(mat_sd)
}
# Let's iterate the model 
tse_list<-list()
for(i in 1:3){
  MatE <- randomE(
    n_species = n_species, n_resources = n_resources,
    mean_consumption = 2, mean_production = 2, maintenance = 0.4
  )
  resources_c <- rep(100, n_resources)
  tse_list[[i]] <- simulateConsumerResource(
    n_species = n_species,
    n_resources = n_resources, names_species = letters[seq_len(n_species)],
    names_resources = paste0("res", LETTERS[seq_len(n_resources)]), E = MatE,
    x0 = rep(0.001, n_species), resources = resources_c,
    growth_rates = runif(n_species),
    t_end = 10000,
    t_step = 1,
    norm = T
  )
}

# Extract abundance tables
tse_models<-list()
for(i in 1:length(tse_list)){
  tse_models[[i]]<-as.matrix(assay(tse_list[[i]]))
}
# Plot standard deviation of the Shannon diversity of the obtanied abundance tables
plot(models_sd(tse_models), type = "l", lwd = 3, col = "brown",
     xlab = "Time steps", ylab = "sd")


# n-iterated model
#CRsims <- list()
#for(i in 1:3){
#  MatE <- randomE(
#    n_species = n_species, n_resources = n_resources,
#    mean_consumption = 2, mean_production = 2, maintenance = 0.4
#  )
#  resources_c <- rep(100, n_resources)
#  CRsims[[i]] <- simulateConsumerResource(
#    n_species = n_species,
#    n_resources = n_resources, names_species = letters[seq_len(n_species)],
#    names_resources = paste0("res", LETTERS[seq_len(n_resources)]), E = matE,
#    x0 = rep(0.001, n_species), resources = resources_c,
#    growth_rates = runif(n_species),
#    error_variance = 0.001,
#    t_end = 10000,
#    t_step = 1,
#    t_external_events = 5000,
#    t_external_durations = 1000,
#    sigma_external = 0.05,
#    stochastic = T,
#    norm = T
#  )
#}

load("~/MicroDynamics/02_Data/crm.RData")

# Function for estimate the diversity variance over time
ab_tables_div<-function(ab_tables_list, diversity_type){
  if(diversity_type == "shannon"){
    div_list<-list()
    for(j in 1:length(ab_tables_list)){
      div_table<-c()
      for(i in 1:ncol(ab_tables_list[[j]])){
        div_table[i]<-Div_Shannon(ab_tables_list[[j]][, i])
      }
      div_list[[j]]<-div_table
    }
  } else
    if(diversity_type == "simpson"){
      div_list<-list()
      for(j in 1:length(ab_tables_list)){
        div_table<-c()
        for(i in 1:ncol(ab_tables_list[[j]])){
          div_table[i]<-Dom_Simpson(ab_tables_list[[j]][, i])
        }
        div_list[[j]]<-div_table
      }
    } else
      if(diversity_type == "pielou"){
        div_list<-list()
        for(j in 1:length(ab_tables_list)){
          div_table<-c()
          for(i in 1:ncol(ab_tables_list[[j]])){
            div_table[i]<-Eq_Pielou(ab_tables_list[[j]][, i])
          }
          div_list[[j]]<-div_table
        }
      } else
        if(diversity_type == "ginisimpson"){
          div_list<-list()
          for(j in 1:length(ab_tables_list)){
            div_table<-c()
            for(i in 1:ncol(ab_tables_list[[j]])){
              div_table[i]<-1-Dom_Simpson(ab_tables_list[[j]][, i])
            }
            div_list[[j]]<-div_table
          }
        }
  return(div_list)
}

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

# Extract abundance tables
CRsims_T<-list()
for(i in 1:length(CRsims)){
  CRsims_T[[i]]<-as.matrix(assay(CRsims[[i]]))
}

# Diversity over time
library(ggplot2)
pielou_df <- data.frame(
  "Time" = 1:length(ab_tables_div(CRsims_T, "pielou")[[1]]), 
  "Diversity" = ab_tables_div(CRsims_T, "pielou")[[1]]
)
ggplot(pielou_df, aes(x = Time, y = Diversity)) + 
  geom_point() + theme_classic()

# Separate the data
sep_tables<-function(ab_tables, start_p, end_p){
  sep_tables<-list()
  for(i in 1:length(ab_tables)){
    sep_tables[[i]]<-ab_tables[[i]][,start_p:end_p]
  }
  return(sep_tables)
}

basal_tables <- sep_tables(CRsims_T, 400, 499)
pert_tables <- sep_tables(CRsims_T, 500, 599)
post_tables <- sep_tables(CRsims_T, 600, 699)
recovered_tables <- sep_tables(CRsims_T, 700, 1000)

# Pielou diversity
pielou_basal<-unlist(ab_tables_div(basal_tables, "pielou"))
pielou_pert<-unlist(ab_tables_div(pert_tables, "pielou"))
pielou_post<-unlist(ab_tables_div(post_tables, "pielou"))
pielou_recovered<-unlist(ab_tables_div(recovered_tables, "pielou"))

# Simpson dominance
simpson_basal<-unlist(ab_tables_div(basal_tables, "simpson"))
simpson_pert<-unlist(ab_tables_div(pert_tables, "simpson"))
simpson_post<-unlist(ab_tables_div(post_tables, "simpson"))
simpson_recovered<-unlist(ab_tables_div(recovered_tables, "simpson"))

df_diversity<-data.frame(
  Period = c(rep(1, length(pielou_basal)), rep(2, length(pielou_pert)),
             rep(3, length(pielou_post)), rep(4, length(pielou_recovered))),
  Pielou =  c(pielou_basal, pielou_pert, pielou_post, pielou_recovered),
  Simpson = c(simpson_basal, simpson_pert, simpson_post, simpson_recovered)
)

# Evenness violin plot
library(ggstatsplot)
ggbetweenstats(
  data  = df_diversity,
  x     = Period,
  y     = Pielou,
  title = "Pielou evenness"
)
# Dominance violin plot
ggbetweenstats(
  data  = df_diversity,
  x     = Period,
  y     = Simpson,
  title = "Simpson dominance"
)

# Rank abundance analysis
rank_abs<-function(ab_tables, start_p, end_p){
  rank_abs<-list()
  for(j in 1:length(ab_tables)){
    sum_abs<-c()
    for(i in 1:nrow(ab_tables[[1]])){
      sum_abs[i]<-sum(ab_tables[[j]][,start_p:end_p][i,])
    }
    rank_abs[[j]]<-sort(sum_abs/nrow(ab_tables[[1]]), decreasing = T)
  }
  rank_abs_mat<-matrix(unlist(rank_abs), nrow=length(ab_tables),
                       ncol=nrow(ab_tables[[1]]), byrow=TRUE)
  return(rank_abs_mat)
}

basal_rank <- rank_abs(CRsims_T, 400, 499)
pert_rank <- rank_abs(CRsims_T, 500, 599)
post_rank <- rank_abs(CRsims_T, 600, 699)
recovered_rank <- rank_abs(CRsims_T, 700, 1000)
full_rank<-rbind(basal_rank, pert_rank, post_rank, recovered_rank)
full_rank<-full_rank/full_rank[,1][which.max(full_rank[,1])]

# RAD-analysis
library(RADanalysis)

sample_classes <- c(rep(1, length(CRsims_T)),rep(2, length(CRsims_T)),
                    rep(3, length(CRsims_T)), rep(4, length(CRsims_T)))
line_cols <- c("green3","red3","darkorange1", "dodgerblue4")
# Plot the axis
plot(1e10,xlim = c(1,50),ylim = c(0,0.5),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,5,10,15,20,25,30,35,40,45,50))
sfsmisc::eaxis(side = 2,at = c(0,0.1,0.2,0.3,0.4,0.5),las = 0)

# Plot the curves
a <- representative_RAD(norm_rad = full_rank,
                        sample_ids = which(sample_classes == 1),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[1],0.5),
                        border = NA)
a <- representative_RAD(norm_rad = full_rank,
                        sample_ids = which(sample_classes == 2),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[2],0.5),
                        border = NA)
a <- representative_RAD(norm_rad = full_rank,
                        sample_ids = which(sample_classes == 3),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[3],0.5),
                        border = NA)
a <- representative_RAD(norm_rad = full_rank,
                        sample_ids = which(sample_classes == 4),
                        plot = TRUE,confidence = 0.9,with_conf = TRUE,
                        col = scales::alpha(line_cols[4],0.5),
                        border = NA)
legend("topright",bty = "n",
       legend = c("Basal","Perturbated","Post-perturbated",
                  "Recovered"), col = line_cols,lwd = 3)

# Distance matrix using Manhattan distance
d <- dist(x = full_rank,method = "manhattan")
# Ordination using classical multi-dimensional scaling
mds <- cmdscale(d = d,k = 5,eig = TRUE)
# Points' plot
plot(mds$points,xlab = "First coordinate",ylab = "Second coordinate",
     pch = 19, cex =1,col = line_cols[sample_classes],
     main = "Multi-Dimensional Scaling plot")
# Representative points with error bars
a <- representative_point(input = mds$points,
                          ids = which(sample_classes == 1),
                          col = scales::alpha(line_cols[1],0.5),
                          plot = TRUE,standard_error_mean = TRUE,
                          pch = 19, cex = 4)
a <- representative_point(input = mds$points,
                          ids = which(sample_classes == 2),
                          col = scales::alpha(line_cols[2],0.5),
                          plot = TRUE,standard_error_mean = TRUE,
                          pch = 19, cex = 4)
a <- representative_point(input = mds$points,
                          ids = which(sample_classes == 3),
                          col = scales::alpha(line_cols[3],0.5),
                          plot = TRUE,standard_error_mean = TRUE,
                          pch = 19, cex = 4)
a <- representative_point(input = mds$points,
                          ids = which(sample_classes == 4),
                          col = scales::alpha(line_cols[4],0.5),
                          plot = TRUE,standard_error_mean = TRUE,
                          pch = 19, cex = 4)
legend("topleft",bty = "n",
       legend = c("Basal","Perturbated", "Post-perturbated",
                  "Recovered"), col = line_cols, pch = 19)



########## p1$layers

otumat<-CRsims_T[[1]][,-c(1:399)]
otumat<-as.matrix(otumat)
rownames(otumat)<-paste0("OTU", 1:nrow(otumat))
colnames(otumat)<-paste0("Sample", 1:ncol(otumat))

library(phyloseq)
otumat <- otu_table(otumat, taxa_are_rows = T)

sampledata <- sample_data(data.frame(
  Samples = c(rep("basal", 100), rep("pert", 100),
              rep("post", 100), rep("rec", 301)),
  Ab_mean = colSums(otu_table(otumat))/length(otu_table(otumat)),
  row.names = colnames(otu_table(otumat)),
  stringsAsFactors = F
))

physeq <- merge_phyloseq(otumat, sampledata)




# Dysbiosis analysis

library(dysbiosisR)
# Bray-Curtis distance matrix
dist.mat <- phyloseq::distance(physeq, "bray")
# Get reference samples
ref.samples <- sample_names(subset_samples(physeq, 
                                           Samples == "basal"))
# Community level variation analysis
dysbiosis_1 <- dysbiosisMedianCLV(physeq,
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
dysbiosis_2 <- euclideanDistCentroids(physeq,
                                      dist_mat = dist.mat,
                                      use_squared = TRUE,
                                      group_col = "Samples",
                                      control_label = "basal",
                                      case_label = "pert")

p2 <- plotDysbiosis(df=dysbiosis_2,
                    xvar="Samples",
                    yvar="CentroidDist_score",
                    colors=c(basal="green4", pert="red3",
                             post="orange3", rec="blue4"),
                    show_points = FALSE) +
  labs(x="", y="Dysbiosis Score (Centroid)") +
  theme_bw(base_size = 14)
p2

# Dysbiosis plot measures according to the Combined alpha-beta diversity based score
dysbiosis_3 <- combinedShannonJSD(physeq,
                                  reference_samples = ref.samples)


p3 <- plotDysbiosis(df=dysbiosis_3,
                    xvar="Samples",
                    yvar="ShannonJSDScore",
                    colors=c(basal="green4", pert="red3",
                             post="orange3", rec="blue4"),
                    show_points = FALSE)+
  labs(x="", y="Shannon-JSD\nDysbiosis Score") +
  theme_bw(base_size = 14)
p3



library(earlywarnings)
# Potential plot applied to the analysis of diversity over time
x <- ab_tables_div(CRsims_T, "pielou")

xmat <- matrix(unlist(x), length(x), length(x[[1]]))
xmean <- c()
for(i in 1:ncol(xmat)){
  xmean[i] <- mean(xmat[,i])
}; xmean

res <- movpotential_ews(xmean, param = NULL)
p <- PlotPotential(res$res, title = '', 
                   xlab.text = '', ylab.text = '', 
                   cutoff = 0.5, plot.contours = TRUE, binwidth = 0.2)
print(p)


# Resources analysis
CRsims_R <- list()
for(i in 1:length(CRsims)){
  CRsims_R[[i]] <- t(as.data.frame(metadata(CRsims[[i]])$resources))
}

# Separate tables
basal_tablesR <- sep_tables(CRsims_R, 400, 499)
pert_tablesR <- sep_tables(CRsims_R, 500, 599)
post_tablesR <- sep_tables(CRsims_R, 600, 699)
recovered_tablesR <- sep_tables(CRsims_R, 700, 1000)

# Pielou diversity
pielou_basalR<-unlist(ab_tables_div(basal_tablesR, "pielou"))
pielou_pertR<-unlist(ab_tables_div(pert_tablesR, "pielou"))
pielou_postR<-unlist(ab_tables_div(post_tablesR, "pielou"))
pielou_recoveredR<-unlist(ab_tables_div(recovered_tablesR, "pielou"))

df_diversityR<-data.frame(
  Period = c(rep(1, length(pielou_basalR)), rep(2, length(pielou_pertR)),
             rep(3, length(pielou_postR)), rep(4, length(pielou_recoveredR))),
  Pielou =  c(pielou_basalR, pielou_pertR, pielou_postR, pielou_recoveredR)
)

# Evenness violin plot
library(ggstatsplot)
ggbetweenstats(
  data  = df_diversityR,
  x     = Period,
  y     = Pielou,
  title = "Pielou evenness of resources"
)


library(earlywarnings)
library(EWS)
library(EWSmethods)
library(codyn)
library(vegan)

CR_data <- data.frame(time = seq(1:length(ab_table_div(CRsims_T[[1]], "shannon"))),
                        abundance = ab_table_div(CRsims_T[[1]], "shannon"))

ews_metrics <- c("SD","ar1","skew")

roll_ews <- uniEWS(data = CR_data, metrics =  ews_metrics, method = "rolling", winsize = 50)

plot(roll_ews,  y_lab = "Abundances")

exp_ews <- uniEWS(data = CR_data, metrics =  ews_metrics, method = "expanding",
                  burn_in = 10, threshold = 2,  tail.direction = "one.tailed")

plot(exp_ews, y_lab = "Abundances")
