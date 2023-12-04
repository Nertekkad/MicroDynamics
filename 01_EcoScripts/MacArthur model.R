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
  stochastic = T
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

# Function for estimate the shannon, pielou and simpson diversity
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
for(i in 1:5){
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
CRsims <- list()
for(i in 1:5){
  CRsims[[i]] <- simulateConsumerResource(
    n_species = n_species,
    n_resources = n_resources, names_species = letters[seq_len(n_species)],
    names_resources = paste0("res", LETTERS[seq_len(n_resources)]), E = matE,
    x0 = rep(0.001, n_species), resources = resources_c,
    growth_rates = runif(n_species),
    error_variance = 0.001,
    t_end = 10000,
    t_step = 1,
    t_external_events = 5000,
    t_external_durations = 50,
    sigma_external = 0.05,
    stochastic = T
  )
}

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

# Tables normalization
norm_models<-function(ab_tables){
  normalized_tables<-list()
  for(i in 1:length(ab_tables)){
    normalized_tables[[i]]<-ab_tables[[i]]/colSums(ab_tables[[i]])
  }
  return(normalized_tables)
}

CRsims_T<-norm_models(CRsims_T)

# Diversity over time
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
plot(1e10,xlim = c(1,20),ylim = c(0,1),
     xlab = "Species rank",ylab = "Abundance",cex.lab = 1.5,axes = FALSE)
sfsmisc::eaxis(side = 1,at = c(1,10))
sfsmisc::eaxis(side = 2,at = c(0,2,4,6,8,10,12,14,16,18,20),las = 0)

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
