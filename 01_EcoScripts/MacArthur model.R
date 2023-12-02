library(miaSim)
# Number of species ans resources
n_species <- 6
n_resources <- 3
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
  t_end = 5000,
  t_step = 1,
  t_external_events = 3000,
  t_external_durations = 50,
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
    t_end = 5000,
    t_step = 1,
  )
}

# Extract abundance tables
tse_models<-list()
for(i in 1:length(tse_list)){
  tse_models[[i]]<-as.matrix(assay(tse_list[[i]]))
}
# Plot standard deviation of the shannon diversity of the obtanied abundance tables
plot(models_sd(tse_models), type = "l", lwd = 3, col = "brown",
     xlab = "Time steps", ylab = "sd")


