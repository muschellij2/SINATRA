#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

### Ensure the proper arguments are inputted
if (length(args) != 3) {
  stop("At least one argument must be supplied", call.=FALSE)
}

### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))

### Load in the R Libraries ###
library(truncnorm)
library(doParallel)
library(svd)
library(numbers)

### Load SINATRA functions ###
source("SINATRA_Code/ec_computation.R")
source("SINATRA_Code/generate_directions.R")
source("SINATRA_Code/gp_inference.R")
source("SINATRA_Code/metric_curve_simulation.R")
source("SINATRA_Code/plotting_functions.R")
source("SINATRA_Code/roc_curve_simulation.R")
source("SINATRA_Code/shape_reconstruction.R")
source("SINATRA_Code/simulated_data_generation.R")
source("SINATRA_Code/RATEv2.R")
sourceCpp("SINATRA_Code/BAKRGibbs.cpp")

### Set the parameters for the analysis ###
set.seed(4913, kind = "L'Ecuyer-CMRG")
n.simulations <- 100


# take this input from command line
num_causal_region <- args[1]
num_shared_region <- args[2]
causal_points <- args[3]
shared_points <- args[3]


######################################################################################
######################################################################################
######################################################################################

### Setup DoParallel ###
no_cores <- detectCores() - 1 
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl)  

### Run the analysis in Parallel ###
simulation_results <- foreach(i=1:n.simulations, .combine = 'rbind', .noexport = c('GaussKernel')) %:% 
  foreach(j=c(1,5,10,15,20), .combine = 'rbind', .noexport = c('GaussKernel')) %dopar% {
    
    set.seed(3*i+j)
    
    res <- tryCatch( generate_ROC_with_coned_directions(nsim = 50, curve_length = 50, grid_size = 25, distance_to_causal_point = 0.1, 
                                                        causal_points = causal_points,shared_points = shared_points, desired_num_cones = j, eta = 0.1, 
                                                        truncated = -1, two_curves = TRUE, ball = TRUE, ball_radius = 2.5, type = 'vertex',
                                                        min_points = 3,directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                                        mode = 'sphere', fpr = 0.05, start = 1, cusps = 50,
                                                        subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region),
                     error = function(x) {
                       return(matrix(nrow = 0,ncol = 3))
                     }
    )
    ### Label the results for each trial and directions ###
    rdf <- cbind(res, rep(j, dim(res)[1]) )
    rdf <- data.frame(cbind(rdf,rep(i,dim(rdf)[1])))
    rdf <- plyr::rename(rdf,c("X1" = "FPR","X2" = "TPR","X3" = "Class","X4" = "Index","X5" = "Num_Cones","X6" = "Trial"))
  }

stopCluster(cl)

######################################################################################
######################################################################################
######################################################################################
### Aggregate results ###
rdfmeans <- aggregate(simulation_results[c("FPR","TPR")],
                      by = list("Num_Cones" = simulation_results$Num_Cones,
                                "Index" = simulation_results$Index,
                                "Class" = simulation_results$Class), mean)

rdfmeans$Num_Cones <- as.factor(rdfmeans$Num_Cones)

### Plot results ###
#ROC_curve_plt <- ggplot(data <- rdfmeans[rdfmeans$Class == 1,],aes(x = FPR, y = TPR, color = Num_Cones)) +
#  geom_line(stat = "identity") +
#  labs(x = "FPR", y = "TPR") +
#  ggtitle(sprintf("causal:%d,shared:%d,type:%s",causal_points,shared_points,reconstruction_type)) +
#  geom_abline(intercept = 0, slope = 1)

######################################################################################
######################################################################################
######################################################################################
### save the results ###

# save the raw results
sim_results_file = sprintf("~/data/Results/SINATRA/Sphere_Simulation/data_ROC_causal%d_shared%d_%s.RData",
                           num_causal_region,num_shared_region,reconstruction_type)
save(simulation_results, file = sim_results_file)

# save the dataframe
df_results_file = sprintf("~/data/Results/SINATRA/Sphere_Simulation/df_ROC_causal%d_shared%d_%s.RData",
                          num_causal_region,num_shared_region,reconstruction_type)
save(rdfmeans, file = df_results_file)

