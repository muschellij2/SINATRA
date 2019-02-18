### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))

### Load in the R Libraries ###
library(truncnorm)
library(doParallel)
library(svd)

### Load SINATRA functions ###
source("Final Code for Functions (Really)/RATEv2.R")
sourceCpp("Final Code for Functions (Really)/sinatra/R/BAKRGibbs.cpp")

### Set the parameters for the analysis ###
set.seed(4913, kind = "L'Ecuyer-CMRG")
n.simulations <- 2
reconstruction_type <- "feature" #set how we assess each reconstruction
causal_points <- 15
shared_points <- 5


######################################################################################
######################################################################################
######################################################################################

### Setup DoParallel ###
no_cores <- detectCores() - 2 
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl)  

### Run the analysis in Parallel ###
simulation_results <- foreach(i=1:n.simulations, .combine = 'rbind', .noexport = c('GaussKernel')) %:% 
  foreach(j=c(1,5), .combine = 'rbind', .noexport = c('GaussKernel')) %dopar% {
    
    set.seed(3*i+j)
    
    res <- generate_ROC_with_coned_directions(nsim = 20, curve_length = 15, grid_size = 25, distance_to_causal_point = 0.1, 
                                              causal_points = causal_points,shared_points = shared_points, desired_num_cones = j, eta = 0.1, 
                                              truncated = FALSE, two_curves = TRUE, ball = TRUE, ball_radius = 2.5, type = reconstruction_type,
                                              min_points = 3,directions_per_cone = 4, cap_radius = 0.15, radius = 1)
    
    ### Label the results for each trial and directions ###
    rdf <- cbind(res, rep(j, dim(res)[1]) )
    rdf <- data.frame(cbind(rdf,rep(i,dim(rdf)[1])))
    rdf <- plyr::rename(rdf,c("X1" = "FPR","X2" = "TPR","X3" = "Index","X4" = "Num_Directions","X5" = "Trial"))
  }

stopCluster(cl)

######################################################################################
######################################################################################
######################################################################################
### Aggregate results ###
rdfmeans <- aggregate(simulation_results[c("FPR","TPR")],
                      by = list("Num_Directions" = simulation_results$Num_Directions,
                                "Index" = simulation_results$Index), mean)

rdfmeans$Num_Directions <- as.factor(rdfmeans$Num_Directions)

### Plot results ###
ROC_curve_plt <- ggplot(data <- rdfmeans,aes(x = FPR, y = TPR, color = Num_Directions)) +
  geom_line(stat = "identity") +
  labs(x = "FPR", y = "TPR") +
  ggtitle(sprintf("causal:%d,shared:%d,type:%s",causal_points,shared_points,reconstruction_type)) +
  geom_abline(intercept = 0, slope = 1)

print(ROC_curve_plt)
######################################################################################
######################################################################################
######################################################################################
### save the results ###

# save the raw results
sim_results_file = "~/Desktop/Simulations/Results/data_metric_curve_causal15_shared5_vertex.RData"
save(simulation_results, file = sim_results_file)

# save the final plot
avg_curve_file = "~/Desktop/Simulations/Results/plot_metric_curve_causal15_shared5_vertex.pdf"
pdf(file = avg_curve_file)
plot(avg_curve_plt)
dev.off()