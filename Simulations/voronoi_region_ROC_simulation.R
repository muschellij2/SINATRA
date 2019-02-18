#################################################################
#################################################################
#################################################################
### Simulation Purpose ###
# Analyze how power is affected by number of initial cones and directions per cone.
#

#################################################################
#################################################################
#################################################################

library(R.utils)
library(purrr)
library(functional)
library(dichromat)
library(ggplot2)
library(Rcpp)
library(foreach)

sourceDirectory("~/projects/Research/SINATRA/oscar_simulations/SINATRA_Code/")
sourceCpp("~/projects/Research/SINATRA/oscar_simulations/SINATRA_Code/BAKRGibbs.cpp")
source("~/projects/Research/SINATRA/oscar_simulations/R_Scripts/Package_Setup.R")

#################################################################
#################################################################
#################################################################
### Function to generate the ROC below ###


GenerateVoronoiRegionROC <- function(initial_cones, directions_per_cone){
  
  set.seed(288)
  
  ### Parameters to analyze ###
  initial_cones <- initial_cones
  directions_per_cone <- directions_per_cone
  curve_length <- 100
  desired_num_cones <- 1
  sphere_smoothness <- 3
  ball_radius <- 2
  num_cusps <- 150
  class <- 1
  
  
  ### Generate durections ###
  directions <- generate_equidistributed_cones(initial_cones,0.1,directions_per_cone)
  current_num_cones = dim(directions)[1]/(directions_per_cone)
  
  ### Generate Data ###
  sphere_data <- generate_data_sphere_simulation(nsim = 50, curve_length = curve_length,directions, noise_points = 5, causal_points = 5,
                                                 ball_radius = ball_radius, ec_type = "SECT", subdivision = sphere_smoothness, cusps = num_cusps,
                                                 causal_regions_1 = c(1, 45, 55, 95, 100), causal_regions_2 = c(25, 145, 80), shared_regions = c(120))
  
  ### Prune the directions ###
  
  temp <- prune_directions_to_desired_number(data = sphere_data$data[,-1],directions, current_num_cones,curve_length,
                                             directions_per_cone,desired_num_cones)
  pruned_directions <- temp[[1]]
  ec_curve_data <- temp[[2]]
  ec_curve_data <- cbind(sphere_data$data[,1],ec_curve_data) # add the labels back on
  
  pruned_rate_values <- find_rate_variables_with_other_sampling_methods(ec_curve_data, radius = 0, bandwidth = 0.01,
                                                                        weights = TRUE, type = 'ESS')[, 2]
  

  ### create the ROC curve ###
  total_rate_roc <- matrix(0, nrow = length(pruned_rate_values), ncol = 2)
  vertex_region_dictionary <- sphere_data$vertex_region_dict
  
  # analyze only the class one
  for (j in seq(class, length(sphere_data$complex_points), 2)){
    
    sphere <- vcgSphere(subdivision = sphere_smoothness)
    sphere$vb[1:3, ] <- t(sphere_data$complex_points[[j]]) # put the cusps on the sphere
    complex <- convert_off_file(sphere)
    
    
    rate_ROC <- matrix(nrow = 0, ncol = 2)
    for (threshold in quantile(pruned_rate_values, probs = seq(1, 0, length.out = length(pruned_rate_values))) ){
      
      rate_positive_vertices <- compute_selected_vertices_cones(dir = pruned_directions, complex = complex, rate_vals = pruned_rate_values,
                                                                len = curve_length, threshold = threshold,
                                                                cone_size = directions_per_cone, ball_radius = 2)
      
      true_regions <- unique(vertex_region_dictionary[sphere_data$causal_points1])
      
      
      TPR_FPR <- calculate_TPR_FPR_closest_regions(rate_positive_vertices, true_regions, num_cusps, 
                                                   vertex_region_dictionary)
      
      rate_ROC <- rbind(rate_ROC, TPR_FPR)
    }
    
    total_rate_roc <- total_rate_roc + rate_ROC
    print(j)
    
  }
  
  ### standardize the averaged ROC ###
  total_rate_roc <- (total_rate_roc / max(total_rate_roc))
  total_rate_roc <- data.frame(total_rate_roc)
  colnames(total_rate_roc) <- c("FPR","TPR")
  
  ### Plot the ROC curve ###
  ROC_curve_plt <- ggplot(data <- data.frame(total_rate_roc), aes(x = FPR, y = TPR)) +
    geom_line(stat = "identity") +
    labs(x = "FPR", y = "TPR") +
    ggtitle(sprintf("initial_cones:%d,dirpercone:%d", initial_cones, directions_per_cone)) +
    geom_abline(intercept = 0, slope = 1)
  
  plot(ROC_curve_plt)
  
  ### Save the plots ###
  save_dir <-  sprintf("/Users/sianamuljadi/projects/Research/SINATRA/Results/Sphere_Sim_Results/voronoi_region_idea/parameter_variation/initial_cones:%d,dirpercone:%d.pdf",
                       initial_cones, directions_per_cone)
  ggsave(filename = save_dir)
  
  total_rate_roc
}

#################################################################
#################################################################
#################################################################
### Run the ROC parameter analysis ###

initial.cone.list <- seq(75, 125, by = 25)
dirpercone.list <- seq(5, 20, by = 5)

ROCs <- foreach(i = 1:length(initial.cone.list)) %:% 
          foreach(j = 1:length(dirpercone.list)) %do%{
              print(sprintf("%s,%s",initial.cone.list[i],dirpercone.list[j]))
              GenerateVoronoiRegionROC(initial.cone.list[i], dirpercone.list[j])
          }

