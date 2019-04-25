library(Rcpp)
library(parallel)
library(Rvcg)
library(svd)

source("../SINATRA_Pipeline_Branch/ec_computation.R")
source("../SINATRA_Pipeline_Branch/generate_directions.R")
source("../SINATRA_Pipeline_Branch/mesh_functions.R")
source("../SINATRA_Pipeline_Branch/RATEv2.R")
source("../SINATRA_Pipeline_Branch/shape_reconstruction.R")
source("../SINATRA_Pipeline_Branch/gp_inference.R")
source("../SINATRA_Pipeline_Branch/plotting_functions.R")
sourceCpp("../SINATRA_Pipeline_Branch/BAKRGibbs.cpp") 
# Need gcc; fortran included as part of gcc in new mac OS # solved by creating .R/Makevars and 
#updating flibs variable = FLIBS=-L/opt/local/lib/gcc48/ - solves Rcpp Armadillo issue

### Compare really similar shapes
# Two classes of spheres with random noise on the surface, assign random classes


desired_num_cones <- 15
cap_radius <- 0.10
directions_per_cone <- 4


### Generate directions ###
dir <- generate_equidistributed_cones(desired_num_cones,cap_radius,directions_per_cone)

### Generate Data ###

nsim <- 50
curve_length <- 50
ball_radius <- 1.5
subdivision <- 3

data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]))

for (i in 1:nsim){
  
  sphere1 = vcgSphere(subdivision = subdivision)
  sphere2 = vcgSphere(subdivision = subdivision)
  
  # Add noise to the sphere
  sphere1$vb[1:3,] = sphere1$vb[1:3,]  * rnorm(dim(sphere1$vb)[2], mean = 1, sd = 0.02)
  sphere2$vb[1:3,] = sphere2$vb[1:3,]  * rnorm(dim(sphere2$vb)[2], mean = 1, sd = 0.02)
  
  sphere_mesh1 = convert_off_file(sphere1)
  sphere_mesh2 = convert_off_file(sphere2)

  ec_curve_class1 <- matrix(NA,nrow = 1,ncol=0)
  ec_curve_class2 <- matrix(NA,nrow = 1,ncol=0)
  
  ### compute EC curves for both classes of curves
  for (j in 1:dim(dir)[1]){
    
    vertex_function_class_1 <- sphere_mesh1$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])
    vertex_function_class_2 <- sphere_mesh2$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])
    
    curve1 <- compute_standardized_ec_curve(sphere_mesh1, vertex_function_class_1, curve_length-1, first_column_index = FALSE,ball_radius)
    curve2 <- compute_standardized_ec_curve(sphere_mesh2, vertex_function_class_2, curve_length-1, first_column_index = FALSE,ball_radius)
    
    # transform the ECT as desired
    curve1 <- update_ec_curve(curve1, "ECT")
    curve2 <- update_ec_curve(curve2, "ECT")
    
    # omit the length data, for now
    ec_curve_class1 <- c(ec_curve_class1,curve1[,2])
    ec_curve_class2 <- c(ec_curve_class2,curve2[,2])
  }
  
  data <- rbind(data,c(1,ec_curve_class1))
  data <- rbind(data,c(-1,ec_curve_class2))
  
}


### Run the model + select features with RATE
# how does bandwidth impact reconstruction? 
rate_values <- find_rate_variables_with_other_sampling_methods(data,radius = 0, bandwidth = 0.1,
                                                               weights = TRUE, type = 'ESS')[,2]


### Plot it back onto shape, and make rotating plot
sphere1 <- vcgSphere(subdivision = subdivision)
sphere1$vb[1:3,] <- sphere1$vb[1:3,]  * rnorm(dim(sphere1$vb)[2], mean = 1, sd = 0.035)
complex <- convert_off_file(sphere1)

# reconstruct birth times of vertices
quantiles <- seq(1,0,length.out = length(rate_values))
quantiled_rate <- quantile(rate_values,probs = quantiles)
order <- rep(1,length(complex$Vertices))
for (i in 1:length(quantiles)){

  reconstructed_vertices <- compute_selected_vertices_cones(dir, complex, rate_values,
                                                           curve_length, quantiled_rate[i],
                                                           directions_per_cone, ball_radius, TRUE)
  
  # calculate order of vertices
  for (v in reconstructed_vertices){
    order[v] <- min(order[v], quantiles[i])
  }
  
}



##########################################
##########################################
##########################################
### Compare really dissimilar shapes
# Class 1: spheres
# Class 2: Cubes?

### Generate directions ###
dir <- generate_equidistributed_cones(desired_num_cones,cap_radius,directions_per_cone)

### Generate Data ###

nsim <- 50
curve_length <- 50
ball_radius <- 1.5
subdivision <- 3

data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]))

for (i in 1:nsim){
  
  obj1 = vcgSphere(subdivision = subdivision)
  obj2 = vcgSphere(subdivision = subdivision)
  
  # Add noise to the sphere
  obj1$vb[1:3,] = obj1$vb[1:3,]  * rnorm(dim(obj1$vb)[2], mean = 1, sd = 0.02)
  obj2$vb[1:3,] = obj2$vb[1:3,]  * rnorm(dim(obj2$vb)[2], mean = 1, sd = 0.02)
  
  obj_mesh1 = convert_off_file(obj1)
  obj_mesh2 = convert_off_file(obj2)
  
  ec_curve_class1 <- matrix(NA,nrow = 1,ncol=0)
  ec_curve_class2 <- matrix(NA,nrow = 1,ncol=0)
  
  ### compute EC curves for both classes of curves
  for (j in 1:dim(dir)[1]){
    
    vertex_function_class_1 <- obj_mesh1$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])
    vertex_function_class_2 <- obj_mesh2$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])
    
    curve1 <- compute_standardized_ec_curve(obj_mesh1, vertex_function_class_1, curve_length-1, first_column_index = FALSE,ball_radius)
    curve2 <- compute_standardized_ec_curve(obj_mesh2, vertex_function_class_2, curve_length-1, first_column_index = FALSE,ball_radius)
    
    # transform the ECT as desired
    curve1 <- update_ec_curve(curve1, "ECT")
    curve2 <- update_ec_curve(curve2, "ECT")
    
    # omit the length data, for now
    ec_curve_class1 <- c(ec_curve_class1,curve1[,2])
    ec_curve_class2 <- c(ec_curve_class2,curve2[,2])
  }
  
  data <- rbind(data,c(1,ec_curve_class1))
  data <- rbind(data,c(-1,ec_curve_class2))
  
}


### Run the model + select features with RATE
# how does bandwidth impact reconstruction? 
rate_values <- find_rate_variables_with_other_sampling_methods(data,radius = 0, bandwidth = 0.1,
                                                               weights = TRUE, type = 'ESS')[,2]


### Plot it back onto shape, and make rotating plot
obj1 <- vcgSphere(subdivision = subdivision)
obj1$vb[1:3,] <- obj1$vb[1:3,]  * rnorm(dim(obj1$vb)[2], mean = 1, sd = 0.035)

obj1 <- vcgSphere(subdivision = subdivision)
obj1$vb[1:3,] <- obj1$vb[1:3,]  * rnorm(dim(obj1$vb)[2], mean = 1, sd = 0.035)

complex1 <- convert_off_file(obj1)
complex2 <- convert_off_file(obj2)

# reconstruct birth times of vertices
quantiled_rate <- quantile(rate_values,probs = seq(1,0,length.out = length(rate_values)))
order <- rep(length(rate_values),length(rate_values))
for (i in 1:len(quantiled_rate)){
  
  reconstructed_vertices1 <- compute_selected_vertices_cones(dir, complex1, rate_values,
                                                            curve_length, quantiled_rate[i],
                                                            directions_per_cone, ball_radius, TRUE)
  
  # calculate order of vertices
  for (v in reconstructed_vertices1){}
  order[]
  
}