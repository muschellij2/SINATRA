
### Take in arguments ###
#arguments <- commandArgs(trailingOnly=TRUE)
#print(arguments)
#
#### Ensure the proper arguments are inputted
#if (length(arguments) != 4) {
#  stop(sprintf("there are %s args",length(arguments)), call.=FALSE)
#}

### Load in the R Libraries ###
library(truncnorm)
library(doParallel)
library(rgl)
library(svd)
library(numbers)
library(sinatra)
library(Rvcg)
library(FNN)

### Load SINATRA functions ###
#source("SINATRA_Code/ec_computation.R")
#source("SINATRA_Code/generate_directions.R")
#source("SINATRA_Code/gp_inference.R")
#source("SINATRA_Code/metric_curve_simulation.R")
#source("SINATRA_Code/plotting_functions.R")
#source("SINATRA_Code/roc_curve_simulation.R")
#source("SINATRA_Code/shape_reconstruction.R")
#source("SINATRA_Code/simulated_data_generation.R")
#source("SINATRA_Code/RATEv2.R")
#sourceCpp("SINATRA_Code/BAKRGibbs.cpp")

### Set the parameters for the analysis ###
set.seed(4913, kind = "L'Ecuyer-CMRG")

arguments <- as.numeric(arguments)
nsim <- arguments[1]
num_cones <- arguments[2]
directions_per_cone <- arguments[3]
curve_length <- arguments[4]
nsim = 10
num_cones = 25
directions_per_cone = 10
curve_length = 50
cap_radius <- 0.15


######################################################
###################### Setup #####################
######################################################


### Generate directions ###
dir <- generate_equidistributed_cones(num_cones,cap_radius,directions_per_cone)

### Generate Data ###
cusps <- 50
subdivision <- 3
ball_radius <- 1.5
ec_type <- 'ECT'

### Create the Cusps on the sphere ###
regions =  generate_equidistributed_points(cusps,cusps)

#Initiate the causal points
sphere = vcgSphere(subdivision = subdivision)
region_vertex_dictionary <- vector("list",dim(regions)[1])

sphere_vertices <- asEuclidean(t(sphere$vb))

#get distances between regions and vertices
distances <- as.matrix(pdist::pdist(regions,sphere_vertices))

for (i in 1:(dim(sphere_vertices))[1]){
  closest_region <- which.min(distances[,i])
  region_vertex_dictionary[[closest_region]] <- c(region_vertex_dictionary[[closest_region]],i)
}

vertex_region_dictionary <- apply(distances,2,FUN = which.min)

causal_regions_1 <- c(1)
causal_regions_2 <- c(50)
shared_regions <- c(25,34)

causal_points <- 10
noise_points <- 10

### Start Timing
times <- 1:5
for(l in 1:5){
  t1 <- Sys.time()

  ######################################################
  ###################### Compute EC #####################
  ######################################################


  data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]))

  ### Create the data ###
  for (i in 1:nsim){
    sphere1 = vcgSphere(subdivision = subdivision)
    sphere2 = vcgSphere(subdivision = subdivision)

    # Add noise to the sphere
    sphere1$vb[1:3,] = sphere1$vb[1:3,]  * rnorm(dim(sphere1$vb)[2], mean = 1, sd = 0.02)
    sphere2$vb[1:3,] = sphere2$vb[1:3,]  * rnorm(dim(sphere2$vb)[2], mean = 1, sd = 0.02)

    # Elevate the causal regions - Needs to be changed
    for (j in 1:length(causal_regions_1)){
      causal_dir1 = regions[causal_regions_1[j],]
      closest_points_class1 = knnx.index(data = t(sphere$vb[-4,]),query = matrix(causal_dir1,ncol = 3), k = causal_points)
      sphere1$vb[1:3,closest_points_class1] = sphere1$vb[1:3,closest_points_class1]  * 1.55 + rnorm(1, mean = 0, sd = 0.1)
    }

    for (j in 1:length(causal_regions_2)){
      causal_dir2 = regions[causal_regions_2[j],]
      closest_points_class2 = knnx.index(data = t(sphere$vb[-4,]),query = matrix(causal_dir2,ncol = 3), k = causal_points)
      sphere2$vb[1:3,closest_points_class2] = sphere2$vb[1:3,closest_points_class2]  * 1.55 + rnorm(1, mean = 0, sd = 0.1)
    }

    # Elevate the shared regions - Needs to be changed
    for (k in 1:length(shared_regions)){
      shared_dir = regions[shared_regions[k],]
      closest_points_shared = knnx.index(data = t(sphere$vb[-4,]),query = matrix(shared_dir,ncol = 3), k = noise_points)
      shared_points = sphere$vb[1:3,closest_points_shared]  * 0.55 + rnorm(1, mean = 0, sd = 0.1)
      sphere1$vb[1:3,closest_points_shared] = shared_points
      sphere2$vb[1:3,closest_points_shared] = shared_points

    }


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
      curve1 <- update_ec_curve(curve1, ec_type)
      curve2 <- update_ec_curve(curve2, ec_type)

      # omit the length data, for now
      ec_curve_class1 <- c(ec_curve_class1,curve1[,2])
      ec_curve_class2 <- c(ec_curve_class2,curve2[,2])
    }

    data <- rbind(data,c(1,ec_curve_class1))
    data <- rbind(data,c(-1,ec_curve_class2))

  }

  ######################################################
  ###################### Apply RATE / Reconstruct #####################
  ######################################################


  ### Run the model + select features with RATE
  # how does bandwidth impact reconstruction?
  rate_values_sim <- find_rate_variables_with_other_sampling_methods(data,bandwidth = 0.1,type = 'ESS')[,2]

  ### Plot it back onto shape, and make rotating plot
  sphere1 <- vcgSphere(subdivision = subdivision)
  sphere1$vb[1:3,] <- sphere1$vb[1:3,]  * rnorm(dim(sphere1$vb)[2], mean = 1, sd = 0.02)

  ### Plot it back onto shape, and make rotating plot
  sphere2 <- vcgSphere(subdivision = subdivision)
  sphere2$vb[1:3,] <- sphere2$vb[1:3,]  * rnorm(dim(sphere2$vb)[2], mean = 1, sd = 0.02)

  for (j in 1:length(causal_regions_1)){
    causal_dir1 = regions[causal_regions_1[j],]
    closest_points_class1 = knnx.index(data = t(sphere$vb[-4,]),query = matrix(causal_dir1,ncol = 3), k = causal_points)
    sphere1$vb[1:3,closest_points_class1] = sphere1$vb[1:3,closest_points_class1]  * 1.55
  }

  for (j in 1:length(causal_regions_2)){
    causal_dir2 = regions[causal_regions_2[j],]
    closest_points_class2 = knnx.index(data = t(sphere$vb[-4,]),query = matrix(causal_dir2,ncol = 3), k = causal_points)
    sphere2$vb[1:3,closest_points_class2] = sphere2$vb[1:3,closest_points_class2]  * 1.55
  }

  for (k in 1:length(shared_regions)){
    shared_dir = regions[shared_regions[k],]
    closest_points_shared = knnx.index(data = t(sphere$vb[-4,]),query = matrix(shared_dir,ncol = 3), k = noise_points)
    shared_points = sphere1$vb[1:3,closest_points_shared]  * 0.55
    sphere1$vb[1:3,closest_points_shared] = shared_points
    shared_points = sphere2$vb[1:3,closest_points_shared]  * 0.55
    sphere2$vb[1:3,closest_points_shared] = shared_points
  }

  complex1<- convert_off_file(sphere1)
  complex2 <- convert_off_file(sphere2)

  # reconstruct birth times of vertices
  vert_matrix1 <- reconstruct_vertices_on_shape(dir, complex1, rate_values_sim, curve_length, cuts = 300,
                                                directions_per_cone, ball_radius, TRUE)

  vert_matrix2 <- reconstruct_vertices_on_shape(dir, complex2, rate_values_sim, curve_length, cuts = 300,
                                                directions_per_cone, ball_radius, TRUE)

  t2 <- Sys.time()

  simulation.time <- as.numeric(t2 - t1, units = "secs")
  times[l] <- simulation.time
}

avg_time <- mean(times)
stdev_time <- sd(times)

file <-"~/data/SINATRA/timingresults_sd.txt"
cat(sprintf("shapes: %d , num_cones: %d , dirpercone: %d , ec_curve: %d , time: %f, std_dev: %f ",
                   nsim, num_cones, directions_per_cone,
                   curve_length, avg_time, stdev_time), file = file, append = TRUE, sep = "\n")
