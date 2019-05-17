library(Rcpp)
library(rgl) # commented out for athena
library(ks)
library(mvtnorm)
library(pdist)
library(MASS)
library(ks)
library(truncnorm)
library(spatstat)
library(Rvcg)
library(doParallel)
library(svd)
library(plyr)
library(reshape2)
library(ggplot2)
library(FNN)

### Get the causal and shared regions on the sphere ###
cusps <- 50
causal_points <- 5
noise_points <- 5
causal_regions_1 <- c(23)
causal_regions_2 <- c(42)
shared_regions <- c(10, 45)
subdivision <- 3



regions <- generate_equidistributed_points(cusps,cusps)
sphere <- vcgSphere(subdivision = subdivision)

complex_points = list()
shared_points_list = list()
total_shared_points = c()
total_closest_points_class1 = c()
total_closest_points_class2 = c()
total_cusps_list = c()


# iterate through class 1
for (j in 1:length(causal_regions_1)){
  causal_dir1 = regions[causal_regions_1[j],]
  closest_points_class1 = knnx.index(data = t(sphere$vb[-4,]),query = matrix(causal_dir1,ncol = 3), k = causal_points)
  total_closest_points_class1 = c(total_closest_points_class1,closest_points_class1)
  total_cusps_list[[length(total_cusps_list) + 1]] = closest_points_class1
  
}

# iterate through class 2
for (j in 1:length(causal_regions_2)){
  causal_dir2 = regions[causal_regions_2[j],]
  closest_points_class2 = knnx.index(data = t(sphere$vb[-4,]),query = matrix(causal_dir2,ncol = 3), k = causal_points)
  total_closest_points_class2 = c(total_closest_points_class2,closest_points_class2)
  total_cusps_list[[length(total_cusps_list) + 1]] = closest_points_class2
}

for (k in 1:length(shared_regions)){
  shared_dir = regions[shared_regions[k],]
  closest_points_shared = knnx.index(data = t(sphere$vb[-4,]),query = matrix(shared_dir,ncol = 3), k = noise_points)
  total_shared_points = c(total_shared_points,closest_points_shared)
}

### Create the data ###
sphere1 = vcgSphere(subdivision = subdivision)
sphere2 = vcgSphere(subdivision = subdivision)

# Add noise to the sphere
sphere1$vb[1:3,] = sphere1$vb[1:3,]  * rnorm(dim(sphere1$vb)[2], mean = 1, sd = 0.035)
sphere2$vb[1:3,] = sphere2$vb[1:3,]  * rnorm(dim(sphere2$vb)[2], mean = 1, sd = 0.035)

# Elevate the causal regions - Needs to be changed
for (j in 1:length(causal_regions_1)){
  causal_dir1 = regions[causal_regions_1[j],]
  closest_points_class1 = knnx.index(data = t(sphere$vb[-4,]),query = matrix(causal_dir1,ncol = 3), k = causal_points)
  sphere1$vb[1:3,closest_points_class1] = sphere1$vb[1:3,closest_points_class1]  * 0.55 + rnorm(1, mean = 0, sd = 0.1)
}

for (j in 1:length(causal_regions_2)){
  causal_dir2 = regions[causal_regions_2[j],]
  closest_points_class2 = knnx.index(data = t(sphere$vb[-4,]),query = matrix(causal_dir2,ncol = 3), k = causal_points)
  sphere2$vb[1:3,closest_points_class2] = sphere2$vb[1:3,closest_points_class2]  * 0.55 + rnorm(1, mean = 0, sd = 0.1)
}

# Elevate the shared regions - Needs to be changed
for (k in 1:length(shared_regions)){
  shared_dir = regions[shared_regions[k],]
  closest_points_shared = knnx.index(data = t(sphere$vb[-4,]),query = matrix(shared_dir,ncol = 3), k = noise_points)
  shared_points = sphere$vb[1:3,closest_points_shared]  * 1.35 + rnorm(1, mean = 0, sd = 0.1)
  sphere1$vb[1:3,closest_points_shared] = shared_points
  sphere2$vb[1:3,closest_points_shared] = shared_points
}
  
mfrow3d(1,2)
cols1 = rep('white', dim(sphere1$vb)[2])
cols1[total_closest_points_class1] = 'red'
cols1[total_shared_points] = 'blue'
plot3d(sphere1,col = cols1,axes = FALSE, box = FALSE, xlab = '', ylab = '',zlab='')

cols2 = rep('white', dim(sphere2$vb)[2])
cols2[total_closest_points_class2] = 'red'
cols2[total_shared_points] = 'blue'
plot3d(sphere2,col = cols2,axes = FALSE, xlab = '', ylab = '',zlab='')

########################################################################
########################################################################
########################################################################
### Run the SINATRA Analysis ###

desired_num_cones <- 20
cap_radius <- 0.15
directions_per_cone <- 5

nsim <- 50
curve_length <- 20
ball_radius <- 1.5
subdivision <- 3
ec_type <- 'ECT'

### Generate directions ###
dir <- generate_equidistributed_cones(desired_num_cones,cap_radius,directions_per_cone)

### Generate Data ###
cusps <- 50

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


### Get the causal and shared regions on the sphere ###


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
    sphere1$vb[1:3,closest_points_class1] = sphere1$vb[1:3,closest_points_class1]  * 0.55 + rnorm(1, mean = 0, sd = 0.1)
  }
  
  for (j in 1:length(causal_regions_2)){
    causal_dir2 = regions[causal_regions_2[j],]
    closest_points_class2 = knnx.index(data = t(sphere$vb[-4,]),query = matrix(causal_dir2,ncol = 3), k = causal_points)
    sphere2$vb[1:3,closest_points_class2] = sphere2$vb[1:3,closest_points_class2]  * 0.55 + rnorm(1, mean = 0, sd = 0.1)
  }
  
  # Elevate the shared regions - Needs to be changed
  for (k in 1:length(shared_regions)){
    shared_dir = regions[shared_regions[k],]
    closest_points_shared = knnx.index(data = t(sphere$vb[-4,]),query = matrix(shared_dir,ncol = 3), k = noise_points)
    shared_points = sphere$vb[1:3,closest_points_shared]  * 1.35 + rnorm(1, mean = 0, sd = 0.1)
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


### Run the model + select features with RATE
# how does bandwidth impact reconstruction? 
rate_values_sim <- find_rate_variables_with_other_sampling_methods(data,radius = 0, bandwidth = 0.1,
                                                                   weights = TRUE, type = 'ESS')[,2]

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
vert_matrix1 <- reconstruct_vertices_on_shape(dir, complex1, rate_values_sim, curve_length, cuts = length(rate_values_sim),
                                              directions_per_cone, ball_radius, TRUE)

vert_matrix2 <- reconstruct_vertices_on_shape(dir, complex2, rate_values_sim, curve_length, cuts = length(rate_values_sim),
                                              directions_per_cone, ball_radius, TRUE)

# define heatmap colors
color1='blue'
color2='lightgreen'
color3='orangered'
color3 = 'red'
col_pal=c(color1,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)

# plot, using absolute birth times
rgl.close()
mfrow3d(2,1)
# vert_heat1 <- colfunc(cuts)[vert_matrix1[,1]] #absolute
vert_heat1 = colfunc(1 + max(vert_matrix1[,1]) - min(vert_matrix1[,1]))[1 + vert_matrix1[,1] - min(vert_matrix1[,1])] # relative
plot3d(sphere1, col = vert_heat1, back="lines", specular="white", axes = FALSE,xlab = '', ylab = '',zlab='')

#vert_heat2 <- colfunc(cuts)[vert_matrix2[,1]]
vert_heat2 = colfunc(1 + max(vert_matrix2[,1]) - min(vert_matrix2[,1]))[1 + vert_matrix2[,1] - min(vert_matrix2[,1])] # relative
plot3d(sphere2, col = vert_heat2, back="lines", specular="white", axes = FALSE,xlab = '', ylab = '',zlab='')

########################################################################
########################################################################
########################################################################
### Misc functions ###

generate_equidistributed_points <- function(desired_number, N){
  a <- 4*pi/N
  d <- sqrt(a)
  points <- matrix(NA,0,3) 
  
  M_theta <- round(pi/d)
  d_theta <- pi/M_theta
  d_phi <- a/d_theta
  for (i in 0:(M_theta-1)){
    theta <- pi*(i + 0.5)/M_theta
    M_phi <- round(2*pi*sin(theta)/d_phi)
    for (j in 0:(M_phi - 1)){
      phi <- 2*pi*j/M_phi
      point <- c( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))
      points <- rbind(points,point)
    }
  }
  if (dim(points)[1] < desired_number){
    return(generate_equidistributed_points(desired_number,N+1))
  }
  else{
    points = matrix(points[1:desired_number,],ncol = 3)
    return(points)
  }
}
  