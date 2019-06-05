library(Rcpp)
library(parallel)
library(Rvcg)
library(rgl)
library(svd)
library(R.utils)

source("../athena_simulations/SINATRA_Code/ec_computation.R")
source("../athena_simulations/SINATRA_Code/generate_directions.R")
source("../athena_simulations/SINATRA_Code/mesh_functions.R")
source("../athena_simulations/SINATRA_Code/RATEv2.R")
source("../athena_simulations/SINATRA_Code/shape_reconstruction.R")
source("../athena_simulations/SINATRA_Code/gp_inference.R")
source("../athena_simulations/SINATRA_Code/plotting_functions.R")
sourceCpp("../athena_simulations/SINATRA_Code/BAKRGibbs.cpp")
# Need gcc; fortran included as part of gcc in new mac OS # solved by creating .R/Makevars and
#updating flibs variable = FLIBS=-L/opt/local/lib/gcc48/ - solves Rcpp Armadillo issue

### Helper functions
difference <- function(x,y){
  return( sqrt(sum((x-y)^2)) )
}


### Compare really similar shapes
# Two classes of spheres with random noise on the surface, assign random classes


desired_num_cones <- 20
cap_radius <- 0.15
directions_per_cone <- 10


### Generate directions ###
dir <- generate_equidistributed_cones(desired_num_cones,cap_radius,directions_per_cone)

### Generate Data ###

nsim <- 10
curve_length <- 25
ball_radius <- 2.5
subdivision <- 4

data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]))

for (i in 1:nsim){

  sphere1 = vcgSphere(subdivision = subdivision)
  sphere2 = vcgSphere(subdivision = subdivision)

  num_v <- (dim(sphere1$vb)[2])
  # Add noise to the sphere

  sphere2$vb[1:3,] = sphere2$vb[1:3,]  * rnorm(num_v, mean = 1, sd = 0.01)

  # draw a parabola on the sphere
  # find closest k points to (0,0,1) on the sphere
  # scale by distance away from the point
  sphere1$vb[1:3,] = sphere1$vb[1:3,] * rnorm(num_v, mean = 1, sd = 0.01)

  for (i in 1:num_v){
    #computes the 2D euclidean distance on the grid between the points
    dist1 = difference(dir[1,], sphere1$vb[1:3,i])
    if (dist1 < 0.2){
      sphere1$vb[1:3,i] <- sphere1$vb[1:3,i]*((1 + 10*(0.2 - dist1)))^0.5
    }
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
    curve1 <- update_ec_curve(curve1, "DECT")
    curve2 <- update_ec_curve(curve2, "DECT")

    # omit the length data, for now
    ec_curve_class1 <- c(ec_curve_class1,curve1[,2])
    ec_curve_class2 <- c(ec_curve_class2,curve2[,2])
  }

  data <- rbind(data,c(1,ec_curve_class1))
  data <- rbind(data,c(0,ec_curve_class2))

}


 ### Run the model + select features with RATE
# how does bandwidth impact reconstruction?
rate_values <- find_rate_variables_with_other_sampling_methods(data, bandwidth = 0.01, type = 'EP')[,2]

plot(rate_values)
### Plot it back onto shape, and make rotating plot
sphere1 <- vcgSphere(subdivision = subdivision)
sphere1$vb[1:3,] = sphere1$vb[1:3,] * rnorm(num_v, mean = 1, sd = 0.01)
neighbors <- c()

for (i in 1:num_v){
  #computes the 2D euclidean distance on the grid between the points
  dist1 = difference(dir[1,], sphere1$vb[1:3,i])
  if (dist1 < 0.2){
    sphere1$vb[1:3,i] <- sphere1$vb[1:3,i]*((1 + 10*(0.2 - dist1)))^0.5
    neighbors <- c(neighbors,i)
  }
}
complex <- convert_off_file(sphere1)

#### Test reconstruction of vertices
reconstructed_vs <- compute_selected_vertices_cones(dir, complex, rate_values, curve_length, 0.001,
                                directions_per_cone, ball_radius,
                                TRUE, 0)
print(length(reconstructed_vs))
cols = rep('white', dim(complex$Vertices)[1])
cols[reconstructed_vs ] <- 'green'
plot3d(sphere1, col = cols, back="lines", specular="white", axes = FALSE,xlab = '', ylab = '',zlab='')



#### reconstruct birth times of vertices
cuts <- 1000
vert_matrix <- reconstruct_vertices_on_shape(dir, complex, rate_values, curve_length, cuts = cuts,
                                             directions_per_cone, ball_radius, TRUE)

# define heatmap colors
color1='blue'
color2='lightgreen'
color3='orangered'
color3 = 'red'
col_pal=c(color1,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)

# plot, using absolute birth times
vert_heat <- colfunc(cuts)[vert_matrix[,1]]
plot3d(sphere1, col = vert_heat, back="lines", specular="white", axes = FALSE,xlab = '', ylab = '',zlab='')


col.table <- colfunc(cuts)
color.bar(col.table,min = 0, max = 100, nticks = 11)


####################################################################################
####################################################################################
####################################################################################
### Compare slightly dissimilar shapes
# Class 1: spheres
# Class 2: spheres with higher noise

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

  obj1 = vcgSphere(subdivision = subdivision)
  obj2 = vcgSphere(subdivision = subdivision)

  # Add noise to the sphere
  obj1$vb[1:3,] = obj1$vb[1:3,]  * rnorm(dim(obj1$vb)[2], mean = 1, sd = 0.02)
  obj2$vb[1:3,] = obj2$vb[1:3,]  * rnorm(dim(obj2$vb)[2], mean = 1, sd = 0.3)

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
obj1$vb[1:3,] <- obj1$vb[1:3,]  * rnorm(dim(obj1$vb)[2], mean = 1, sd = 0.02)

obj2 <- vcgSphere(subdivision = subdivision)
obj2$vb[1:3,] <- obj1$vb[1:3,]  * rnorm(dim(obj1$vb)[2], mean = 1, sd = 0.3)

complex1 <- convert_off_file(obj1)
complex2 <- convert_off_file(obj2)

# reconstruct birth times of vertices
cuts <- 1000
vert_matrix1 <- reconstruct_vertices_on_shape(dir, complex1, rate_values, curve_length, cuts = cuts,
                                             directions_per_cone, ball_radius, TRUE)

vert_matrix2 <- reconstruct_vertices_on_shape(dir, complex2, rate_values, curve_length, cuts = cuts,
                                              directions_per_cone, ball_radius, TRUE)

# define heatmap colors
color1='blue'
color2='lightgreen'
color3='orangered'
color3 = 'red'
col_pal=c(color1,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)

# plot, using absolute birth times
mfrow3d(1,2)
vert_heat1 <- colfunc(cuts)[cuts - vert_matrix1[,1]]
plot3d(obj1, col = vert_heat1, back="lines", specular="white", axes = FALSE,xlab = '', ylab = '',zlab='')

vert_heat2 <- colfunc(cuts)[cuts - vert_matrix2[,1]]
plot3d(obj2, col = vert_heat2, back="lines", specular="white", axes = FALSE,xlab = '', ylab = '',zlab='')


# heatmap
col.table <- colfunc(cuts)
color.bar(col.table,min = 0, max = 100, nticks = 11)




####################################################################################
####################################################################################
####################################################################################
### Compare very dissimilar objects
# Class 1: mushrooms
# Class 2: spheres with higher noise

# read in off file
dragon <- vcgImport("nullcomparisonshapes/dragon.off")
mushroom <- vcgImport("nullcomparisonshapes/mushroom.off")

# perturb surface with some noise
data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]))

for (i in 1:nsim){

  obj1 <- dragon
  obj2 <- mushroom

  # Add noise to the objects
  obj1$vb[1:3,] <- obj1$vb[1:3,] + rnorm(dim(obj1$vb)[2], mean = 1, sd = 0.02)
  obj2$vb[1:3,] <- obj2$vb[1:3,] + rnorm(dim(obj2$vb)[2], mean = 1, sd = 0.02)

  obj_mesh1 <- convert_off_file(obj1)
  obj_mesh2 <- convert_off_file(obj2)

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
obj1 <- dragon
obj2 <- mushroom

complex1 <- convert_off_file(obj1)
complex2 <- convert_off_file(obj2)

# reconstruct birth times of vertices
cuts <- 1000
vert_matrix1 <- reconstruct_vertices_on_shape(dir, complex1, rate_values, curve_length, cuts = cuts,
                                              directions_per_cone, ball_radius, TRUE)

vert_matrix2 <- reconstruct_vertices_on_shape(dir, complex2, rate_values, curve_length, cuts = cuts,
                                              directions_per_cone, ball_radius, TRUE)

# define heatmap colors
color1='blue'
color2='lightgreen'
color3='orangered'
color3 = 'red'
col_pal=c(color1,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)

# plot, using absolute birth times
mfrow3d(1,2)
vert_heat1 <- colfunc(cuts)[cuts - vert_matrix1[,1]]
plot3d(obj1, col = vert_heat1, back="lines", specular="white", axes = FALSE,xlab = '', ylab = '',zlab='')

vert_heat2 <- colfunc(cuts)[cuts - vert_matrix2[,1]]
plot3d(obj2, col = vert_heat2, back="lines", specular="white", axes = FALSE,xlab = '', ylab = '',zlab='')

