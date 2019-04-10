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
causal_regions_1 <- c(1,7,21)
causal_regions_2 <- c(3,10,15)
shared_regions <- c(4,25,30,39)
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
  
cols1 = rep('white', dim(sphere1$vb)[2])
cols1[total_closest_points_class1] = 'red'
cols1[total_shared_points] = 'blue'
plot3d(sphere1,col = cols1,axes = FALSE, box = FALSE, main = "Class 1 Simulated Shape")

cols2 = rep('white', dim(sphere2$vb)[2])
cols2[total_closest_points_class2] = 'red'
cols2[total_shared_points] = 'blue'
plot3d(sphere2,col = cols2,axes = FALSE, main = "Class 2 Simulated Shape")


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
  