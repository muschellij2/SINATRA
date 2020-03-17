set.seed(55)
library(sinatra)
library(FNN)
library(rgl)
library(Rvcg)
library(plyr)
library(pdist)
library(gglasso)

#Parameters for the Analysis

cap_radius = 0.25
num_cones = 5
directions_per_cone = 5
len = 75
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  cap_radius, directions_per_cone = directions_per_cone)
rotation_matrix=matrix(c(0.8065218,0.5911149,0.01028626,0,0.5583186,-0.7558216,0.34207344,0,0.1944301,-0.2816325,-0.93961692,0,0,0,0,1),ncol=4,byrow=TRUE)
ball = TRUE
ball_radius = 0.5
ec_type = 'ECT'

num_causal_region = 3
num_shared_region = 3
causal_points = 10
shared_points = 10

subdivision = 3

nsim = 25

cusps = 2*num_causal_region + num_shared_region + 1
causal_dirs = generate_equidistributed_points(cusps, cusps +1)
causal_regions_1 = sample(1:cusps,num_causal_region)
causal_regions_2 = sample((1:cusps)[-causal_regions_1],num_causal_region)
shared_regions = sample((1:cusps)[-c(causal_regions_1,causal_regions_2)],num_shared_region)


data = generate_data_sphere_simulation(nsim = nsim,dir = dirs, curve_length = len,noise_points = shared_points,
                                       causal_points = causal_points,ball_radius = ball_radius, subdivision = subdivision,
                                       cusps = cusps, causal_regions_1 = causal_regions_1, causal_regions_2 = causal_regions_2,
                                       shared_regions = shared_regions, ec_type = ec_type)


#### Plot Soln ####

mesh = vcgSphere(subdivision = 3)
mesh$vb[1:3,] = t(data$complex_points[[1]])
cols = rep('white', dim(mesh$vb)[2])
cols[data$causal_points1] = 'red'
cols[data$noise] = 'blue'
plot3d(mesh, col = cols)

#### Function ####
mesh_to_matrix = function(mesh){
  v_points = mesh$vb[-4,]
  x = c(t(mesh$vb))
  return(x)
}
generate_data_sphere_simulation_new = function (nsim, curve_length, dir, noise_points = 5, causal_points = 5, 
          ball = TRUE, ball_radius = 2, ec_type = "ECT", subdivision = 3, 
          cusps, causal_regions_1 = c(1), causal_regions_2 = c(3), 
          shared_regions = c(4)) 
{
  regions = generate_equidistributed_points(cusps, cusps)
  sphere = vcgSphere(subdivision = subdivision)
  print(paste("Causal Regions 1: "))
  print(causal_regions_1)
  print("Causal Regions 2: ")
  print(causal_regions_2)
  print("Shared Regions: ")
  print(shared_regions)
  complex_points = list()
  shared_points_list = list()
  total_shared_points = c()
  total_closest_points_class1 = c()
  total_closest_points_class2 = c()
  total_cusps_list = c()
  region_vertex_dictionary <- vector("list", dim(regions)[1])
  sphere_vertices <- asEuclidean(t(sphere$vb))
  distances <- as.matrix(pdist(regions, sphere_vertices))
  for (i in 1:(dim(sphere_vertices))[1]) {
    closest_region <- which.min(distances[, i])
    region_vertex_dictionary[[closest_region]] <- c(region_vertex_dictionary[[closest_region]], 
                                                    i)
  }
  vertex_region_dictionary <- apply(distances, 2, FUN = which.min)
  for (j in 1:length(causal_regions_1)) {
    causal_dir1 = regions[causal_regions_1[j], ]
    closest_points_class1 = knnx.index(data = t(sphere$vb[-4, 
                                                          ]), query = matrix(causal_dir1, ncol = 3), k = causal_points)
    total_closest_points_class1 = c(total_closest_points_class1, 
                                    closest_points_class1)
    total_cusps_list[[length(total_cusps_list) + 1]] = closest_points_class1
  }
  for (j in 1:length(causal_regions_2)) {
    causal_dir2 = regions[causal_regions_2[j], ]
    closest_points_class2 = knnx.index(data = t(sphere$vb[-4, 
                                                          ]), query = matrix(causal_dir2, ncol = 3), k = causal_points)
    total_closest_points_class2 = c(total_closest_points_class2, 
                                    closest_points_class2)
    total_cusps_list[[length(total_cusps_list) + 1]] = closest_points_class2
  }
  for (k in 1:length(shared_regions)) {
    shared_dir = regions[shared_regions[k], ]
    closest_points_shared = knnx.index(data = t(sphere$vb[-4, 
                                                          ]), query = matrix(shared_dir, ncol = 3), k = noise_points)
    total_shared_points = c(total_shared_points, closest_points_shared)
  }
  data <- matrix(NA, nrow = 0, ncol = (1+length(sphere$vb)))
  for (i in 1:nsim) {
    sphere1 = vcgSphere(subdivision = subdivision)
    sphere2 = vcgSphere(subdivision = subdivision)
    sphere1$vb[1:3, ] = sphere1$vb[1:3, ] * rnorm(dim(sphere1$vb)[2], 
                                                  mean = 1, sd = 0.035)
    sphere2$vb[1:3, ] = sphere2$vb[1:3, ] * rnorm(dim(sphere2$vb)[2], 
                                                  mean = 1, sd = 0.035)
    for (j in 1:length(causal_regions_1)) {
      causal_dir1 = regions[causal_regions_1[j], ]
      closest_points_class1 = knnx.index(data = t(sphere$vb[-4, 
                                                            ]), query = matrix(causal_dir1, ncol = 3), k = causal_points)
      sphere1$vb[1:3, closest_points_class1] = sphere1$vb[1:3, 
                                                          closest_points_class1] * 0.55 + rnorm(1, mean = 0, 
                                                                                                sd = 0.1)
    }
    for (j in 1:length(causal_regions_2)) {
      causal_dir2 = regions[causal_regions_2[j], ]
      closest_points_class2 = knnx.index(data = t(sphere$vb[-4, 
                                                            ]), query = matrix(causal_dir2, ncol = 3), k = causal_points)
      sphere2$vb[1:3, closest_points_class2] = sphere2$vb[1:3, 
                                                          closest_points_class2] * 0.55 + rnorm(1, mean = 0, 
                                                                                                sd = 0.1)
    }
    for (k in 1:length(shared_regions)) {
      shared_dir = regions[shared_regions[k], ]
      closest_points_shared = knnx.index(data = t(sphere$vb[-4, 
                                                            ]), query = matrix(shared_dir, ncol = 3), k = noise_points)
      shared_points = sphere$vb[1:3, closest_points_shared] * 
        1.35 + rnorm(1, mean = 0, sd = 0.1)
      sphere1$vb[1:3, closest_points_shared] = shared_points
      sphere2$vb[1:3, closest_points_shared] = shared_points
    }
    sphere_mesh1 = convert_off_file(sphere1)
    sphere_mesh2 = convert_off_file(sphere2)
    complex_points[[(2 * i - 1)]] = t(sphere1$vb[1:3, ])
    complex_points[[2 * i]] = t(sphere2$vb[1:3, ])
    shared_points_list[[i]] = shared_points
    ec_curve_class1 <- matrix(NA, nrow = 1, ncol = 0)
    ec_curve_class2 <- matrix(NA, nrow = 1, ncol = 0)
#    for (j in 1:dim(dir)[1]) {
#      vertex_function_class_1 <- sphere_mesh1$Vertices %*% 
#        c(dir[j, 1], dir[j, 2], dir[j, 3])
#      vertex_function_class_2 <- sphere_mesh2$Vertices %*% 
#        c(dir[j, 1], dir[j, 2], dir[j, 3])
    curve1 = mesh_to_matrix(sphere1)
    curve2 = mesh_to_matrix(sphere2)
    print(length(curve1))
      #if (ball == TRUE) {
      #  curve1 <- compute_standardized_ec_curve(sphere_mesh1, 
      #                                          vertex_function_class_1, curve_length - 1, 
      #                                          first_column_index = FALSE, ball_radius)
      #  curve2 <- compute_standardized_ec_curve(sphere_mesh2, 
      #                                          vertex_function_class_2, curve_length - 1, 
      #                                          first_column_index = FALSE, ball_radius)
      #}
      #else {
      #  curve1 <- compute_discrete_ec_curve(sphere_mesh1, 
      #                                      vertex_function_class_1, curve_length - 1, 
      #                                      first_column_index = FALSE)
      #  curve2 <- compute_discrete_ec_curve(sphere_mesh2, 
      #                                      vertex_function_class_2, curve_length - 1, 
      #                                      first_column_index = FALSE)
      #}
      #curve1 <- update_ec_curve(curve1, ec_type)
      #curve2 <- update_ec_curve(curve2, ec_type)
      #ec_curve_class1 <- c(ec_curve_class1, curve1[, 2])
      #ec_curve_class2 <- c(ec_curve_class2, curve2[, 2])
    ec_curve_class1 <- c(ec_curve_class1, curve1)
    ec_curve_class2 <- c(ec_curve_class2, curve2)
    print(length(ec_curve_class1))
    #}
    data <- rbind(data, c(1, ec_curve_class1))
    data <- rbind(data, c(-1, ec_curve_class2))
  }
  data_list = list(data = data, noise = total_shared_points, 
                   causal_points1 = total_closest_points_class1, causal_points2 = total_closest_points_class2, 
                   complex_points = complex_points, shared_points_list = shared_points_list, 
                   total_cusps_list = total_cusps_list, region_vertex_dict = region_vertex_dictionary, 
                   vertex_region_dict = vertex_region_dictionary)
  return(data_list)
}

m = mesh_to_matrix(mesh)

#### Test New Function ####
library(glmnet)
data2 = generate_data_sphere_simulation_new(nsim = nsim,dir = dirs, curve_length = len,noise_points = shared_points,
                                       causal_points = causal_points,ball_radius = ball_radius, subdivision = subdivision,
                                       cusps = cusps, causal_regions_1 = causal_regions_1, causal_regions_2 = causal_regions_2,
                                       shared_regions = shared_regions, ec_type = ec_type)
groups = rep(1:(dim(data2$data)[2]/3),each = 3)
#lasso = cv.gglasso(x = data2$data[,-1], y = data2$data[,1], group = groups, loss = "logit")
lasso = cv.glmnet(x = data2$data[,-1], y = data2$data[,1], alpha = 0.0,  family = "binomial")
coefs  = coef(lasso, s = 'lambda.min')

#Binary Option
nonzero = which(abs(coefs[-1])>0.8)
results = unique(groups[nonzero])


mesh = vcgSphere(subdivision = 3)
mesh$vb[1:3,] = t(data2$complex_points[[1]])
cols = rep('white', dim(mesh$vb)[2])
cols[data$causal_points1] = 'red'
cols[data$causal_points2] = 'red'
cols[data$noise] = 'blue'
cols[results] = 'green'
plot3d(mesh, col = cols)
plot3d(mesh, col = heat_colors)
#### Sanity Checks ####

j = cor(t(data_summary$data[,-1]))
heatmap(j)



#### Comp 2 ####
indices = c(261,295,833,224,939,278,252,593,177,357,4845,4215,4619,162,244,4618,4225,31,4705,4524,219,4301,4898,262,4860,69,64,4478,4452,4224,102,179,4176)

color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)

col_pal2=c(color1,color1,color2,color2,color3)
colfunc2 <- colorRampPalette(col_pal2)

