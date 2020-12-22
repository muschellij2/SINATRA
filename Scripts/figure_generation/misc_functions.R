set.seed(55)
library(devtools)
install()

library(sinatra)
library(FNN)
library(rgl)
library(Rvcg)
library(plyr)
library(pdist)
library(gglasso)

#Parameters for the Analysis


mesh_to_matrix = function(mesh){
  v_points = mesh$vb[-4,]
  x = c(t(mesh$vb))
  return(x)
}
generate_data_sphere_simulation_baseline = function (nsim, curve_length, dir, noise_points = 5, causal_points = 5, 
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

generate_ROC_with_coned_directions_baseline <- function(nsim = 10, curve_length = 25, grid_size = 25, distance_to_causal_point = 0.1,
                                               causal_points = 10,shared_points = 3, num_cones = 5, eta = 0.1,
                                               truncated = 300, two_curves = TRUE, ball = TRUE, ball_radius = 2.5, type = 'vertex',
                                               min_points = 3,directions_per_cone = 4, cap_radius = 0.15, radius = 0,ec_type = 'ECT',
                                               mode = 'sphere',
                                               subdivision = 3,num_causal_region = 5, num_shared_region = 5){
  print("generating directions")
  
  
  print("generating data")
  # generate data
  
  if (mode == 'sphere'){
    #cusps = 2*num_causal_region + num_shared_region + 1
    cusps = 2*num_causal_region + num_shared_region + 1
    causal_dirs = generate_equidistributed_points(cusps,cusps)
    causal_regions_1 = sample(1:cusps,num_causal_region)
    causal_regions_2 = sample((1:cusps)[-causal_regions_1],num_causal_region)
    shared_regions = sample((1:cusps)[-c(causal_regions_1,causal_regions_2)],num_shared_region)
    directions <- generate_equidistributed_cones(num_cones,cap_radius,directions_per_cone)
    data <- generate_data_sphere_simulation_baseline(nsim = nsim,dir = directions, curve_length = curve_length,noise_points = shared_points,
                                            causal_points = causal_points, ball_radius = ball_radius, subdivision = subdivision,
                                            cusps = cusps, causal_regions_1 = causal_regions_1, causal_regions_2 = causal_regions_2,
                                            shared_regions = shared_regions, ec_type = ec_type)
    directions <- directions
    ec_curve_data <- data$data
    
  }
  num_cones <- dim(directions)[1]/directions_per_cone
  
  print("getting rate values")
  rate_values <- find_rate_variables_with_other_sampling_methods(ec_curve_data, bandwidth = 0.01, type = 'ESS')[,2]
  #Indices for Two Classes
  index1 = seq(1,nsim,2)
  complex_points1 = data$complex_points[index1]
  
  index2 = seq(2,nsim,2)
  complex_points2 = data$complex_points[index2]
  #Compute ROC using training data
  if (type == 'vertex'){
    if (two_curves == TRUE){
      roc_curve1 =  compute_roc_curve_vertex(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                             curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                             eta = eta, directions_per_cone = directions_per_cone, directions = directions, class = 1,truncated = truncated,
                                             ball_radius = ball_radius, radius = radius, mode = mode,subdivision = subdivision)
      roc_curve1 = cbind(roc_curve1, rep(1,dim(roc_curve1)[1]))
      roc_curve1 = cbind(roc_curve1,(1:dim(roc_curve1)[1]))
      
      roc_curve2 =  compute_roc_curve_vertex(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                             curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                             eta = eta, directions_per_cone = directions_per_cone, directions = directions,class = 2,truncated = truncated,
                                             ball_radius = ball_radius, radius = radius, mode = mode, subdivision = subdivision)
      roc_curve2 = cbind(roc_curve2, rep(2,dim(roc_curve2)[1]))
      roc_curve2 = cbind(roc_curve2,(1:dim(roc_curve2)[1]))
      
      roc_curve = rbind(roc_curve1,roc_curve2)
    }
    else{
      roc_curve = compute_roc_curve_vertex(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                           curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                           eta = eta, directions_per_cone = directions_per_cone, directions = directions,  class = 0, truncated = truncated,
                                           ball = ball, ball_radius = ball_radius, radius = radius, mode = mode, subdivision = subdivision)
      roc_curve = cbind(roc_curve,(1:dim(roc_curve)[1]))
    }
    return(roc_curve)
  }
  if (type == 'feature'){
    roc_curve = compute_roc_curve_features(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                           curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                           eta = eta, directions_per_cone = directions_per_cone, class = 0, truncated = truncated,
                                           ball = ball,ball_radius = ball_radius,
                                           dir = directions, min_points = min_points,mode = mode, subdivision = subdivision)
    roc_curve = cbind(roc_curve,(1:dim(roc_curve)[1]))
    return(roc_curve)
  }
  if (type == 'cone'){
    print('cone')
    roc_curve = compute_roc_curve_cones(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                        curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                        eta = eta, directions_per_cone = directions_per_cone, class = 0, truncated = truncated,
                                        ball = ball, ball_radius = ball_radius,
                                        dir = directions,  min_points = min_points, radius = radius,mode = mode, subdivision = subdivision)
    roc_curve = cbind(roc_curve,(1:dim(roc_curve)[1]))
    return(roc_curve)
  }
  if (type == 'cusp'){
    print('cusp')
    roc_curve = compute_roc_curve_modified_vertex(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                                  curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                                  eta = eta, directions_per_cone = directions_per_cone, class = 0, truncated = truncated,
                                                  ball = ball, ball_radius = ball_radius,
                                                  dir = directions, radius = radius,mode = mode, subdivision = subdivision)
    roc_curve = cbind(roc_curve,(1:dim(roc_curve)[1]))
    return(roc_curve)
  }
}