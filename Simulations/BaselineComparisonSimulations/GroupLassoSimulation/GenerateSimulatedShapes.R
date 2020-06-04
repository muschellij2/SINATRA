library(rgl)
library(Rvcg)
library(FNN)
library(rdist)

get_euclidean_fps_landmarks = function(mesh, num_landmarks){
  ### Gets landmarks using farthest point sampling
  landmarks = rdist::farthest_point_sampling(mesh$Vertices, metric = 'Euclidean')
  landmarks
}

# not implemented yet
get_geodesic_fps_landmarks = function(mesh, num_landmarks){
  
}

mesh_to_matrix = function(mesh){
  v_points = mesh$vb[-4,]
  x = c(t(v_points))
  return(x)
}

### Need to make the return a list:
# Need data, complex points, causal_points1, causal_points2
generate_landmark_data_sphere_simulation = function (nsim, noise_region_size = 5, causal_region_size = 5, 
                                                     subdivision = 3, 
                                                     num_shape_scaffold, num_landmarks,
                                                     causal_regions_1 = c(1), causal_regions_2 = c(3), 
                                                     shared_regions = c(4)) 
{
  landmarks = generate_equidistributed_points(num_landmarks, num_landmarks) # sample landmarks on the sphere. If the landmark
    # is contained in a causal region, it's a true positive.
  shape_scaffold = generate_equidistributed_points(num_shape_scaffold, num_shape_scaffold) ## cusps is the number of  points to generate the shape

  print(paste("Causal Regions 1: "))
  print(causal_regions_1)
  
  print("Causal Regions 2: ")
  print(causal_regions_2)
  
  print("Shared Regions: ")
  print(shared_regions)
  
  #initialize vector to hold the data 
  data <- matrix(NA, nrow = 0, ncol = (1+3*num_landmarks))
  
  
  
  # generate the shapes
  for (i in 1:nsim) {
    sphere = vcgSphere(subdivision = subdivision)
    sphere1 = vcgSphere(subdivision = subdivision)
    sphere2 = vcgSphere(subdivision = subdivision)
    
    data_class1_indices <- knnx.index(data = t(sphere1$vb[-4,]),
                                      query = matrix(landmarks, ncol = 3), k = 1) 
    # for each region, find the closest point on the sphere
    
    data_class2_indices <- knnx.index(data = t(sphere2$vb[-4,]),
                                      query = matrix(landmarks, ncol = 3), k = 1)
    
    # which one of these regions are the causal points?
    
    
    # perturb surface of the sphere with some noise
    sphere1$vb[1:3, ] = sphere1$vb[1:3, ] * rnorm(dim(sphere1$vb)[2], 
                                                  mean = 1, sd = 0.035)
    sphere2$vb[1:3, ] = sphere2$vb[1:3, ] * rnorm(dim(sphere2$vb)[2], 
                                                  mean = 1, sd = 0.035)
    
    class1_positive_landmarks = c()
    class2_positive_landmarks = c()
    
    for (j in 1:length(causal_regions_1)) {
      causal_region1 = shape_scaffold[causal_regions_1[j], ]
      
      # find the points on the  sphere complex to perturb (in order to create the causal region)
      closest_points_class1 = knnx.index(data = t(sphere$vb[-4,]),
                                         query = matrix(causal_region1, ncol = 3), k = causal_region_size)
      
      # use the above to find the landmarks which are true positives for class1
      class1_positive_landmarks = c(class1_positive_landmarks,
                                    which(data_class1_indices %in% closest_points_class1))
      
      
      sphere1$vb[1:3, closest_points_class1] = sphere$vb[1:3,
                                                          closest_points_class1] * 0.55 + rnorm(1, mean = 0, 
                                                                                                sd = 0.1)
    }
    
    for (j in 1:length(causal_regions_2)) {
      causal_region2 = shape_scaffold[causal_regions_2[j], ]
      closest_points_class2 = knnx.index(data = t(sphere$vb[-4, 
                                                            ]), query = matrix(causal_region2, ncol = 3), k = causal_region_size)
      
      
      class2_positive_landmarks = c(class2_positive_landmarks,
                                    which(data_class2_indices %in% closest_points_class2))
      
      sphere2$vb[1:3, closest_points_class2] = sphere$vb[1:3, 
                                                          closest_points_class2] * 0.55 + rnorm(1, mean = 0, 
                                                                                                sd = 0.1) 
      
    }
    
    for (j in 1:length(shared_regions)) {
      shared_dir = shape_scaffold[shared_regions[j], ]
      closest_points_shared = knnx.index(data = t(sphere$vb[-4, 
                                                            ]), query = matrix(shared_dir, ncol = 3), k = noise_region_size)
      shared_points = sphere$vb[1:3, closest_points_shared] * 
        1.35 + rnorm(1, mean = 0, sd = 0.1) #1.35; this noise should be different for each row?
      sphere1$vb[1:3, closest_points_shared] = shared_points
      sphere2$vb[1:3, closest_points_shared] = shared_points
    }
    
    # basically consider the shapes generated at this point.
    # change data_class1_indices to landmarks that are contained in data_class1.
    
    # get perturbed positions of points
    data_class1 <- c(sphere1$vb[-4,data_class1_indices])
    data_class2 <- c(sphere2$vb[-4,data_class2_indices])
    
    data <- rbind(data, c(1, data_class1))
    data <- rbind(data, c(-1, data_class2))
  }
  
  # should return landmarks, data, which landmarks are positives.
  return(list(data,class1_positive_landmarks, class2_positive_landmarks ))
  # return list of landmarks that are true positive as well.
}

### Also generate a function that creates an averaged ROC curve,



### Generates the ROC function with the lasso procedure.
# Need to hack the rate_values so that it can take in the group lasso coefs.
generate_ROC_baseline = function(nsim = 20, num_shape_scaffold = 100, num_landmarks = 1000, causal_region_size = 40, noise_region_size = 40, 
                                    subdivision = 4, num_causal_region = 2, num_shared_region = 3){
  print("generating data")
  
  causal_regions_1 = sample(1:num_shape_scaffold,num_causal_region)
  causal_regions_2 = sample((1:num_shape_scaffold)[-causal_regions_1],num_causal_region)
  shared_regions = sample((1:num_shape_scaffold)[-c(causal_regions_1,causal_regions_2)],num_shared_region)
  
  shape_data = generate_landmark_data_sphere_simulation(nsim = nsim, subdivision = subdivision, 
                                                  causal_region_size = causal_region_size, noise_region_size = noise_region_size, 
                                                  num_shape_scaffold = num_shape_scaffold, num_landmarks = num_landmarks,
                                                  causal_regions_1 = causal_regions_1, causal_regions_2 = causal_regions_2,
                                                   shared_regions = shared_regions)

  
  data = shape_data[[1]]
  class1_positive_landmarks = shape_data[[2]]
  class2_positive_landmarks = shape_data[[3]]
  
  print("Getting Lasso Coefficients")
  groups = rep(1:(dim(data[,-1])[2]/3), each = 3)
  lasso = cv.glmnet(x = data[,-1], y = data[,1], alpha = 1.0,  family = "binomial") # how do these work?
  coefs = coef(lasso, s = 'lambda.min') # it might have  to do with  this choice of lasso coefficient.
  
  
  ### Need to organize the data to get the proper landmark chosen.
  #Compute ROC using training data
  

  roc_curve = compute_roc_curve_landmark(lasso_coefs = coefs,
                                         causal_points1 = class1_positive_landmarks,
                                         causal_points2 = class2_positive_landmarks)
  roc_curve = cbind(roc_curve,(1:dim(roc_curve)[1]))

  
  return(roc_curve)
}  


compute_roc_curve_landmark = function(lasso_coefs, causal_points1, causal_points2){
  
  num_landmarks = as.integer(length(lasso_coefs)/3)
  groups <- rep(1:(length(lasso_coefs)/3), each = 3)
  
  ROC <- matrix(0,nrow = 0,ncol = 2)
  
  for (threshold in quantile(abs(lasso_coefs),probs = seq(1,0,length.out = length(lasso_coefs))) ){

    chosen_coordinates = which(abs(lasso_coefs) >= threshold)
    
    positive_landmarks = unique(groups[chosen_coordinates])
    negative_landmarks = setdiff(1:num_landmarks,positive_landmarks)
    
    combined_true_landmarks = union(causal_points1, causal_points2)
    combined_false_landmarks = setdiff(1:num_landmarks, combined_true_landmarks)

    ROC <- rbind(ROC, calculate_TPR_FPR(positive_landmarks,negative_landmarks,
                                                  combined_true_landmarks,combined_false_landmarks))

  }
  
  return(ROC)
}

generate_averaged_ROC_baseline = function(runs = 10, nsim = 50, num_shape_scaffold = 100, num_landmarks = 1000, causal_region_size = 40, noise_region_size = 40, 
                                          subdivision = 4, num_causal_region = 2, num_shared_region = 3){

  total_roc = matrix(0, 1+3*num_landmarks, ncol = 2)
  
  for (j in 1:runs){
    roc_curve = generate_ROC_baseline(nsim = nsim, num_shape_scaffold = num_shape_scaffold , num_landmarks = num_landmarks,
                                      causal_region_size = causal_region_size, noise_region_size = noise_region_size,
                                      subdivision = subdivision,
                                      num_causal_region = num_causal_region, num_shared_region = num_shared_region)
    
    total_roc = total_roc + as.matrix(roc_curve)[,1:2]
  }
  
  total_roc = total_roc/runs
  return(total_roc)
  
}
