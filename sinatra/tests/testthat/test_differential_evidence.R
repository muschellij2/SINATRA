context("Tests differential evidence implementation")
########################################################################
########################################################################
########################################################################
  ### Setup ###

  desired_num_cones <- 20
  cap_radius <- 0.15
  directions_per_cone <- 5

  nsim <- 15 # number of shapes in the data set to generate
  curve_length <- 25
  ball_radius <- 1.5
  ec_type <- 'ECT'

  # generate the set of directions to on which measure EC
  dir <- generate_equidistributed_cones(desired_num_cones,cap_radius,directions_per_cone)

  cusps <- 50
  subdivision <- 3 # granularity of the generated shape

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
  causal_regions_2 <- c(40)
  shared_regions <- c(25,10)

  # set the size of the causal / shared regions.
  causal_points <- 10
  noise_points <- 10


  data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]))

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

  rate_values_sim <- find_rate_variables_with_other_sampling_methods(data, bandwidth = 0.1, type = 'ESS')[,2]

  sphere1 <- vcgSphere(subdivision = subdivision)
  sphere1$vb[1:3,] <- sphere1$vb[1:3,]  * rnorm(dim(sphere1$vb)[2], mean = 1, sd = 0.02)

  for (j in 1:length(causal_regions_1)){
    causal_dir1 = regions[causal_regions_1[j],]
    closest_points_class1 = knnx.index(data = t(sphere$vb[-4,]),query = matrix(causal_dir1,ncol = 3), k = causal_points)
    sphere1$vb[1:3,closest_points_class1] = sphere1$vb[1:3,closest_points_class1]  * 1.55
  }
  for (k in 1:length(shared_regions)){
    shared_dir = regions[shared_regions[k],]
    closest_points_shared = knnx.index(data = t(sphere$vb[-4,]),query = matrix(shared_dir,ncol = 3), k = noise_points)
    shared_points = sphere1$vb[1:3,closest_points_shared]  * 0.55
    sphere1$vb[1:3,closest_points_shared] = shared_points
  }

  complex1<- convert_off_file(sphere1)

  cuts <-length(rate_values_sim)
  vert_matrix1 <- reconstruct_vertices_on_shape(dir, complex1, rate_values_sim, curve_length, cuts = length(rate_values_sim),
                                                directions_per_cone, ball_radius, TRUE)

  ########################################################################
  ########################################################################
  ########################################################################
  ### Test ###

  color1='blue'
  color2='lightgreen'
  color3='orangered'
  color3 = 'red'
  col_pal=c(color1,color2,color2,color3)
  colfunc <- colorRampPalette(col_pal)

  # plot, using absolute birth times
  vert_heat1 <- colfunc(cuts)[vert_matrix1[,1]] #absolute
  #vert_heat1 = colfunc(1 + max(vert_matrix1[,1]) - min(vert_matrix1[,1]))[1 + vert_matrix1[,1] - min(vert_matrix1[,1])] # relative
  plot3d(sphere1, col = vert_heat1, back="lines", specular="white", axes = FALSE,xlab = '', ylab = '',zlab='')

  num_null_regions <- 500
  x <- compute_differential_evidence(sphere1, vert_matrix1[,1], 15, 1, num_null_regions)

test_that("Differential Evidence Output is okay",{
  expect_equal(length(x), 1 + num_null_regions)
})


hist(x)
abline(v = x[1], col = 'red')
