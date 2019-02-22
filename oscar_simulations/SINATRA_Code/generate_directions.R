######################################################################################
######################################################################################
######################################################################################

######## Direction Pruning Code ########

# assume that the cones include their central direction. Outputs the inter-class correlations and intra-class correlations
compute_total_correlations <- function(data,num_cones, curve_length, dir_per_cone){
  median_correlations <- 1:num_cones
  for (i in 0:(num_cones-1)){
    # test the correlations of only the central directions
    central_dir_index <- i*curve_length*dir_per_cone
    ec_curves <- data[,(central_dir_index+1):(central_dir_index + curve_length)]
    median_correlations[i+1] <- median(cor(t(ec_curves)))
  }
  return(median_correlations)
}

# input which cone to get- the function then computes the correlation matrix, both the interclass and intraclass correlation matrices
# works with central cone direction but can be easily modified.
compute_cone_class_correlations <- function(data, cone, curve_length, dir_per_cone){
  central_dir_index <- cone*curve_length*dir_per_cone
  num_data <- dim(data)[1]
  class_1_ec_curves <- data[seq(1,num_data,2),(central_dir_index+1):(central_dir_index + curve_length)]
  class_2_ec_curves <- data[seq(2,num_data,2),(central_dir_index+1):(central_dir_index + curve_length)]
  
  
  class_1_correlations <- cor(t(class_1_ec_curves))
  class_2_correlations <- cor(t(class_2_ec_curves))
  interclass_correlations <- cor(t(class_1_ec_curves),t(class_2_ec_curves))
  
  return(list(class1 = class_1_correlations, class2 = class_2_correlations, inter = interclass_correlations))
}

prune_directions_to_desired_number <- function(data, directions, num_cones, curve_length, dir_per_cone,desired_number){
  cors <- compute_total_correlations(data, num_cones, curve_length,dir_per_cone)
  
  # get the desired number of best cones by correlation, by removing the ones with high correlation
  idxs <- order(cors)[(desired_number+1):num_cones]
  idxs <- (idxs - 1)*dir_per_cone + 1
  
  return(update_data_and_directions(idxs,data,directions,curve_length,dir_per_cone))
}

# prune directions + data with correlations greater than 0.98, say.
# We might want to keep only a certain number; pick the k smallest ones?
prune_directions <- function(data, directions, num_cones, curve_length, dir_per_cone){
  cors <- compute_total_correlations(data, num_cones, curve_length,dir_per_cone)
  
  # get the central directions to remove, and then the associated directions in the cone
  idxs <- which(cors > 0.98)
  idxs <- (idxs - 1)*dir_per_cone + 1
  
  return(update_data_and_directions(idxs,data,directions,curve_length,dir_per_cone))
}

# prune directions that repeat already observed information, along with low variance directions - use interclass/intraclass metrics
# rank the directions with the most intraclass variance, least interclass variance
prune_low_var_and_collinear_directions <- function(data, directions, num_cones,curve_length, dir_per_cone){
  # rank the directions by least variance in a class, and most variance between classes. We can summarize the
  # associated correlation matrices by some statistic?
  
  
  
  # add variables that aren't collinear. Work with the central directions
}

# Helper function for updating data / directions, given the central cone indices to prune
update_data_and_directions <- function(idxs, data, directions, curve_length,dir_per_cone){
  # get the associated directions by shifting the idxs of the central directions 
  direction_idxs <- idxs
  for (i in 1:(dir_per_cone-1)){
    direction_idxs <- c(direction_idxs, idxs + i)
  }
  # remove the directions associated with each of these indices
  pruned_directions <- directions[-direction_idxs,]
  
  # Remove the data corresponding to these cones as well
  temp_idxs <- (direction_idxs - 1)*curve_length+1
  data_idxs <- temp_idxs
  for (i in 1:(curve_length-1)){
    data_idxs <- c(data_idxs, temp_idxs + i)
  }
  pruned_data <- data[,-data_idxs]
  
  return(list(pruned_directions,pruned_data))
}

######################################################################################
######################################################################################
######################################################################################

####### Generate Directions #######
generate_equidistributed_cones <- function(num_directions, cap_radius, directions_per_cone){
  # generate central directions that are equidistributed around the sphere.
  cones <- generate_equidistributed_points(num_directions)
  
  # renormalize these directions
  cones <- t(apply(cones, 1, function(x) x/sqrt(sum(x*x))))
  
  
  # generate directions for each cone
  directions <- matrix(0,ncol=3,nrow=0)
  for (i in 1:(dim(cones)[1]) ){
    directions <- rbind(directions,cones[i,])
    directions <- rbind(directions,rodriq(cones[i,],cap_radius,directions_per_cone-1))
  }
  
  directions
}
generate_random_cones <- function(num_directions,cap_radius,directions_per_cone){
  # uniformly generate random directions on the sphere
  cones <- matrix(rnorm(3*num_directions),ncol = 3)
  
  # renormalize these directions
  cones <- t(apply(cones, 1, function(x) x/sqrt(sum(x*x))))
  
  # generate directions for each cone
  directions <- matrix(0,ncol=3,nrow=0)
  for (i in 1:num_directions){
    directions <- rbind(directions,rodriq(cones[i,],cap_radius,directions_per_cone))
  }
  
  
  directions
}

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

# z is the direction around which we want to rotate
# z0 is a vector in the cone around z.
rodriq<-function(z,r,j){
  z<-z/sqrt(z[1]^2+z[2]^2+z[3]^2)
  if( z[1]*z[2]*z[3]==0){
    z0<-c(0,0,0)
    z0<-c(1,1,1)*(z[]==0)
  } 
  else{
    z0<-z
    z0[1]<- 2/z[1]
    z0[2]<- -1/z[2]
    z0[3]<- -1/z[3]
  }
  z0<-r*(z0/sqrt(z0[1]^2+z0[2]^2+z0[3]^2))
  z0<-z+z0
  z0<-z0/sqrt(z0[1]^2+z0[2]^2+z0[3]^2)
  B<-cross(z,z0)
  C<-z*(as.vector(z%*%z0))
  #print(B)
  #print(C)
  dir<-matrix(0,ncol=3,nrow=j)
  for (i in 1:j) {
    dir[i,]<-z0*cos(2*pi*i/j)+B*sin(2*pi*i/j)+C*(1-cos(2*pi*i/j))
    
  }
  return(dir)
}
cross<-function(x,y){
  a<-x
  a[1]<-x[2]*y[3]-x[3]*y[2]
  a[2]<-x[3]*y[1]-x[1]*y[3]
  a[3]<-x[1]*y[2]-x[2]*y[1]
  return(a)
}