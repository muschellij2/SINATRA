library(RColorBrewer)
library(Rvcg)
library(rgl)
library(FNN)
library(pracma)
library(matlib)
library(far)
library(rgl)

convert_off_file = function(mesh){
  vertices=as.matrix(t(mesh$vb)[,1:3])
  faces=as.matrix(t(mesh$it))
  edges=vcgGetEdge(mesh)
  edges=as.matrix(edges[,1:2])
  complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
  return(complex)
}


##########################################################################################
##########################################################################################
##########################################################################################
### Heatmap Code ###

# map rate_values to vertices; output is a list.
vertex_to_feature_dict <- function(complex, dirs, curve_length, ball_radius){
  num_vertices <- dim(complex$Vertices)[1]
  f <- functional::Curry(vertex_to_projections, dirs = dirs, curve_length = curve_length, ball_radius = ball_radius)
  
  t(apply(complex$Vertices, 1, f))
}

# maps a vertex to the projections in all directions 
vertex_to_projections <- function(vertex, dirs, curve_length,ball_radius){
  f <- functional::Curry(vertex_to_dir_projection, vertex = vertex, dirs = dirs,
                         curve_length = curve_length, ball_radius = ball_radius)
  
  purrr::map_int(1:(dim(dirs)[1]), f)
}

# idx represents the index in the set of directions
# returns for each vertex in the complex, the index of the projection in the given direction
vertex_to_dir_projection <- function(vertex, idx, dirs, curve_length, ball_radius){
  vtx_projection <- vertex%*%dirs[idx,]
  buckets <- seq(-ball_radius,ball_radius,length.out = curve_length+1)
  
  # map vertex projection to the feature index
  projection_bucket <- cut(vtx_projection, buckets, labels = FALSE)
  
  # update index to reflect rate values
  projection_bucket <- as.integer(projection_bucket + (idx - 1)*curve_length)
  
  projection_bucket
}


# color the vertices of the complex based on how many RATE-selected projections it lies in
get_projection_heatmap <- function(complex, directions, rate_values, 
                                   curve_length,ball_radius, threshold = 1/length(rate_values)){
  rate_selected_features <- which(rate_values > threshold)
  vertex_to_feature <- vertex_to_feature_dict(complex, directions, curve_length, ball_radius)
  
  # each entry in the list is
  vertex_to_num_projections <- matrix(0,ncol = dim(complex$Vertices)[1])
  for (i in 1:dim(complex$Vertices)[1]){
    vertex_to_num_projections[i] = sum(vertex_to_feature[i,] %in% rate_selected_features)
  }
  vertex_to_num_projections
}


#Figures out at which importance are the vertices reconstructed
reconstruct_vertices_on_shape = function(dir, complex, rate_vals, len, cuts=10, cone_size, ball_radius,ball = TRUE){
  vert_matrix = matrix(0,nrow = dim(complex$Vertices)[1], ncol = 2)
  cut = cuts
  reconstructed_vertices = c()
  if (ball == TRUE){
    for (threshold in quantile(rate_vals,probs = seq(1,0,length.out = cuts)) ){
      selected_vertices = compute_selected_vertices_cones(dir = dir, complex = complex, rate_vals = rate_vals, len = len, threshold = threshold,
                                                          cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
      selected_vertices = setdiff(selected_vertices,reconstructed_vertices)
      vert_matrix[selected_vertices,1] = cut
      vert_matrix[selected_vertices,2] = threshold
      cut = cut-1
      reconstructed_vertices = c(reconstructed_vertices,selected_vertices)
    }}
  else{
    for (threshold in quantile(rate_vals,probs = seq(1,0,length.out = cuts)) ){
      selected_vertices = compute_selected_vertices_cones(dir = dir, complex = complex, rate_vals = rate_vals, len = len, threshold = threshold,
                                                          cone_size = directions_per_cone, ball = FALSE)
      selected_vertices = setdiff(selected_vertices,reconstructed_vertices)
      vert_matrix[selected_vertices,1] = cut
      vert_matrix[selected_vertices,2] = threshold
      cut = cut-1
      reconstructed_vertices = c(reconstructed_vertices,selected_vertices)
    }
  }
  return(vert_matrix)
}

# color the vertices of the complex based on sum of RATE values on each projection it lies in
get_RATE_weighted_heatmap <- function(complex, directions, rate_values,
                                      curve_length,ball_radius, threshold = 1/length(rate_values)){
  rate_selected_features <- which(rate_values > threshold)
  vertex_to_feature <- vertex_to_feature_dict(complex, directions, curve_length, ball_radius)
  
  # each entry in the list is
  vertex_to_rate_weighted_projections <- matrix(0,ncol = dim(complex$Vertices)[1])
  for (i in 1:dim(complex$Vertices)[1]){
    vertex_to_rate_weighted_projections[i] = sum(rate_values[which(vertex_to_feature[i,] %in% rate_selected_features)])
  }
  vertex_to_rate_weighted_projections
}

# Wrapper for plotting teeth within each group #

plot_results_teeth_simple=function(files, features1, features2, color1, color2, alpha1=0.65, alpha2=0.65,
                                   dir1, dir2, len=65, level=20, slices=25, n=10, thresh = 1.00,
                                   directions_per_cone = 5){
  #Loop through files in the directory
  for (i in 1:length(files)){
    #We set direction as (0,0.75,1)
    file1=vcgImport(files[i])
    file_1=process_off_file_v3(files[i])
    #Try to see if this file has any critical points, if it doesn't just plot the tooth.
    vert1 = compute_selected_vertices_cones(dir = dir1, complex =file_1, rate_vals = features1, len = len, threshold = (thresh/(length(features1))),
                                            cone_size = directions_per_cone,ball_radius = ball_radius, ball = FALSE)
    vert2 =  compute_selected_vertices_cones(dir = dir2, complex =file_1, rate_vals = features2, len = len, threshold = (thresh/(length(features2))),
                                             cone_size = directions_per_cone,ball_radius = ball_radius, ball = FALSE)
    intersected = intersect(vert1,vert2)
    fc3 <- colorRampPalette(c(color1,color2))
    colors = rep('white', dim(veg1$vb)[2])
    colors[setdiff(vert1,vert2)] =fc3(10)[2]
    colors[setdiff(vert2,vert1)] =fc3(10)[9]
    colors[intersected] = fc3(10)[6]
    plot3d(file1, colors = colors,axes = FALSE, xlab = '',ylab = '', zlab = '')
    #Rotate the tooth for a view of the 'bottom'
    rotation_matrix=matrix(c(0.99972576,0.02127766,0.00978078,0,0.01589701,-0.92330807,0.38373107,0,0.01719557,-0.38347024,-0.92339307,0,0,0,0,1),ncol=4,byrow=TRUE)
    rgl.viewpoint(userMatrix = rotation_matrix)
  }
}

##########################################################################################
##########################################################################################
##########################################################################################
### Get Vertex Persistence for reconstruction ###


### Code for smoothening functions ###


