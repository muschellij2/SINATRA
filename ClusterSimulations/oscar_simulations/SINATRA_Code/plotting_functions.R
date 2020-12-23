library(RColorBrewer)

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

##########################################################################################
##########################################################################################
##########################################################################################
### Get Vertex Persistence for reconstruction ###


### Code for smoothening functions ###


