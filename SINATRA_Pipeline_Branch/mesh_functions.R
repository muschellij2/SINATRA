# for every direction, compute the euler chacteristic of a 3D complex, filtered by direction. 
#write dumb version that recomputes ec for every complex in the filtration
# more advanced version that computes the ec based on updates to the complex in the filtration.

# Complex (for simplicial complex) should contain vertices, faces, edges, and t? as in the EC3D.R file.
# We assume that Complex has a field Vertices such that complex$Vertices is a n x 3 matrix such that every row 
# gives the coordinates of the vertices. Complex$Edges should be a list of the edges, m x 2, where each row contains
# the vertices belonging to each edge. Complex$Faces has the same format.

# How to specify direction on the sphere? potentially use spherical coordinates, or unit vector coordinates.
# function should just be dot product of each point's coordinates with this unit vector, perhaps upon centering first.

# integrate the euler characteristic to get a curve.

library(geometry)
library(misc3d)
library(rgl)
library(Rvcg)

# For point cloud data; one option for manual nearest neighbor type thing, or just use vietoris rips
# check https://stackoverflow.com/questions/22630226/3d-surface-plot-with-xyz-coordinates
create_simplicial_complex <- function(){
  # rely on alpha shape package; many other functions for surface reconstruction.
  # might need to play around with alpha
  ashape3d.obj <- ashape3d(x,alpha = 0.5)
  plot(ashape3d.obj) # to visualize
  
  vertices <- 1:length(x)
  edges <- ashape3d.obj$edge[which(ashape3d.obj$edge[,4] == 1),1:2] # column 4 of this object returns which edges are attached
  faces <- ashape3d.obj$triang[which(ashape3d.obj$triang[,4] == 1),1:3] 
  
  
  complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
}


# We assume the input is a text file in a similar format as an OFF file (not the same)
# The first line should say #vertices #faces #edges; the following lines are the vertices coordinates in each new line;
# following are the edges, with the indices of the vertices in the edges. The Following is are the vertices in each face, 
#preceded by the number of vertices in that face. See the example text file.
process_text_file <- function(fileName){
  con = file(fileName, "r") # con for connection
  lines = readLines(con) # apparently this is dangerous with big files
  close(con)
  split_line = strsplit(lines[1]," ")[[1]] #this is a list, we index the first element (a vector) with the double bracket 
  num_vertices = strtoi(split_line[1])
  num_edges = strtoi(split_line[3])
  num_faces = strtoi(split_line[2])
  n = num_vertices + num_edges + num_faces
  dim = length(strsplit(lines[2]," ")[[1]]) # get the dimension of the points
  face_size = strtoi(strsplit(lines[2 + num_vertices + num_edges], " ")[[1]])[1]
  
  
  # store the vertices
  vertices <- matrix(NA,nrow = num_vertices,ncol = dim)
  for (i in 2:(num_vertices+1)){
    vertices[i-1,] <- strtoi(strsplit(lines[i], " ")[[1]]) # coordinates of each vertex
  }
  
  #store the edges
  edges <- matrix(NA,nrow = num_edges,ncol = 2)
  for (i in (num_vertices+2):(num_vertices+num_edges+1)){
    edges[i - num_vertices - 1,] <- strtoi(strsplit(lines[i], " ")[[1]]) + 1# + 1 if the vertices are zero indexed
  }
  
  #store the faces
  faces <- matrix(NA, nrow = num_faces, ncol = face_size)
  for (i in (num_vertices+num_edges+2):(n+1)){
    faces[i - num_vertices - num_edges - 1,] <- 
      strtoi(strsplit(lines[i], " ")[[1]])[1:face_size+1]  +1 # + 1 if the vertices are zero indexed
  }
  
  complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
}
#Processing OFF files specifically

process_off_file_v3=function(input_dir){
  off=vcgImport(input_dir,silent = TRUE)
  vertices=as.matrix(t(off$vb)[,1:3])
  faces=as.matrix(t(off$it))
  edges=vcgGetEdge(off)
  edges=as.matrix(edges[,1:2])
  complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
  return(complex)
}

# The complex is filtered by some function on the vertices of the simplicial complex; in this case, we use the distance
# to some plane in a certain direction in S2 (the two dimensional sphere). To compute this distance in R2, we simply multiply 
# coordinates of the vertices by a rotation matrix. In three dimensions, we dot product with the unit vector in the desired 
# direction.

#step size is how fine we want the ec curve to be.
# at the last step, we center the curve and integrate it to obtain the smooth euler characetristic curve transform (SECT)
# gEuler achieves this.

#The Euler characteristic of a subcomplex which in our case is a a surface of a polyhedra has a simple form: V - E + F
#curve_length outputs the ec curve as a vector of the same size. 
#vertex_function are the function values on the vertices
#SECT - smooth Euler Characteristic Transform
# returns approximately num_directions
#### Processing OFF files from a Directory

create_ec_matrix_mult_d=function(directory,directions,len,ball_radius = 1, ball = FALSE,ec_type = 'ECT'){
  curve_length = len
  file_names=list.files(path=directory,full.names = TRUE)
  number_files=length(file_names)
  data <- matrix(NA,nrow=number_files,ncol = curve_length*dim(directions)[1])
  for (i in 1:number_files){
    print(paste('On File', i))
    off=process_off_file_v3(file_names[i])
    curve_mult_d <- matrix(NA,nrow = 1,ncol=0)
    if (ball == FALSE){
      for (j in 1:dim(directions)[1]){
        vertex_function=off$Vertices%*%directions[j,]
        curve <- compute_discrete_ec_curve(off, vertex_function, len-1, first_column_index = FALSE)
        if (ec_type == 'ECT'){
          curve = curve
        }
        if (ec_type == 'SECT'){
          curve <- integrate_ec_curve(curve)
        }
        curve_mult_d <- cbind(curve_mult_d,t(curve[,2]))
      }
      data[i,] <- curve_mult_d
    }
    if (ball == TRUE){
      for (j in 1:dim(directions)[1]){
        vertex_function=off$Vertices%*%directions[j,]
        curve <- compute_standardized_ec_curve(off, vertex_function, len-1, first_column_index = FALSE,ball_radius = ball_radius)
        if (ec_type == 'ECT'){
          curve = curve
        }
        if (ec_type == 'SECT'){
          curve <- integrate_ec_curve(curve)
        }
        if (ec_type == 'DECT'){
          curve = differentiate_ec_curve(curve)
        }
        curve_mult_d <- cbind(curve_mult_d,t(curve[,2]))
        # omit the length data, for now
      }
      data[i,] <- curve_mult_d
    }
  }
  return(data)
}
create_comparison_matrix_mult_d=function(directory1,directory2,directions,len, ball = FALSE, ball_radius = 1, ec_type = 'ECT'){
  matrix_1=create_ec_matrix_mult_d(directory1,directions,len,ball_radius = ball_radius, ball = ball, ec_type = ec_type)
  matrix_2=create_ec_matrix_mult_d(directory2,directions,len,ball_radius = ball_radius, ball = ball, ec_type = ec_type)
  matrix_1=cbind(rep(1,dim(matrix_1)[1]),matrix_1)
  matrix_2=cbind(rep(0,dim(matrix_2)[1]),matrix_2)
  final_matrix=rbind(matrix_1,matrix_2)
  return(final_matrix)
}


#### Summary Function for the Real Data ####

real_data_summary = function(dir1,dir2,direction=dir,len=len,class1='Class 1', 
                             class2='Class 2',radius=2,accuracy=FALSE,ec_type = 'ECT', ball = TRUE, ball_radius = 1.5){
  #Generate Matrix of EC curves and labels
  comparison_matrix=create_comparison_matrix_mult_d(dir1,dir2,direction,len,ec_type = ec_type, ball = ball, ball_radius = ball_radius)
  if (accuracy==TRUE){
    accuracy=kfoldcvgp_multiple(iter=100,k=5,comparison_matrix)
  }
  else{
    accuracy=0
  }
  #Feature Selection
  #want_indices=find_rate_variables(comparison_matrix,radius=radius)
  want_indices = 0
  want_indices_rate_total=find_rate_variables_with_other_sampling_methods(comparison_matrix,radius = radius, bandwidth = 0.01,
                                                                          weights = TRUE, type = 'ESS')
  #want_indices_bayesian=find_bayesian_variables(comparison_matrix,param = 0.1,radius = radius)
  #want_indices_lasso=find_lasso_variables(comparison_matrix,radius = radius)
  #want_indices_elastic=find_elastic_variables(comparison_matrix,radius=radius)
  final_list=list(data = comparison_matrix, Rate=want_indices,Rate2=want_indices_rate_total)
  return(final_list)
}
