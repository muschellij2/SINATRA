compute_standardized_ec_curve <- function(complex, vertex_function, curve_length, first_column_index = FALSE, ball_radius){
  
  V = complex$Vertices
  E = complex$Edges
  F = complex$Faces
  
  
  if ( dim(V)[1] != dim(vertex_function)[1] ) {
    print('The size of function should be same as the number of vertices')
    return(0)
  }
  
  # In the case where the vertex_indices are given in the first column. We pick out the birth times of each vertex in every
  # edge and face of the complex.
  if (first_column_index == TRUE){
    
    edge_birth_times = sapply(1:dim(E)[1],function(ind) max(vertex_function[ which(vertex_function[,1]%in%E[ind,]), 2 ] ))
    face_birth_times = sapply(1:dim(F)[1],function(ind) max(vertex_function[ which(vertex_function[,1]%in%F[ind,]), 2 ] ))
  }else{
    
    #Apply a function on indices on this vector? Count the number of edges that are present given
    # the value of the filtration (get max vertex value at each edge, face).
    edge_birth_times = sapply(1:dim(E)[1],function(ind) max(vertex_function[ E[ind,] ]))
    face_birth_times = sapply(1:dim(F)[1],function(ind) max(vertex_function[ F[ind,] ]))
  }
  
  if (first_column_index == FALSE){
    threshold = seq(from=-ball_radius,to=ball_radius,length =curve_length+1)
  } else{
    threshold = seq(from=min(vertex_function[,2])-buffer,to=max(vertex_function[,2])+buffer,length =curve_length+1)
  }
  
  
  ec_curve = matrix(0,curve_length+1,2)
  ec_curve[,1] <- threshold
  
  # count how many of each object is born given the threshold time.
  for ( i in 1:length(threshold) ) {
    v = length(which(vertex_function <= threshold[i]))
    e = length(which(edge_birth_times <= threshold[i]))
    f = length(which(face_birth_times <= threshold[i]))
    ec_curve[i,2] <- v - e + f;
  }
  return(ec_curve)
}

differentiate_ec_curve <- function(ec_curve){
  differences <- diff(ec_curve[,2],lag = 1)
  ec_curve[,2] <- c(ec_curve[1,2], differences)
  
  ec_curve
}

integrate_ec_curve <- function(ec_curve){
  length <- length(ec_curve[,2])
  ec_curve[,2] <- ec_curve[,2] - mean(ec_curve[,2])
  ec_curve[,2] <- cumsum(ec_curve[,2])* (ec_curve[length,1] - ec_curve[1,1])/length
  
  ec_curve
}

# helper function for the sphere data generation function
update_ec_curve = function(curve,ec_type){
  if (ec_type == "SECT"){
    curve <- integrate_ec_curve(curve)
  } else if(ec_type == "DECT"){
    curve <- differentiate_ec_curve(curve)
  } else{
    curve <- curve
  }
  curve
}

