compute_selected_vertices_cones = function(dir, complex, rate_vals, len, threshold=-1, cone_size, ball_radius){
  if (threshold==-1){
    threshold=1/length(rate_vals)
  }
  if ((dim(dir)[1] %% cone_size) != 0){
    print('Number of Cones not a multiple of directions')
    return(0)
  }
  coned_vertices=list()
  for (j in 1:(dim(dir)[1] / cone_size)){
    cone_dirs=dir[((j-1)*(cone_size)+1):(j*cone_size),]
    cone_rate_vals=rate_vals[(j-1)*(cone_size*len)+1:(j*cone_size*len)]
    coned_vertices[[j]]=summarize_vertices(dir = cone_dirs, complex, rate_vals = cone_rate_vals, len,
                                           reduction_operation = intersect, threshold, cone_size, ball_radius)
  }
  total_selected_vertices=Reduce(union,coned_vertices)
  return(total_selected_vertices)
}


summarize_vertices = function(dir,complex,rate_vals,len,reduction_operation=intersect,threshold,cone_size, ball_radius){
  indices=which(rate_vals>threshold)
  selected_vertices=list()
  
  # Count how many projections are selected for
  for(i in 1:dim(dir)[1]){
    
    vtx_projection <- complex$Vertices[,1:3]%*%dir[i,]
    buckets <- seq(-ball_radius,ball_radius,length.out = len + 1)
    
    # map vertex projection to the feature index
    projection_bucket <- cut(vtx_projection, buckets, labels = FALSE)
    
    # update index to reflect rate values
    projection_bucket <- projection_bucket + (i - 1)*len
    
    selected_vertices[[i]] <- which(projection_bucket %in% indices)
  }
  
  
  final_selected_vertices <- Reduce(reduction_operation,selected_vertices)
  
  return(final_selected_vertices)
}