setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
library(R.utils)
sourceDirectory('~/Documents/SINATRA/SINATRA_Pipeline_Branch/')
sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')
#load('Caricature_ROC.Rdata')
#library(svd)
#set.seed(15)
#set.seed(35)
set.seed(15)

#Specifying Directories and the Associated Files


num_causal_region <- 2
num_shared_region <- 1
causal_points <- 10
shared_points <- 10



g = generate_ROC_with_coned_directions(nsim = 25, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                  causal_points = causal_points,shared_points = shared_points, desired_num_cones = 10, eta = 0.1, 
                                  truncated = 100, two_curves = FALSE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                  min_points = 3,directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                  mode = 'sphere', fpr = 0.05, start = 1, cusps = 50,
                                  subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region)



simulation_results <- foreach(i=1:n.simulations, .combine = 'rbind', .noexport = c('GaussKernel')) %:% 
  foreach(j=c(1,5,10,15,20), .combine = 'rbind', .noexport = c('GaussKernel')) %dopar% {
    
    set.seed(5*i+j)
    
    res <- tryCatch( generate_ROC_with_coned_directions(nsim = 25, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                                        causal_points = causal_points,shared_points = shared_points, desired_num_cones = j, eta = 0.1, 
                                                        truncated = -1, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                                        min_points = 3,directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                                        mode = 'sphere', fpr = 0.05, start = 1, cusps = 50,
                                                        subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region),
                     error = function(x) {
                       return(matrix(nrow = 0,ncol = 3))
                     }
    )
    ### Label the results for each trial and directions ###
    rdf <- cbind(res, rep(j, dim(res)[1]) )
    rdf <- data.frame(cbind(rdf,rep(i,dim(rdf)[1])))
    rdf <- plyr::rename(rdf,c("X1" = "FPR","X2" = "TPR","X3" = "Class","X4" = "Index","X5" = "Num_Cones","X6" = "Trial"))
  }
