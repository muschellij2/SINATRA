rm(list=ls())
set.seed(4913, kind = "L'Ecuyer-CMRG")
library(sinatra)
library(FNN)
library(rgl)
library(Rvcg)
library(plyr)
library(pdist)
library(gglasso)
library(numbers)
library(data.table)
library(stringr)
library(ggplot2)
#Parameters for the Analysis

cap_radius = 0.15
num_cones = 5
directions_per_cone = 5
len = 75
num_vertices = 5131

#### Function ####
summarize_list = function(roc_curve_list){
  num_curves = length(roc_curve_list)
  curve = roc_curve_list[[1]]
  
  for (i in 2:num_curves){
    curve = curve+roc_curve_list[[i]]
  }
  
  curve = curve/num_curves
  return(curve)
  
}


scale_and_normalize = function(x){
  x = abs(x)
  x = x/max(x)
  return(x)
}

pval_correct = function(x){
  if (x > 0.1){
    return(0)
  }
  else{
    new_val = -exp(1) * x * log(x)
    return(1/new_val)
  }
}

find_pval_matrices = function(region_sizes,file_paths, file_names,indices,num_test_regions,method,reduce = max, functions_to_summarize = 10, stem = '~/Documents/real_data_comparisons'){
  tarsius_microcebus_total_matrix = list()
  tarsius_saimiri_total_matrix = list()
  tarsius_mirza_total_matrix = list()
  
  
  tarsius_mirza_pval_matrix = matrix(0,nr = length(region_sizes) *  length(file_names),nc = 2)
  tarsius_microcebus_pval_matrix = matrix(0,nr = length(region_sizes) *  length(file_names),nc = 2)
  tarsius_saimiri_pval_matrix = matrix(0,nr = length(region_sizes) *  length(file_names),nc = 2)
  
  for (j in 1:length(region_sizes)){
    region_size = region_sizes[j]
    print(region_size)
    tarsius_microcebus_matrix = matrix(0,nr = length(meshes),nc = (num_test_regions +1))
    tarsius_mirza_matrix = matrix(0,nr = length(meshes),nc = (num_test_regions +1))
    tarsius_saimiri_matrix = matrix(0,nr = length(meshes),nc = (num_test_regions +1))
    
    tarsius_microcebus_pvals = c()
    tarsius_mirza_pvals = c()
    tarsius_saimiri_pvals = c()
    for (i in 1:length(meshes)){
      file_name = file_names[i]
      file_path = file_paths[i]
      #mkdirs(folder_name)
      print(paste('On File', file_name))
      mesh = process_off_file_v3(file_path)
      index = indices[i]
      
      vertices = rep(0,dim(mesh$Vertices)[1])
      mirza_weight_path = paste(stem,'/Tarsius_Mirza/map_',file_names[i],'.csv',sep = '')
      microcebus_weight_path = paste(stem, '/Tarsius_Microcebus/map_',file_names[i],'.csv',sep = '')
      saimiri_weight_path = paste(stem,'/Tarsius_Saimiri/map_',file_names[i],'.csv',sep = '')
      
      g1 = read.csv(mirza_weight_path, header = FALSE)
      g2 = read.csv(microcebus_weight_path, header = FALSE)
      g3 = read.csv(saimiri_weight_path, header = FALSE)
      
      g1 = g1[,1:functions_to_summarize]
      g2 = g3[,1:functions_to_summarize]
      g3 = g2[,1:functions_to_summarize]
      
      g1 = apply(g1, MARGIN = 2, FUN = scale_and_normalize)
      g2 = apply(g2, MARGIN = 2, FUN = scale_and_normalize)
      g3 = apply(g3, MARGIN = 2, FUN = scale_and_normalize)
      
      tarsius_mirza_rates = apply(g1, MARGIN = 1, FUN = reduce)
      tarsius_saimiri_rates = apply(g2, MARGIN = 1, FUN = reduce)
      tarsius_microcebus_rates = apply(g3, MARGIN = 1, FUN = reduce)
      
      
      mesh = vcgImport(file_path)
      
      tarsius_mirza_values = compute_differential_evidence(complex = mesh, mesh_fcn = tarsius_mirza_rates,vertex = index,region_size = region_size,num_test_regions = num_test_regions,method = method)
      tarsius_saimiri_values = compute_differential_evidence(complex = mesh, mesh_fcn = tarsius_saimiri_rates,vertex = index,region_size = region_size,num_test_regions = num_test_regions,method = method)
      tarsius_microcebus_values = compute_differential_evidence(complex = mesh, mesh_fcn = tarsius_microcebus_rates,vertex = index,region_size = region_size,num_test_regions = num_test_regions,method = method)
      
      tarsius_microcebus_matrix[i,] = tarsius_microcebus_values
      tarsius_mirza_matrix[i,]      = tarsius_mirza_values
      tarsius_saimiri_matrix[i,]    = tarsius_saimiri_values
      
      
      tarsius_microcebus_pval = length(which(tarsius_microcebus_values[1] < tarsius_microcebus_values[-1]))/(num_test_regions+1)
      tarsius_mirza_pval = length(which(tarsius_mirza_values[1] < tarsius_mirza_values[-1]))/(num_test_regions+1)
      tarsius_saimiri_pval = length(which(tarsius_saimiri_values[1] < tarsius_saimiri_values[-1]))/(num_test_regions+1)
      
      tarsius_microcebus_pvals = c(tarsius_microcebus_pvals,tarsius_microcebus_pval)
      tarsius_mirza_pvals = c(tarsius_mirza_pvals,tarsius_mirza_pval)
      tarsius_saimiri_pvals = c(tarsius_saimiri_pvals,tarsius_saimiri_pval)
      
    }
    min_val = ((j-1) * length(file_names)) + 1
    max_val = ((j) * length(file_names)) 
    tarsius_mirza_pval_matrix[min_val:max_val,1] = tarsius_mirza_pvals
    tarsius_mirza_pval_matrix[min_val:max_val,2] = region_size
    tarsius_microcebus_pval_matrix[min_val:max_val,1] = tarsius_microcebus_pvals
    tarsius_microcebus_pval_matrix[min_val:max_val,2] = region_size
    tarsius_saimiri_pval_matrix[min_val:max_val,1] = tarsius_saimiri_pvals
    tarsius_saimiri_pval_matrix[min_val:max_val,2] = region_size
    
    tarsius_microcebus_total_matrix[[j]] = tarsius_microcebus_matrix
    tarsius_saimiri_total_matrix[[j]] = tarsius_saimiri_matrix
    tarsius_mirza_total_matrix[[j]] = tarsius_mirza_matrix
  }
  
  final_list = list(tarsius_mirza_pval_matrix = tarsius_mirza_pval_matrix, tarsius_saimiri_pval_matrix = tarsius_saimiri_pval_matrix,
                    tarsius_microcebus_pval_matrix = tarsius_microcebus_pval_matrix, tarsius_microcebus_total_matrix = tarsius_microcebus_total_matrix,
                    tarsius_saimiri_total_matrix = tarsius_saimiri_total_matrix, tarsius_mirza_total_matrix = tarsius_mirza_total_matrix)
  
  return(final_list)
  
  
}

data_path = '~/Dropbox (Princeton)/Data + Experiments Tim Sudijono/LimitShapeCodeDemo/cleaned_real_data/Tarsius/'
file_paths = list.files(path = data_path,full.names = TRUE)
file_names = list.files(path = data_path,full.names = FALSE)

region_sizes = c(10,50,100,150,200)
num_test_regions = 500
method = 'knn'
meshes = (1:33)
indices = c(261,295,833,224,939,278,252,593,177,357,4845,4215,4619,162,244,4618,4225,31,4705,4524,219,4301,4898,262,4860,69,64,4478,4452,4224,102,179,4176)
functions_to_summarize = 15
reduce = mean

knn_pvals_mean = find_pval_matrices(region_sizes = region_sizes,file_names = file_names,file_paths = file_paths, 
                                  indices = indices,num_test_regions = num_test_regions,method = method, functions_to_summarize = functions_to_summarize,reduce = reduce)
save.image('~/Documents/real_data_comparisons/2020-0626.RData')
reduce = max
knn_pvals_max = find_pval_matrices(region_sizes = region_sizes,file_names = file_names,file_paths = file_paths, 
                                  indices = indices,num_test_regions = num_test_regions,method = method, functions_to_summarize = functions_to_summarize,reduce = reduce)
save.image('~/Documents/real_data_comparisons/2020-0626.RData')

tarsius_microcebus_frame = data.frame(knn_pvals_mean$tarsius_microcebus_pval_matrix)
tarsius_mirza_frame = data.frame(knn_pvals_mean$tarsius_mirza_pval_matrix)
tarsius_saimiri_frame = data.frame(knn_pvals_mean$tarsius_saimiri_pval_matrix)

write.csv(tarsius_microcebus_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/knn_method/tarius_microcebus_limit_shapes_mean.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_mirza_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/knn_method/tarius_mirza_limit_shapes_mean.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_saimiri_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/knn_method/tarius_saimiri_limit_shapes_mean.csv',col.names = FALSE,row.names = FALSE)

knn_mean_mean = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=mean)
knn_mean_mean[,1] = region_sizes
colnames(knn_mean_mean) = c('Region Size', 'Tarsius-Microcebus')
knn_mean_mean['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=mean)[2]
knn_mean_mean['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=mean)[2]
knn_mean_mean

knn_mean_median = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=median)
knn_mean_median[,1] = region_sizes
colnames(knn_mean_median) = c('Region Size', 'Tarsius-Microcebus')
knn_mean_median['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=median)[2]
knn_mean_median['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=median)[2]
knn_mean_median


tarsius_microcebus_frame = data.frame(knn_pvals_max$tarsius_microcebus_pval_matrix)
tarsius_mirza_frame = data.frame(knn_pvals_max$tarsius_mirza_pval_matrix)
tarsius_saimiri_frame = data.frame(knn_pvals_max$tarsius_saimiri_pval_matrix)

write.csv(tarsius_microcebus_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/knn_method/tarius_microcebus_limit_shapes_max.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_mirza_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/knn_method/tarius_mirza_limit_shapes_max.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_saimiri_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/knn_method/tarius_saimiri_limit_shapes_max.csv',col.names = FALSE,row.names = FALSE)

knn_max_mean = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=mean)
knn_max_mean[,1] = region_sizes
colnames(knn_max_mean) = c('Region Size', 'Tarsius-Microcebus')
knn_max_mean['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=mean)[2]
knn_max_mean['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=mean)[2]
knn_max_mean

knn_max_median = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=median)
knn_max_median[,1] = region_sizes
colnames(knn_max_median) = c('Region Size', 'Tarsius-Microcebus')
knn_max_median['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=median)[2]
knn_max_median['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=median)[2]
knn_max_median

knn_mean_median_corrected = apply(X = knn_mean_median,MARGIN = c(1,2),FUN = pval_correct)
knn_mean_median_corrected[,1] = knn_mean_median[,1]
save.image('~/Documents/real_data_comparisons/2020-0626.RData')
load('~/Documents/real_data_comparisons/2020-0626.RData')

#### Area based ####

region_sizes = c(10,50,100,150,200)
num_test_regions = 500
method = 'area'
meshes = (1:33)
indices = c(261,295,833,224,939,278,252,593,177,357,4845,4215,4619,162,244,4618,4225,31,4705,4524,219,4301,4898,262,4860,69,64,4478,4452,4224,102,179,4176)
functions_to_summarize = 15
reduce = mean

area_pvals_mean = find_pval_matrices(region_sizes = region_sizes,file_names = file_names,file_paths = file_paths, 
                                    indices = indices,num_test_regions = num_test_regions,method = method, functions_to_summarize = functions_to_summarize,reduce = reduce)
save.image('~/Documents/real_data_comparisons/2020-0626-area.RData')
load('~/Documents/real_data_comparisons/2020-0626-area.RData')
reduce = max
area_pvals_max = find_pval_matrices(region_sizes = region_sizes,file_names = file_names,file_paths = file_paths, 
                                   indices = indices,num_test_regions = num_test_regions,method = method, functions_to_summarize = functions_to_summarize,reduce = reduce)
save.image('~/Documents/real_data_comparisons/2020-0626-area.RData')

tarsius_microcebus_frame = data.frame(area_pvals_mean$tarsius_microcebus_pval_matrix)
tarsius_mirza_frame = data.frame(area_pvals_mean$tarsius_mirza_pval_matrix)
tarsius_saimiri_frame = data.frame(area_pvals_mean$tarsius_saimiri_pval_matrix)

write.csv(tarsius_microcebus_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/area_method/tarius_microcebus_limit_shapes_mean.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_mirza_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/area_method/tarius_mirza_limit_shapes_mean.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_saimiri_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/area_method/tarius_saimiri_limit_shapes_mean.csv',col.names = FALSE,row.names = FALSE)


area_mean_mean = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=mean)
area_mean_mean[,1] = region_sizes
colnames(area_mean_mean) = c('Region Size', 'Tarsius-Microcebus')
area_mean_mean['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=mean)[2]
area_mean_mean['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=mean)[2]
area_mean_mean

area_mean_median = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=median)
area_mean_median[,1] = region_sizes
colnames(area_mean_median) = c('Region Size', 'Tarsius-Microcebus')
area_mean_median['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=median)[2]
area_mean_median['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=median)[2]
area_mean_median


tarsius_microcebus_frame = data.frame(area_pvals_max$tarsius_microcebus_pval_matrix)
tarsius_mirza_frame = data.frame(area_pvals_max$tarsius_mirza_pval_matrix)
tarsius_saimiri_frame = data.frame(area_pvals_max$tarsius_saimiri_pval_matrix)

write.csv(tarsius_microcebus_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/area_method/tarius_microcebus_limit_shapes_max.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_mirza_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/area_method/tarius_mirza_limit_shapes_max.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_saimiri_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/area_method/tarius_saimiri_limit_shapes_max.csv',col.names = FALSE,row.names = FALSE)

area_max_mean = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=mean)
area_max_mean[,1] = region_sizes
colnames(area_max_mean) = c('Region Size', 'Tarsius-Microcebus')
area_max_mean['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=mean)[2]
area_max_mean['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=mean)[2]
area_max_mean

area_max_median = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=median)
area_max_median[,1] = region_sizes
colnames(area_max_median) = c('Region Size', 'Tarsius-Microcebus')
area_max_median['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=median)[2]
area_max_median['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=median)[2]
area_max_median

area_mean_median_corrected = apply(X = area_mean_median,MARGIN = c(1,2),FUN = pval_correct)
area_mean_median_corrected[,1] = area_mean_median[,1]

area_mean_median_corrected
#### Write ####
save.image('~/Documents/real_data_comparisons/2020-0626-area.RData')
load('~/Documents/real_data_comparisons/2020-0626-area.RData')


area_max_mean
area_max_median
area_mean_mean
area_mean_median
knn_max_mean
knn_max_median
knn_mean_mean
knn_mean_median

write.csv(area_mean_median,'~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/limit_shapes_area_mean_median_pvals.csv')
write.csv(knn_mean_median,'~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/limit_shapes_knn_mean_median_pvals.csv')

#### Misspecified ####
data_path = '~/Dropbox (Princeton)/Data + Experiments Tim Sudijono/LimitShapeCodeDemo/cleaned_real_data/Tarsius/'
file_paths = list.files(path = data_path,full.names = TRUE)
file_names = list.files(path = data_path,full.names = FALSE)
stem = '~/Documents/real_data_comparisons_misspecified'

region_sizes = c(10,50,100,150,200)
num_test_regions = 500
method = 'knn'
meshes = (1:33)
indices = c(261,295,833,224,939,278,252,593,177,357,4845,4215,4619,162,244,4618,4225,31,4705,4524,219,4301,4898,262,4860,69,64,4478,4452,4224,102,179,4176)
functions_to_summarize = 15
reduce = mean

misspecified_knn_pvals_mean = find_pval_matrices(region_sizes = region_sizes,file_names = file_names,file_paths = file_paths, 
                                                 indices = indices,num_test_regions = num_test_regions,method = method, functions_to_summarize = functions_to_summarize,reduce = reduce, stem = stem)
save.image('~/Documents/real_data_comparisons_misspecified/2020-0626.RData')
reduce = max
misspecified_knn_pvals_max = find_pval_matrices(region_sizes = region_sizes,file_names = file_names,file_paths = file_paths, 
                                                indices = indices,num_test_regions = num_test_regions,method = method, functions_to_summarize = functions_to_summarize,reduce = reduce, stem = stem)
save.image('~/Documents/real_data_comparisons_misspecified/2020-0626.RData')

tarsius_microcebus_frame = data.frame(misspecified_knn_pvals_mean$tarsius_microcebus_pval_matrix)
tarsius_mirza_frame = data.frame(misspecified_knn_pvals_mean$tarsius_mirza_pval_matrix)
tarsius_saimiri_frame = data.frame(misspecified_knn_pvals_mean$tarsius_saimiri_pval_matrix)

write.csv(tarsius_microcebus_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/misspecified_knn_method/tarius_microcebus_limit_shapes_mean.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_mirza_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/misspecified_knn_method/tarius_mirza_limit_shapes_mean.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_saimiri_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/misspecified_knn_method/tarius_saimiri_limit_shapes_mean.csv',col.names = FALSE,row.names = FALSE)

misspecified_knn_mean_mean = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=mean)
misspecified_knn_mean_mean[,1] = region_sizes
colnames(misspecified_knn_mean_mean) = c('Region Size', 'Tarsius-Microcebus')
misspecified_knn_mean_mean['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=mean)[2]
misspecified_knn_mean_mean['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=mean)[2]
misspecified_knn_mean_mean

misspecified_knn_mean_median = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=median)
misspecified_knn_mean_median[,1] = region_sizes
colnames(misspecified_knn_mean_median) = c('Region Size', 'Tarsius-Microcebus')
misspecified_knn_mean_median['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=median)[2]
misspecified_knn_mean_median['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=median)[2]
misspecified_knn_mean_median


tarsius_microcebus_frame = data.frame(misspecified_knn_pvals_max$tarsius_microcebus_pval_matrix)
tarsius_mirza_frame = data.frame(misspecified_knn_pvals_max$tarsius_mirza_pval_matrix)
tarsius_saimiri_frame = data.frame(misspecified_knn_pvals_max$tarsius_saimiri_pval_matrix)

write.csv(tarsius_microcebus_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/misspecified_knn_method/tarius_microcebus_limit_shapes_max.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_mirza_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/misspecified_knn_method/tarius_mirza_limit_shapes_max.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_saimiri_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/misspecified_knn_method/tarius_saimiri_limit_shapes_max.csv',col.names = FALSE,row.names = FALSE)

misspecified_knn_max_mean = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=mean)
misspecified_knn_max_mean[,1] = region_sizes
colnames(misspecified_knn_max_mean) = c('Region Size', 'Tarsius-Microcebus')
misspecified_knn_max_mean['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=mean)[2]
misspecified_knn_max_mean['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=mean)[2]
misspecified_knn_max_mean

misspecified_knn_max_median = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=median)
misspecified_knn_max_median[,1] = region_sizes
colnames(misspecified_knn_max_median) = c('Region Size', 'Tarsius-Microcebus')
misspecified_knn_max_median['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=median)[2]
misspecified_knn_max_median['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=median)[2]
misspecified_knn_max_median

save.image('~/Documents/real_data_comparisons_misspecified/2020-0626.RData')
load('~/Documents/real_data_comparisons_misspecified/2020-0626.RData')

#### Area based ####

region_sizes = c(10,50,100,150,200)
num_test_regions = 500
method = 'area'
meshes = (1:33)
indices = c(261,295,833,224,939,278,252,593,177,357,4845,4215,4619,162,244,4618,4225,31,4705,4524,219,4301,4898,262,4860,69,64,4478,4452,4224,102,179,4176)
functions_to_summarize = 15
reduce = mean

misspecified_area_pvals_mean = find_pval_matrices(region_sizes = region_sizes,file_names = file_names,file_paths = file_paths, 
                                                  indices = indices,num_test_regions = num_test_regions,method = method, functions_to_summarize = functions_to_summarize,reduce = reduce)
save.image('~/Documents/real_data_comparisons_misspecified/2020-0626-area.RData')
load('~/Documents/real_data_comparisons_misspecified/2020-0626-area.RData')
reduce = max
misspecified_area_pvals_max = find_pval_matrices(region_sizes = region_sizes,file_names = file_names,file_paths = file_paths, 
                                                 indices = indices,num_test_regions = num_test_regions,method = method, functions_to_summarize = functions_to_summarize,reduce = reduce)
save.image('~/Documents/real_data_comparisons_misspecified/2020-0626-area.RData')

tarsius_microcebus_frame = data.frame(misspecified_knn_pvals_mean$tarsius_microcebus_pval_matrix)
tarsius_mirza_frame = data.frame(misspecified_knn_pvals_mean$tarsius_mirza_pval_matrix)
tarsius_saimiri_frame = data.frame(misspecified_knn_pvals_mean$tarsius_saimiri_pval_matrix)

write.csv(tarsius_microcebus_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/misspecified_area_method/tarius_microcebus_limit_shapes_mean.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_mirza_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/misspecified_area_method/tarius_mirza_limit_shapes_mean.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_saimiri_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/misspecified_area_method/tarius_saimiri_limit_shapes_mean.csv',col.names = FALSE,row.names = FALSE)


misspecified_area_mean_mean = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=mean)
misspecified_area_mean_mean[,1] = region_sizes
colnames(misspecified_area_mean_mean) = c('Region Size', 'Tarsius-Microcebus')
misspecified_area_mean_mean['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=mean)[2]
misspecified_area_mean_mean['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=mean)[2]
misspecified_area_mean_mean

misspecified_area_mean_median = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=median)
misspecified_area_mean_median[,1] = region_sizes
colnames(misspecified_area_mean_median) = c('Region Size', 'Tarsius-Microcebus')
misspecified_area_mean_median['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=median)[2]
misspecified_area_mean_median['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=median)[2]
misspecified_area_mean_median


tarsius_microcebus_frame = data.frame(misspecified_area_pvals_max$tarsius_microcebus_pval_matrix)
tarsius_mirza_frame = data.frame(misspecified_area_pvals_max$tarsius_mirza_pval_matrix)
tarsius_saimiri_frame = data.frame(misspecified_area_pvals_max$tarsius_saimiri_pval_matrix)

write.csv(tarsius_microcebus_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/misspecified_area_method/tarius_microcebus_limit_shapes_max.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_mirza_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/misspecified_area_method/tarius_mirza_limit_shapes_max.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_saimiri_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/misspecified_area_method/tarius_saimiri_limit_shapes_max.csv',col.names = FALSE,row.names = FALSE)

misspecified_area_max_mean = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=mean)
misspecified_area_max_mean[,1] = region_sizes
colnames(misspecified_area_max_mean) = c('Region Size', 'Tarsius-Microcebus')
misspecified_area_max_mean['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=mean)[2]
misspecified_area_max_mean['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=mean)[2]
misspecified_area_max_mean

misspecified_area_max_median = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=median)
misspecified_area_max_median[,1] = region_sizes
colnames(misspecified_area_max_median) = c('Region Size', 'Tarsius-Microcebus')
misspecified_area_max_median['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=median)[2]
misspecified_area_max_median['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=median)[2]
misspecified_area_max_median


#### Write ####
save.image('~/Documents/real_data_comparisons_misspecified/2020-0626-area.RData')
load('~/Documents/real_data_comparisons_misspecified/2020-0626-area.RData')


misspecified_area_max_mean
misspecified_area_max_median
misspecified_area_mean_mean
misspecified_area_mean_median
misspecified_knn_max_mean
misspecified_knn_max_median
misspecified_knn_mean_mean
misspecified_knn_mean_median

write.csv(misspecified_area_mean_median,'~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/limit_shapes_misspecified_area_mean_median_pvals.csv')
write.csv(misspecified_knn_mean_median,'~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/limit_shapes_misspecified_knn_mean_median_pvals.csv')


mesh = vcgImport(data_files[i])
cols = rep('white',dim(mesh$vb)[2])
cols[261] = 'red'
plot3d(mesh, col = cols)
weight_path = paste('~/Documents/real_data_comparisons/Tarsius_',comparison,'/map_',file_names[i],'.csv',sep = '')
g = read.csv(weight_path, header = FALSE)
g = g[,1:functions_to_summarize]
g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
rate_values = apply(g, MARGIN = 1, FUN = max)


color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)
cuts = cut(rate_values,20,labels = FALSE)

heat_colors=colfunc(max(cuts) - min(cuts))[cuts - min(cuts)]

plot3d(mesh, col = heat_colors)

#### See what's going on ####
path = '~/Dropbox (Princeton)/Data + Experiments Tim Sudijono/LimitShapeCodeDemo/cleaned_real_data/Mirza/'
teeth = list.files(path, full.names = TRUE)
for (tooth in teeth){
  mesh = vcgImport(tooth)
  print(dim(mesh$vb))
}
