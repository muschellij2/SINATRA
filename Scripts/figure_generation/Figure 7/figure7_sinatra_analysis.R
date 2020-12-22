setwd('/Users/brucewang/Dropbox (Princeton)/Data + Experiments Tim Sudijono/')
library(R.utils)
#sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')
set.seed(55)
library(sinatra)
library(FNN)
library(rgl)
library(Rvcg)
library(plyr)

#Parameters for the Analysis
pset1 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Microcebus/',dir2 = '~/Documents/doug_new_teeth_by_species/Mirza/',
             num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
                dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset2 =  list(dir1 = '~/Documents/doug_new_teeth_by_species/Microcebus/',dir2 = '~/Documents/doug_new_teeth_by_species/Tarsius/',
              num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
              dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset3 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Microcebus/',dir2 = '~/Documents/doug_new_teeth_by_species/Saimiri/',
             num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset4 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Tarsius/',dir2 = '~/Documents/doug_new_teeth_by_species/Saimiri/',
             num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset5 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Tarsius/',dir2 = '~/Documents/doug_new_teeth_by_species/Mirza/',
             num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset6 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Saimiri/',dir2 = '~/Documents/doug_new_teeth_by_species/Mirza/',
             num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))

cap_radius = 0.25
num_cones = 35
directions_per_cone = 5
len = 75
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  cap_radius, directions_per_cone = directions_per_cone)
rotation_matrix=matrix(c(0.8065218,0.5911149,0.01028626,0,0.5583186,-0.7558216,0.34207344,0,0.1944301,-0.2816325,-0.93961692,0,0,0,0,1),ncol=4,byrow=TRUE)
ball = TRUE
ball_radius = 0.5
ec_type = 'ECT'
#### Start Comparison ####

comp1 = real_data_summary(dir1=pset1$dir1,dir2 = pset1$dir2,direction=pset1$dirs,len = pset1$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
comp2 = real_data_summary(dir1=pset2$dir1,dir2 = pset2$dir2,direction=pset2$dirs,len = pset2$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
comp3 = real_data_summary(dir1=pset3$dir1,dir2 = pset3$dir2,direction=pset3$dirs,len = pset3$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
comp4 = real_data_summary(dir1=pset4$dir1,dir2 = pset4$dir2,direction=pset4$dirs,len = pset4$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
comp5 = real_data_summary(dir1=pset5$dir1,dir2 = pset5$dir2,direction=pset5$dirs,len = pset5$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
comp6 = real_data_summary(dir1=pset6$dir1,dir2 = pset6$dir2,direction=pset6$dirs,len = pset6$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)





color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)

col_pal2=c(color1,color1,color2,color2,color3)
colfunc2 <- colorRampPalette(col_pal2)


#### Helper Function ####
get_vertex_rate = function(dir_name, cuts, pset, comp, ball = TRUE, ball_radius = 0.5){
  mesh_list = list()
  for (k in 1:length(list.files(dir_name,full.names = TRUE))){
    file_name = list.files(dir_name,full.names = TRUE)[k]
    print(paste('On File', file_name))
    mesh = process_off_file_v3(file_name)
    heat = reconstruct_vertices_on_shape(dir = pset$dirs,complex = mesh,rate_vals = comp$Rate2[,2],
                                         len = pset$len,cuts = cuts,cone_size = pset$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)
    mesh_list[[k]] = heat[,2]
  }
  return(mesh_list)
  
}

#### Use Comp 2,4,5 ####

tarsius_microcebus_rates = get_vertex_rate(dir_name = pset2$dir2,cuts = (5*75*35),pset = pset2,comp = comp2)
tarsius_saimiri_rates = get_vertex_rate(dir_name = pset4$dir1,cuts = (5*75*35),pset = pset4,comp = comp4)
tarsius_mirza_rates = get_vertex_rate(dir_name = pset5$dir1,cuts = (5*75*35),pset = pset5,comp = comp5)
save.image('real_data_20190613.Rdata')

#### Write to Directory ####

load('real_data_20190613.Rdata')



#1 : 340
#2 : 237
#3 : 813
#4 : 249
#6 : 291
find_pval_matrices = function(region_sizes,file_names,meshes,indices,num_test_regions,method, tarsius_mirza,tarsius_saimiri,tarsius_microcebus){
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
      k = meshes[i]
      file_name = file_names[k]
      #mkdirs(folder_name)
      print(paste('On File', file_name))
      mesh = process_off_file_v3(file_name)
      index = indices[i]
      
      vertices = rep(0,dim(mesh$Vertices)[1])
      
      tarsius_microcebus_rates = tarsius_microcebus[k][[1]]
      tarsius_saimiri_rates = tarsius_saimiri[k][[1]]
      tarsius_mirza_rates = tarsius_mirza[k][[1]]
      
      #write.csv(tarsius_microcebus_rates, file = paste(folder_name,'/tarsius_microcebus.csv', sep = ''),col.names = FALSE)
      #write.csv(tarsius_saimiri_rates, file = paste(folder_name,'/tarsius_saimiri.csv', sep = ''),col.names = FALSE)
      #write.csv(tarsius_mirza_rates, file = paste(folder_name,'/tarsius_mirza.csv', sep = ''),col.names = FALSE)
      
      mesh = vcgImport(file_name)
      #open3d()
      #plot3d(mesh,color = 'white')
      #rgl.points(matrix(t(mesh$vb[-4,index]),ncol = 3), col = 'blue', size = 10)
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

region_sizes = c(10,50,100,150,200)
num_test_regions = 500
method = 'knn'
meshes = (1:33)
indices = c(261,295,833,224,939,278,252,593,177,357,4845,4215,4619,162,244,4618,4225,31,4705,4524,219,4301,4898,262,4860,69,64,4478,4452,4224,102,179,4176)

knn_pvals = find_pval_matrices(region_sizes = region_sizes,file_names = file_names,meshes = meshes,indices = indices,num_test_regions = num_test_regions,method = method,tarsius_mirza = tarsius_mirza_rates,tarsius_microcebus= tarsius_microcebus_rates,tarsius_saimiri = tarsius_saimiri_rates)

tarsius_microcebus_frame = data.frame(knn_pvals$tarsius_microcebus_pval_matrix)
tarsius_mirza_frame = data.frame(knn_pvals$tarsius_mirza_pval_matrix)
tarsius_saimiri_frame = data.frame(knn_pvals$tarsius_saimiri_pval_matrix)

knn_mean = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=mean)
knn_mean[,1] = region_sizes
colnames(knn_mean) = c('Region Size', 'Tarsius-Microcebus')
knn_mean['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=mean)[2]
knn_mean['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=mean)[2]
knn_mean

knn_median = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=median)
knn_median[,1] = region_sizes
colnames(knn_median) = c('Region Size', 'Tarsius-Microcebus')
knn_median['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=median)[2]
knn_median['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=median)[2]
knn_median

write.csv(tarsius_microcebus_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/knn_method/tarius_microcebus.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_mirza_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/knn_method/tarius_mirza.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_saimiri_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/knn_method/tarius_saimiri.csv',col.names = FALSE,row.names = FALSE)

tarsius_mirza_total = knn_pvals$tarsius_mirza_total_matrix
name = 'tarsius_mirza'
stem = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_rate_sums/knn_method/'
for (k in 1:length(region_sizes)){
  reg_size = region_sizes[k]
  comp_matrix = tarsius_mirza_total[[k]]
  new_name = paste(name,'_roi_size_',reg_size,'.csv',sep = '')
  file_name = paste(stem,new_name,sep = '')
  write.csv(comp_matrix,file = file_name,col.names = FALSE,row.names = FALSE)
  
}
tarsius_microcebus_total = knn_pvals$tarsius_microcebus_total_matrix
name = 'tarsius_microcebus'
stem = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_rate_sums/knn_method/'
for (k in 1:length(region_sizes)){
  reg_size = region_sizes[k]
  comp_matrix = tarsius_microcebus_total[[k]]
  new_name = paste(name,'_roi_size_',reg_size,'.csv',sep = '')
  file_name = paste(stem,new_name,sep = '')
  write.csv(comp_matrix,file = file_name,col.names = FALSE,row.names = FALSE)
  
}

tarsius_saimiri_total = knn_pvals$tarsius_saimiri_total_matrix
name = 'tarsius_saimiri'
stem = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_rate_sums/knn_method/'
for (k in 1:length(region_sizes)){
  reg_size = region_sizes[k]
  comp_matrix = tarsius_saimiri_total[[k]]
  new_name = paste(name,'_roi_size_',reg_size,'.csv',sep = '')
  file_name = paste(stem,new_name,sep = '')
  write.csv(comp_matrix,file = file_name,col.names = FALSE,row.names = FALSE)
  
}


save.image('real_data_20190616.Rdata')
load('real_data_20190616.Rdata')
file_names = list.files(pset2$dir2,full.names = TRUE)
#[21:30]
method = 'area'
region_sizes = c(10,50,100,150,200)
num_test_regions = 500
indices = c(261,295,833,224,939,278,252,593,177,357,4845,4215,4619,162,244,4618,4225,31,4705,4524,219,4301,4898,262,4860,69,64,4478,4452,4224,102,179,4176)
indices = indices
#[21:30]
meshes = (1:33)
#[21:30]
area_pvals = find_pval_matrices(region_sizes = region_sizes,file_names = file_names,meshes = meshes,indices = indices,num_test_regions = num_test_regions,method = method,tarsius_mirza = tarsius_mirza_rates,tarsius_microcebus= tarsius_microcebus_rates,tarsius_saimiri = tarsius_saimiri_rates)

tarsius_microcebus_frame = data.frame(area_pvals$tarsius_microcebus_pval_matrix)
tarsius_mirza_frame = data.frame(area_pvals$tarsius_mirza_pval_matrix)
tarsius_saimiri_frame = data.frame(area_pvals$tarsius_saimiri_pval_matrix)

area_mean = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=mean)
area_mean[,1] = region_sizes
colnames(area_mean) = c('Region Size', 'Tarsius-Microcebus')
area_mean['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=mean)[2]
area_mean['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=mean)[2]

area_median = aggregate(tarsius_microcebus_frame$X1, by=list(tarsius_microcebus_frame$X2), FUN=median)
area_median[,1] = region_sizes
colnames(area_median) = c('Region Size', 'Tarsius-Microcebus')
area_median['Tarsius-Mirza'] = aggregate(tarsius_mirza_frame$X1, by=list(tarsius_mirza_frame$X2), FUN=median)[2]
area_median['Tarsius-Saimiri'] = aggregate(tarsius_saimiri_frame$X1, by=list(tarsius_saimiri_frame$X2), FUN=median)[2]

save.image('real_data_20190616.Rdata')
load('Old Stuff/real_data_20190616.Rdata')

area_mean
area_median
knn_mean
knn_median

write.csv(tarsius_microcebus_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/area_method/tarius_microcebus.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_mirza_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/area_method/tarius_mirza.csv',col.names = FALSE,row.names = FALSE)
write.csv(tarsius_saimiri_frame,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_pval_data/area_method/tarius_saimiri.csv',col.names = FALSE,row.names = FALSE)

write.csv(knn_median,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/knn_median_pvals.csv',row.names = FALSE)
write.csv(area_median,file = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/area_median_pvals.csv',row.names = FALSE)

tarsius_mirza_total = area_pvals$tarsius_mirza_total_matrix
name = 'tarsius_mirza'
stem = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_rate_sums/area_method/'
for (k in 1:length(region_sizes)){
  reg_size = region_sizes[k]
  comp_matrix = tarsius_mirza_total[[k]]
  new_name = paste(name,'_roi_size_',reg_size,'.csv',sep = '')
  file_name = paste(stem,new_name,sep = '')
  write.csv(comp_matrix,file = file_name,col.names = FALSE,row.names = FALSE)
  
}
tarsius_microcebus_total = area_pvals$tarsius_microcebus_total_matrix
name = 'tarsius_microcebus'
stem = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_rate_sums/area_method/'
for (k in 1:length(region_sizes)){
  reg_size = region_sizes[k]
  comp_matrix = tarsius_microcebus_total[[k]]
  new_name = paste(name,'_roi_size_',reg_size,'.csv',sep = '')
  file_name = paste(stem,new_name,sep = '')
  write.csv(comp_matrix,file = file_name,col.names = FALSE,row.names = FALSE)
  
}

tarsius_saimiri_total = area_pvals$tarsius_saimiri_total_matrix
name = 'tarsius_saimiri'
stem = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/real_data_results/raw_rate_sums/area_method/'
for (k in 1:length(region_sizes)){
  reg_size = region_sizes[k]
  comp_matrix = tarsius_saimiri_total[[k]]
  new_name = paste(name,'_roi_size_',reg_size,'.csv',sep = '')
  file_name = paste(stem,new_name,sep = '')
  write.csv(comp_matrix,file = file_name,col.names = FALSE,row.names = FALSE)
  
}
