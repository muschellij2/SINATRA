set.seed(4913, kind = "L'Ecuyer-CMRG")
library(sinatra)
library(FNN)
library(rgl)
library(Rvcg)
library(plyr)
library(pdist)
library(numbers)
library(data.table)
library(tidyverse)

setwd('~/Documents')

#Prelim Function
scale_and_normalize = function(x){
  x = abs(x)
  x = x/max(x)
  return(x)
}


#Parameters for the Analysis
path = '~/Documents/spheres2020-04-23'
dirs = list.dirs(path, full.names = TRUE, FALSE)
roc_curves1.1 = list()
roc_curves1.2 = list()
num_vertices = 642
functions_to_summarize = 10
for (i in 1:length(dirs)){
  dir = dirs[i]
  print(dir)
  total_rate_roc1 = matrix(0, nrow = num_vertices,ncol = 2)
  total_rate_roc2 = matrix(0, nrow = num_vertices,ncol = 2)
  dir1 = paste(dir,'/v1/', sep = '')
  dir2 = paste(dir,'/v2/', sep = '')
  causal_points1_vec = read.csv(paste(dir,'/gp1.csv',sep = ''))
  causal_points2_vec = read.csv(paste(dir,'/gp2.csv',sep = ''))
  
  class_1_true_vertices = which(causal_points1_vec == 1)
  class_2_true_vertices = which(causal_points2_vec == 1)
  
  class_1_false_vertices = setdiff(1:num_vertices,class_1_true_vertices)
  class_2_false_vertices = setdiff(1:num_vertices,class_2_true_vertices)
  
  limit_shapes_output_files1 = list.files(path = dir1, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs1 = limit_shapes_output_files1[str_detect(limit_shapes_output_files1,'csv')]
  
  limit_shapes_output_files2 = list.files(path = dir2, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs2 = limit_shapes_output_files2[str_detect(limit_shapes_output_files2,'csv')]
  
  #ROC curve for class 1
  
  counter = 0 
  for (file in limit_shapes_csvs1){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){
      
      #sink("/dev/null")
      rate_positive_vertices = which(rate_values >= threshold)
      
      rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)
      rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                    class_1_true_vertices,class_1_false_vertices))
      
      true_vertices = class_1_true_vertices
      false_vertices = class_1_false_vertices
    }
    total_rate_roc1 = total_rate_roc1 + rate_ROC
      
  }
  total_rate_roc1 = (total_rate_roc1 / counter)
  #ROC Curve for Class 2
  counter = 0 
  for (file in limit_shapes_csvs2){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){
      
      #sink("/dev/null")
      rate_positive_vertices = which(rate_values >= threshold)
      
      rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)
      rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                    class_2_true_vertices,class_2_false_vertices))
      
      true_vertices = class_2_true_vertices
      false_vertices = class_2_false_vertices
    }
    total_rate_roc2 = total_rate_roc2 + rate_ROC
      
  }
  total_rate_roc2 = (total_rate_roc2 / counter)
  
  roc_curves1.1[[i]] = total_rate_roc1
  roc_curves1.2[[i]] = total_rate_roc2
}
#roc_curves1.1 = roc_curves1
#roc_curves1.2 = roc_curves2

total_roc1 = matrix(0, nrow = dim(roc_curves1.1[[1]]), ncol = dim(roc_curves1.2[[1]])[2])
for (j in 1:length(roc_curves1.1)){
  total_roc1 = total_roc1 + roc_curves1.1[[j]]
}
total_roc1 = total_roc1/length(roc_curves1.1)
total_roc2 = matrix(0, nrow = dim(roc_curves1.1[[1]]), ncol = dim(roc_curves1.2[[1]])[2])
for (j in 1:length(roc_curves1.2)){
  total_roc2 = total_roc2 + roc_curves1.2[[j]]
}
total_roc2 = total_roc2/length(roc_curves1.2)
total_roc1 =  cbind(total_roc1,rep('Limit Shapes Class 1', dim(total_roc1)[1]))
total_roc2 =  cbind(total_roc2,rep('Limit Shapes Class 2', dim(total_roc2)[1]))
roc_curve_frame_limit = as.data.frame(rbind(total_roc1, total_roc2))
roc_curve_frame_limit$V1 = as.numeric(as.character(roc_curve_frame_limit$V1))
roc_curve_frame_limit$V2 = as.numeric(as.character(roc_curve_frame_limit$V2))

ggplot(roc_curve_frame_limit, aes(x = V1,y = V2,group = V3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(V3) )) +
  geom_line(stat = "identity") +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Cones") +
  ggtitle(sprintf("2 Causal Regions, 1 Shared Regions, Size 10")) +
  geom_abline(intercept = 0, slope = 1) + 
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

save.image('~/Documents/spheres2020-04-23/roc_curves_10_function.Rdata')
load('~/Dropbox (Princeton)/SINATRA_Data/spheres2020-04-23/roc_curves_10_function.Rdata')

## Scrambled
path = '~/Documents/spheres2020-04-23'
dirs = list.dirs(path, full.names = TRUE, FALSE)
roc_curves1.1 = list()
roc_curves1.2 = list()
num_vertices = 642
functions_to_summarize = 10
for (i in 1:length(dirs)){
  dir = dirs[i]
  print(dir)
  total_rate_roc1 = matrix(0, nrow = num_vertices,ncol = 2)
  total_rate_roc2 = matrix(0, nrow = num_vertices,ncol = 2)
  dir1 = paste(dir,'/v1/', sep = '')
  dir2 = paste(dir,'/v2/', sep = '')
  causal_points1_vec = read.csv(paste(dir,'/gp1.csv',sep = ''))
  causal_points2_vec = read.csv(paste(dir,'/gp2.csv',sep = ''))
  
  class_1_true_vertices = which(causal_points1_vec == 1)
  class_2_true_vertices = which(causal_points2_vec == 1)
  
  class_1_false_vertices = setdiff(1:num_vertices,class_1_true_vertices)
  class_2_false_vertices = setdiff(1:num_vertices,class_2_true_vertices)
  
  limit_shapes_output_files1 = list.files(path = dir1, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs1 = limit_shapes_output_files1[str_detect(limit_shapes_output_files1,'csv')]
  limit_shapes_csvs1 = limit_shapes_csvs1[str_detect(limit_shapes_csvs1,'scrambled')]
  
  limit_shapes_output_files2 = list.files(path = dir2, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs2 = limit_shapes_output_files2[str_detect(limit_shapes_output_files2,'csv')]
  limit_shapes_csvs2 = limit_shapes_csvs2[str_detect(limit_shapes_csvs2,'scrambled')]
  
  #ROC curve for class 1
  
  counter = 0 
  for (file in limit_shapes_csvs1){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){
      
      #sink("/dev/null")
      rate_positive_vertices = which(rate_values >= threshold)
      
      rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)
      rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                    class_1_true_vertices,class_1_false_vertices))
      
      true_vertices = class_1_true_vertices
      false_vertices = class_1_false_vertices
    }
    total_rate_roc1 = total_rate_roc1 + rate_ROC
    
  }
  total_rate_roc1 = (total_rate_roc1 / counter)
  #ROC Curve for Class 2
  counter = 0 
  for (file in limit_shapes_csvs2){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){
      
      #sink("/dev/null")
      rate_positive_vertices = which(rate_values >= threshold)
      
      rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)
      rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                    class_2_true_vertices,class_2_false_vertices))
      
      true_vertices = class_2_true_vertices
      false_vertices = class_2_false_vertices
    }
    total_rate_roc2 = total_rate_roc2 + rate_ROC
    
  }
  total_rate_roc2 = (total_rate_roc2 / counter)
  
  roc_curves1.1[[i]] = total_rate_roc1
  roc_curves1.2[[i]] = total_rate_roc2
}
#roc_curves1.1 = roc_curves1
#roc_curves1.2 = roc_curves2

total_roc1 = matrix(0, nrow = dim(roc_curves1.1[[1]]), ncol = dim(roc_curves1.2[[1]])[2])
for (j in 1:length(roc_curves1.1)){
  total_roc1 = total_roc1 + roc_curves1.1[[j]]
}
total_roc1 = total_roc1/length(roc_curves1.1)
total_roc2 = matrix(0, nrow = dim(roc_curves1.1[[1]]), ncol = dim(roc_curves1.2[[1]])[2])
for (j in 1:length(roc_curves1.2)){
  total_roc2 = total_roc2 + roc_curves1.2[[j]]
}
total_roc2 = total_roc2/length(roc_curves1.2)
total_roc1 =  cbind(total_roc1,rep('Limit Shapes Class 1', dim(total_roc1)[1]))
total_roc2 =  cbind(total_roc2,rep('Limit Shapes Class 2', dim(total_roc2)[1]))
roc_curve_frame_limit_scrambled = as.data.frame(rbind(total_roc1, total_roc2))
roc_curve_frame_limit_scrambled$V1 = as.numeric(as.character(roc_curve_frame_limit_scrambled$V1))
roc_curve_frame_limit_scrambled$V2 = as.numeric(as.character(roc_curve_frame_limit_scrambled$V2))

ggplot(roc_curve_frame_limit_scrambled, aes(x = V1,y = V2,group = V3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(V3) )) +
  geom_line(stat = "identity") +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Cones") +
  ggtitle(sprintf("2 Causal Regions, 1 Shared Regions, Size 10")) +
  geom_abline(intercept = 0, slope = 1) + 
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
save.image('~/Documents/spheres2020-04-23/roc_curves_10_function_scrambled.Rdata')
load('~/Dropbox (Princeton)/SINATRA_Data/spheres2020-04-23/roc_curves_10_function_scrambled.Rdata')

#### Sanity check ####
rate_values = log(rate_values)
shape = vcgImport(file = '/Users/brucewang/Documents/spheres2020-04-23/v1/v1//sphere_9.off')

mfrow3d(1,2)
cols = rep('white',dim(shape$vb)[2])
cols[causal_points1] = 'red'

color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)
cuts = cut(rate_values,20,labels = FALSE)

heat_colors=colfunc(max(cuts) - min(cuts))[cuts - min(cuts)]

plot3d(shape, col = heat_colors)
plot3d(shape, col = cols)


#### Create some data ####
#g1 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
#                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
#                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
#                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
#                                   mode = 'sphere_baseline',num_cusps = cusps,
#                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 0.5,reduce = max, write = FALSE)
#

#write.csv(g1,'~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_elastic_net_max.csv')
g1 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_elastic_net_max2.csv')[,-1][,-4]
roc_curve_frame1.1 = data.frame(g1)
library(ggplot2)
roc_curve_frame1.1 = roc_curve_frame1.1[as.numeric(as.character(roc_curve_frame1.1[,3])) == 1,]
roc_curve_frame1.1[,3] = rep('Baseline (EN Max)',dim(roc_curve_frame1.1)[1])
#g1.2 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
#                                                 causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
#                                                 truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
#                                                 min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
#                                                 mode = 'sphere_baseline',num_cusps = cusps,
#                                                 subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 0.5,reduce = mean)
#
#
#write.csv(g1.2,'~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_elastic_net_mean.csv')
g1.2 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_elastic_net_mean2.csv')[,-1][,-4]
roc_curve_frame1.2 = data.frame(g1.2)
library(ggplot2)
roc_curve_frame1.2 = roc_curve_frame1.2[roc_curve_frame1.2[,3] == 1,]
roc_curve_frame1.2[,3] = rep('Baseline (EN Mean)',dim(roc_curve_frame1.2)[1])
#g1.3 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
#                                                 causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
#                                                 truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
#                                                 min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
#                                                 mode = 'sphere_baseline',num_cusps = cusps,
#                                                 subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 1,reduce = max)
#
#
#write.csv(g1.3,'~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_ridge_max.csv')
g1.3 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_ridge_max2.csv')[,-1][,-4]
roc_curve_frame1.3 = data.frame(g1.3)
library(ggplot2)
roc_curve_frame1.3 = roc_curve_frame1.3[roc_curve_frame1.3[,3] == 1,]
roc_curve_frame1.3[,3] = rep('Baseline (RR Max)',dim(roc_curve_frame1.3)[1])
#g1.4 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
#                                                 causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
#                                                 truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
#                                                 min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
#                                                 mode = 'sphere_baseline',num_cusps = cusps,
#                                                 subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 1,reduce = mean)
#

#write.csv(g1.4,'~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_ridge_mean.csv')
g1.4 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_ridge_mean2.csv')[,-1][,-4]
roc_curve_frame1.4 = data.frame(g1.4)
library(ggplot2)
roc_curve_frame1.4 = roc_curve_frame1.4[roc_curve_frame1.4[,3] == 1,]
roc_curve_frame1.4[,3] = rep('Baseline (RR Mean)',dim(roc_curve_frame1.4)[1])
roc_curve_frame = rbind(roc_curve_frame1.1,roc_curve_frame1.2,roc_curve_frame1.3,roc_curve_frame1.4)
library(ggplot2)
load('~/Documents/SINATRA/Simulations/Results/new/df_ROC_causal1_shared2_10.RData')
rdfmeans = rdfmeans[rdfmeans$Class == 1,]
rdfmeans = rdfmeans[,-2]
rdfmeans = rdfmeans[,-2]
rdfmeans = rdfmeans[rdfmeans$Num_Cones==60,]
rdfmeans = rdfmeans[c('FPR','TPR')]
rdfmeans[,3] = rep('SINATRA',dim(rdfmeans)[1])
names(rdfmeans) = names(roc_curve_frame)
roc_curve_frame = rbind(roc_curve_frame,rdfmeans)
roc_curve_frame_limit2 = roc_curve_frame_limit[roc_curve_frame_limit[,3] == 'Limit Shapes Class 1',]
roc_curve_frame_limit2[,3] = 'Limit Shapes'
roc_curve_frame_limit2 = rbind(c(0,0,'Limit Shapes'),roc_curve_frame_limit2 )
roc_curve_frame_limit_scrambled2 = roc_curve_frame_limit_scrambled[roc_curve_frame_limit_scrambled[,3] == 'Limit Shapes Class 1',]
roc_curve_frame_limit_scrambled2[,3] = 'Limit Shapes (Misspecified)'
roc_curve_frame_limit_scrambled2 = rbind(c(0,0,'Limit Shapes (Misspecified)'),roc_curve_frame_limit_scrambled2 )
roc_curve_frame = rbind(roc_curve_frame,roc_curve_frame_limit2)
roc_curve_frame = rbind(roc_curve_frame,roc_curve_frame_limit_scrambled2)
roc_curve_frame$V1 = as.numeric(as.character(roc_curve_frame$V1))
roc_curve_frame$V2 = as.numeric(as.character(roc_curve_frame$V2))
ggplot(roc_curve_frame, aes(x = V1,y = V2,group = V3)) + 
  geom_line(data = subset(roc_curve_frame, V3 == "SINATRA"), aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 2) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "Limit"), aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5, linetype = 4) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "RR"), aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "EN"), aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5, linetype = 4) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Method", size = 20) +
  ggtitle(sprintf("2 Causal Regions, 1 Shared Regions, Size 10")) +
  coord_cartesian(xlim=c(0, 1.0)) + 
#  geom_abline(intercept = 0, slope = 1,alpha = 0.5) + 
#  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 14), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)
#ggsave('~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_limit_shapes2_truncated.pdf')
ggsave('~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_limit_shapes2.pdf')
write.csv(roc_curve_frame,file = '~/Documents/SINATRA/Scripts/Data/sphere_roc_2causal_1shared2.csv',row.names = FALSE)
#### Fig 2.b. ####
num_causal_region = 6
num_shared_region = 3
causal_points = 10
shared_points = 10
cusps = 50
#Generate Data
#g2 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
#                                                 causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
#                                                 truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
#                                                 min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
#                                                 mode = 'sphere_baseline',num_cusps = cusps,
#                                                 subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 0.5,reduce = max, write = TRUE)
#

#write.csv(g2,'~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_elastic_net_max.csv')
path = '~/Documents/spheres6_3_2020-05-15'
dirs = list.dirs(path, full.names = TRUE, FALSE)
roc_curves2.1 = list()
roc_curves2.2 = list()
num_vertices = 642
functions_to_summarize = 10
for (i in 1:length(dirs)){
  dir = dirs[i]
  print(dir)
  total_rate_roc1 = matrix(0, nrow = num_vertices,ncol = 2)
  total_rate_roc2 = matrix(0, nrow = num_vertices,ncol = 2)
  dir1 = paste(dir,'/v1/', sep = '')
  dir2 = paste(dir,'/v2/', sep = '')
  causal_points1_vec = read.csv(paste(dir,'/gp1.csv',sep = ''))
  causal_points2_vec = read.csv(paste(dir,'/gp2.csv',sep = ''))
  
  class_1_true_vertices = which(causal_points1_vec == 1)
  class_2_true_vertices = which(causal_points2_vec == 1)
  
  class_1_false_vertices = setdiff(1:num_vertices,class_1_true_vertices)
  class_2_false_vertices = setdiff(1:num_vertices,class_2_true_vertices)
  
  limit_shapes_output_files1 = list.files(path = dir1, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs1 = limit_shapes_output_files1[str_detect(limit_shapes_output_files1,'csv')]
  limit_shapes_csvs1 = limit_shapes_csvs1[str_detect(limit_shapes_csvs1,'scrambled') == FALSE]
  
  limit_shapes_output_files2 = list.files(path = dir2, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs2 = limit_shapes_output_files2[str_detect(limit_shapes_output_files2,'csv')]
  limit_shapes_csvs2 = limit_shapes_csvs2[str_detect(limit_shapes_csvs2,'scrambled_sphere') == FALSE]
#  print(limit_shapes_csvs2)
  
  
  #ROC curve for class 1
  
  counter = 0 
  for (file in limit_shapes_csvs1){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){
      
      #sink("/dev/null")
      rate_positive_vertices = which(rate_values >= threshold)
      
      rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)
      rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                    class_1_true_vertices,class_1_false_vertices))
      
      true_vertices = class_1_true_vertices
      false_vertices = class_1_false_vertices
    }
    total_rate_roc1 = total_rate_roc1 + rate_ROC
    
  }
  total_rate_roc1 = (total_rate_roc1 / counter)
  #ROC Curve for Class 2
  counter = 0 
  for (file in limit_shapes_csvs2){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){
      
      #sink("/dev/null")
      rate_positive_vertices = which(rate_values >= threshold)
      
      rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)
      rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                    class_2_true_vertices,class_2_false_vertices))
      
      true_vertices = class_2_true_vertices
      false_vertices = class_2_false_vertices
    }
    total_rate_roc2 = total_rate_roc2 + rate_ROC
    
  }
  total_rate_roc2 = (total_rate_roc2 / counter)
  
  roc_curves2.1[[i]] = total_rate_roc1
  roc_curves2.2[[i]] = total_rate_roc2
}

total_roc2.1 = matrix(0, nrow = dim(roc_curves2.1[[1]]), ncol = dim(roc_curves2.2[[1]])[2])
for (j in 1:length(roc_curves2.1)){
  total_roc2.1 = total_roc2.1 + roc_curves2.1[[j]]
}
total_roc2.1 = total_roc2.1/length(roc_curves2.1)
total_roc2.2 = matrix(0, nrow = dim(roc_curves2.1[[1]]), ncol = dim(roc_curves2.2[[1]])[2])
for (j in 1:length(roc_curves2.2)){
  total_roc2.2 = total_roc2.2 + roc_curves2.2[[j]]
}
total_roc2.2 = total_roc2.2/length(roc_curves2.2)
total_roc2.1 =  cbind(total_roc2.1,rep('Limit Shapes Class 1', dim(total_roc2.1)[1]))
total_roc2.2 =  cbind(total_roc2.2,rep('Limit Shapes Class 2', dim(total_roc2.2)[1]))
roc_curve_frame_limit2 = as.data.frame(rbind(total_roc2.1, total_roc2.2))
roc_curve_frame_limit2$V1 = as.numeric(as.character(roc_curve_frame_limit2$V1))
roc_curve_frame_limit2$V2 = as.numeric(as.character(roc_curve_frame_limit2$V2))

ggplot(roc_curve_frame_limit2, aes(x = V1,y = V2,group = V3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(V3) )) +
  geom_line(stat = "identity") +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Cones") +
  ggtitle(sprintf("6 Causal Regions, 3 Shared Regions, Size 10")) +
  geom_abline(intercept = 0, slope = 1) + 
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

save.image('~/Documents/spheres6_3_2020-05-15/roc_curves_10_function.Rdata')
load('~/Dropbox (Princeton)/SINATRA_Data/spheres6_3_2020-05-15/roc_curves_10_function.Rdata')

## Scrambled
path = '~/Documents/spheres6_3_2020-05-15'
dirs = list.dirs(path, full.names = TRUE, FALSE)
roc_curves2.1 = list()
roc_curves2.2 = list()
num_vertices = 642
functions_to_summarize = 10
for (i in 1:length(dirs)){
  dir = dirs[i]
  print(dir)
  total_rate_roc1 = matrix(0, nrow = num_vertices,ncol = 2)
  total_rate_roc2 = matrix(0, nrow = num_vertices,ncol = 2)
  dir1 = paste(dir,'/v1/', sep = '')
  dir2 = paste(dir,'/v2/', sep = '')
  causal_points1_vec = read.csv(paste(dir,'/gp1.csv',sep = ''))
  causal_points2_vec = read.csv(paste(dir,'/gp2.csv',sep = ''))
  
  class_1_true_vertices = which(causal_points1_vec == 1)
  class_2_true_vertices = which(causal_points2_vec == 1)
  
  class_1_false_vertices = setdiff(1:num_vertices,class_1_true_vertices)
  class_2_false_vertices = setdiff(1:num_vertices,class_2_true_vertices)
  
  limit_shapes_output_files1 = list.files(path = dir1, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs1 = limit_shapes_output_files1[str_detect(limit_shapes_output_files1,'csv')]
  limit_shapes_csvs1 = limit_shapes_csvs1[str_detect(limit_shapes_csvs1,'scrambled')]
  
  limit_shapes_output_files2 = list.files(path = dir2, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs2 = limit_shapes_output_files2[str_detect(limit_shapes_output_files2,'csv')]
  limit_shapes_csvs2 = limit_shapes_csvs2[str_detect(limit_shapes_csvs2,'scrambled')]
  
  #ROC curve for class 1
  
  counter = 0 
  for (file in limit_shapes_csvs1){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){
      
      #sink("/dev/null")
      rate_positive_vertices = which(rate_values >= threshold)
      
      rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)
      rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                    class_1_true_vertices,class_1_false_vertices))
      
      true_vertices = class_1_true_vertices
      false_vertices = class_1_false_vertices
    }
    total_rate_roc1 = total_rate_roc1 + rate_ROC
    
  }
  total_rate_roc1 = (total_rate_roc1 / counter)
  #ROC Curve for Class 2
  counter = 0 
  for (file in limit_shapes_csvs2){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){
      
      #sink("/dev/null")
      rate_positive_vertices = which(rate_values >= threshold)
      
      rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)
      rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                    class_2_true_vertices,class_2_false_vertices))
      
      true_vertices = class_2_true_vertices
      false_vertices = class_2_false_vertices
    }
    total_rate_roc2 = total_rate_roc2 + rate_ROC
    
  }
  total_rate_roc2 = (total_rate_roc2 / counter)
  
  roc_curves2.1[[i]] = total_rate_roc1
  roc_curves2.2[[i]] = total_rate_roc2
}
#roc_curves1.1 = roc_curves1
#roc_curves1.2 = roc_curves2

total_roc1 = matrix(0, nrow = dim(roc_curves2.1[[1]]), ncol = dim(roc_curves2.2[[1]])[2])
for (j in 1:length(roc_curves2.1)){
  total_roc1 = total_roc1 + roc_curves2.1[[j]]
}
total_roc1 = total_roc1/length(roc_curves2.1)
total_roc2 = matrix(0, nrow = dim(roc_curves2.1[[1]]), ncol = dim(roc_curves2.2[[1]])[2])
for (j in 1:length(roc_curves2.2)){
  total_roc2 = total_roc2 + roc_curves2.2[[j]]
}
total_roc2 = total_roc2/length(roc_curves2.2)
total_roc1 =  cbind(total_roc1,rep('Limit Shapes Class 1', dim(total_roc1)[1]))
total_roc2 =  cbind(total_roc2,rep('Limit Shapes Class 2', dim(total_roc2)[1]))
roc_curve_frame_limit_scrambled2 = as.data.frame(rbind(total_roc1, total_roc2))
roc_curve_frame_limit_scrambled2$V1 = as.numeric(as.character(roc_curve_frame_limit_scrambled2$V1))
roc_curve_frame_limit_scrambled2$V2 = as.numeric(as.character(roc_curve_frame_limit_scrambled2$V2))

save.image('~/Documents/spheres6_3_2020-05-15/scrambled_roc_curves_10_function.Rdata')
load('~/Dropbox (Princeton)/SINATRA_Data/spheres6_3_2020-05-15/scrambled_roc_curves_10_function.Rdata')
g2 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_elastic_net_max2.csv')[,-1][,-4]
roc_curve_frame2.1 = data.frame(g2)
library(ggplot2)
roc_curve_frame2.1 = roc_curve_frame2.1[as.numeric(as.character(roc_curve_frame2.1[,3])) == 1,]
roc_curve_frame2.1[,3] = rep('Baseline (EN Max)',dim(roc_curve_frame2.1)[1])
#g2.2 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
#                                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
#                                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
#                                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
#                                                   mode = 'sphere_baseline',num_cusps = cusps,
#                                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 0.5,reduce = mean)
#
#
#write.csv(g2.2,'~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_elastic_net_mean.csv')
g2.2 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_elastic_net_mean2.csv')[,-1][,-4]
roc_curve_frame2.2 = data.frame(g2.2)
library(ggplot2)
roc_curve_frame2.2 = roc_curve_frame2.2[roc_curve_frame2.2[,3] == 1,]
roc_curve_frame2.2[,3] = rep('Baseline (EN Mean)',dim(roc_curve_frame2.2)[1])
#g2.3 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
#                                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
#                                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
#                                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
#                                                   mode = 'sphere_baseline',num_cusps = cusps,
#                                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 1,reduce = max)
#
#
#write.csv(g2.3,'~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_ridge_max.csv')
g2.3 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_ridge_max2.csv')[,-1][,-4]
roc_curve_frame2.3 = data.frame(g2.3)
library(ggplot2)
roc_curve_frame2.3 = roc_curve_frame2.3[roc_curve_frame2.3[,3] == 1,]
roc_curve_frame2.3[,3] = rep('Baseline (RR Max)',dim(roc_curve_frame2.3)[1])
#g2.4 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
#                                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
#                                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
#                                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
#                                                   mode = 'sphere_baseline',num_cusps = cusps,
#                                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 1,reduce = mean)
#
#
#write.csv(g2.4,'~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_ridge_mean.csv')
g2.4 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_ridge_mean2.csv')[,-1][,-4]
roc_curve_frame2.4 = data.frame(g2.4)
library(ggplot2)
roc_curve_frame2.4 = roc_curve_frame2.4[roc_curve_frame2.4[,3] == 1,]
roc_curve_frame2.4[,3] = rep('Baseline (RR Mean)',dim(roc_curve_frame2.4)[1])
roc_curve_frame = rbind(roc_curve_frame2.1,roc_curve_frame2.2,roc_curve_frame2.3,roc_curve_frame2.4)
library(ggplot2)
load('~/Documents/SINATRA/Simulations/Results/new/df_ROC_causal3_shared6_10.RData')
rdfmeans = rdfmeans[rdfmeans$Class == 1,]
rdfmeans = rdfmeans[,-2]
rdfmeans = rdfmeans[,-2]
rdfmeans = rdfmeans[rdfmeans$Num_Cones==60,]
rdfmeans = rdfmeans[c('FPR','TPR')]
rdfmeans[,3] = rep('SINATRA',dim(rdfmeans)[1])
roc_curve_frame = rbind(roc_curve_frame2.1,roc_curve_frame2.2,roc_curve_frame2.3,roc_curve_frame2.4)
names(rdfmeans) = names(roc_curve_frame)
roc_curve_frame = rbind(roc_curve_frame,rdfmeans)
roc_curve_frame_limit2 = roc_curve_frame_limit2[roc_curve_frame_limit2[,3] == 'Limit Shapes Class 1',]
roc_curve_frame_limit2[,3] = 'Limit Shapes'
#roc_curve_frame_limit_scrambled2 = roc_curve_frame_limit_scrambled2[roc_curve_frame_limit_scrambled2[,3] == 'Limit Shapes Class 1',]
roc_curve_frame_limit_scrambled2[,3] = 'Limit Shapes (Misspecified)'
roc_curve_frame = rbind(roc_curve_frame,roc_curve_frame_limit2)
roc_curve_frame = rbind(roc_curve_frame,roc_curve_frame_limit_scrambled2)
ggplot(roc_curve_frame, aes(x = V1,y = V2,group = V3))  +
  geom_line(data = subset(roc_curve_frame, V3 == "SINATRA"), aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 2) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "Limit"), aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5, linetype = 4) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "RR"), aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "EN"), aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5, linetype = 4) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Method") +
  coord_cartesian(xlim=c(0, 1.0)) + 
  ggtitle(sprintf("6 Causal Regions, 3 Shared Regions, Size 10")) +
  #geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
#  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 14), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)
#ggsave('~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared2_truncated.pdf')
ggsave('~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared2.pdf')
write.csv(roc_curve_frame,file = '~/Documents/SINATRA/Scripts/Data/sphere_roc_6causal_3shared2.csv',row.names = FALSE)
#### Fig 2.c. ####
num_causal_region = 10
num_shared_region = 5
causal_points = 10
shared_points = 10
#g3 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
#                                                 causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
#                                                 truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
#                                                 min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
#                                                 mode = 'sphere_baseline',num_cusps = cusps,
#                                                 subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 0.5,reduce = max, write = TRUE)

path = '~/Documents/spheres10_5_2020-05-15'
dirs = list.dirs(path, full.names = TRUE, FALSE)
roc_curves3.1 = list()
roc_curves3.2 = list()
num_vertices = 642
functions_to_summarize = 10
for (i in 1:length(dirs)){
  dir = dirs[i]
  print(dir)
  total_rate_roc1 = matrix(0, nrow = num_vertices,ncol = 2)
  total_rate_roc2 = matrix(0, nrow = num_vertices,ncol = 2)
  dir1 = paste(dir,'/v1/', sep = '')
  dir2 = paste(dir,'/v2/', sep = '')
  causal_points1_vec = read.csv(paste(dir,'/gp1.csv',sep = ''))
  causal_points2_vec = read.csv(paste(dir,'/gp2.csv',sep = ''))
  
  class_1_true_vertices = which(causal_points1_vec == 1)
  class_2_true_vertices = which(causal_points2_vec == 1)
  
  class_1_false_vertices = setdiff(1:num_vertices,class_1_true_vertices)
  class_2_false_vertices = setdiff(1:num_vertices,class_2_true_vertices)
  
  limit_shapes_output_files1 = list.files(path = dir1, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs1 = limit_shapes_output_files1[str_detect(limit_shapes_output_files1,'csv')]
  limit_shapes_csvs1 = limit_shapes_csvs1[str_detect(limit_shapes_csvs1,'scrambled')== FALSE]
  
  limit_shapes_output_files2 = list.files(path = dir2, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs2 = limit_shapes_output_files2[str_detect(limit_shapes_output_files2,'csv')]
  limit_shapes_csvs2 = limit_shapes_csvs2[str_detect(limit_shapes_csvs2,'scrambled') == FALSE]
  
  #ROC curve for class 1
  
  counter = 0 
  for (file in limit_shapes_csvs1){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){
      
      #sink("/dev/null")
      rate_positive_vertices = which(rate_values >= threshold)
      
      rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)
      rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                    class_1_true_vertices,class_1_false_vertices))
      
      true_vertices = class_1_true_vertices
      false_vertices = class_1_false_vertices
    }
    total_rate_roc1 = total_rate_roc1 + rate_ROC
    
  }
  total_rate_roc1 = (total_rate_roc1 / counter)
  #ROC Curve for Class 2
  counter = 0 
  for (file in limit_shapes_csvs2){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){
      
      #sink("/dev/null")
      rate_positive_vertices = which(rate_values >= threshold)
      
      rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)
      rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                    class_2_true_vertices,class_2_false_vertices))
      
      true_vertices = class_2_true_vertices
      false_vertices = class_2_false_vertices
    }
    total_rate_roc2 = total_rate_roc2 + rate_ROC
    
  }
  total_rate_roc2 = (total_rate_roc2 / counter)
  
  roc_curves3.1[[i]] = total_rate_roc1
  roc_curves3.2[[i]] = total_rate_roc2
}

total_roc3.1 = matrix(0, nrow = dim(roc_curves3.1[[1]]), ncol = dim(roc_curves3.2[[1]])[2])
for (j in 1:length(roc_curves3.1)){
  total_roc3.1 = total_roc3.1 + roc_curves3.1[[j]]
}
total_roc3.1 = total_roc3.1/length(roc_curves3.1)
total_roc3.2 = matrix(0, nrow = dim(roc_curves3.1[[1]]), ncol = dim(roc_curves3.2[[1]])[2])
for (j in 1:length(roc_curves3.2)){
  total_roc3.2 = total_roc3.2 + roc_curves3.2[[j]]
}
total_roc3.2 = total_roc3.2/length(roc_curves3.2)
total_roc3.1 =  cbind(total_roc3.1,rep('Limit Shapes Class 1', dim(total_roc3.1)[1]))
total_roc3.2 =  cbind(total_roc3.2,rep('Limit Shapes Class 2', dim(total_roc3.2)[1]))
roc_curve_frame_limit3 = as.data.frame(rbind(total_roc3.1, total_roc3.2))
roc_curve_frame_limit3$V1 = as.numeric(as.character(roc_curve_frame_limit3$V1))
roc_curve_frame_limit3$V2 = as.numeric(as.character(roc_curve_frame_limit3$V2))

#ggplot(roc_curve_frame_limit3, aes(x = V1,y = V2,group = V3)) +
#  geom_line(data = subset(roc_curve_frame, V3 == "SINATRA"), aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 2) +
#  geom_line(data = subset(roc_curve_frame, V3 != "SINATRA"), aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 2, linetype = 2) +
#  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Cones") +
#  ggtitle(sprintf("10 Causal Regions, 5 Shared Regions, Size 10")) +
#  geom_abline(intercept = 0, slope = 1) + 
#  coord_equal(ratio=1) +
#  theme_bw() +
#  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

save.image('~/Documents/spheres10_5_2020-05-15/roc_curves_10_function.Rdata')
load('~/Dropbox (Princeton)/SINATRA_Data/spheres10_5_2020-05-15/roc_curves_10_function.Rdata')

path = '~/Documents/spheres10_5_2020-05-15'
dirs = list.dirs(path, full.names = TRUE, FALSE)
roc_curves3.1 = list()
roc_curves3.2 = list()
num_vertices = 642
functions_to_summarize = 10
for (i in 1:length(dirs)){
  dir = dirs[i]
  print(dir)
  total_rate_roc1 = matrix(0, nrow = num_vertices,ncol = 2)
  total_rate_roc2 = matrix(0, nrow = num_vertices,ncol = 2)
  dir1 = paste(dir,'/v1/', sep = '')
  dir2 = paste(dir,'/v2/', sep = '')
  causal_points1_vec = read.csv(paste(dir,'/gp1.csv',sep = ''))
  causal_points2_vec = read.csv(paste(dir,'/gp2.csv',sep = ''))
  
  class_1_true_vertices = which(causal_points1_vec == 1)
  class_2_true_vertices = which(causal_points2_vec == 1)
  
  class_1_false_vertices = setdiff(1:num_vertices,class_1_true_vertices)
  class_2_false_vertices = setdiff(1:num_vertices,class_2_true_vertices)
  
  limit_shapes_output_files1 = list.files(path = dir1, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs1 = limit_shapes_output_files1[str_detect(limit_shapes_output_files1,'csv')]
  limit_shapes_csvs1 = limit_shapes_csvs1[str_detect(limit_shapes_csvs1,'scrambled')]
  
  limit_shapes_output_files2 = list.files(path = dir2, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs2 = limit_shapes_output_files2[str_detect(limit_shapes_output_files2,'csv')]
  limit_shapes_csvs2 = limit_shapes_csvs2[str_detect(limit_shapes_csvs2,'scrambled')]
  
  #ROC curve for class 1
  
  counter = 0 
  for (file in limit_shapes_csvs1){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){
      
      #sink("/dev/null")
      rate_positive_vertices = which(rate_values >= threshold)
      
      rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)
      rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                    class_1_true_vertices,class_1_false_vertices))
      
      true_vertices = class_1_true_vertices
      false_vertices = class_1_false_vertices
    }
    total_rate_roc1 = total_rate_roc1 + rate_ROC
    
  }
  total_rate_roc1 = (total_rate_roc1 / counter)
  #ROC Curve for Class 2
  counter = 0 
  for (file in limit_shapes_csvs2){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = length(rate_values))) ){
      
      #sink("/dev/null")
      rate_positive_vertices = which(rate_values >= threshold)
      
      rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)
      rate_ROC <- rbind(rate_ROC, calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                                    class_2_true_vertices,class_2_false_vertices))
      
      true_vertices = class_2_true_vertices
      false_vertices = class_2_false_vertices
    }
    total_rate_roc2 = total_rate_roc2 + rate_ROC
    
  }
  total_rate_roc2 = (total_rate_roc2 / counter)
  
  roc_curves3.1[[i]] = total_rate_roc1
  roc_curves3.2[[i]] = total_rate_roc2
}

total_roc3.1 = matrix(0, nrow = dim(roc_curves3.1[[1]]), ncol = dim(roc_curves3.2[[1]])[2])
for (j in 1:length(roc_curves3.1)){
  total_roc3.1 = total_roc3.1 + roc_curves3.1[[j]]
}
total_roc3.1 = total_roc3.1/length(roc_curves3.1)
total_roc3.2 = matrix(0, nrow = dim(roc_curves3.1[[1]]), ncol = dim(roc_curves3.2[[1]])[2])
for (j in 1:length(roc_curves3.2)){
  total_roc3.2 = total_roc3.2 + roc_curves3.2[[j]]
}
total_roc3.2 = total_roc3.2/length(roc_curves3.2)
total_roc3.1 =  cbind(total_roc3.1,rep('Limit Shapes Class 1', dim(total_roc3.1)[1]))
total_roc3.2 =  cbind(total_roc3.2,rep('Limit Shapes Class 2', dim(total_roc3.2)[1]))
roc_curve_frame_limit_scrambled3 = as.data.frame(rbind(total_roc3.1, total_roc3.2))
roc_curve_frame_limit_scrambled3$V1 = as.numeric(as.character(roc_curve_frame_limit_scrambled3$V1))
roc_curve_frame_limit_scrambled3$V2 = as.numeric(as.character(roc_curve_frame_limit_scrambled3$V2))

ggplot(roc_curve_frame_limit3, aes(x = V1,y = V2,group = V3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(V3) )) +
  geom_line(stat = "identity") +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Cones") +
  ggtitle(sprintf("10 Causal Regions, 5 Shared Regions, Size 10")) +
  geom_abline(intercept = 0, slope = 1) + 
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 

save.image('~/Documents/spheres10_5_2020-05-15/scrambled_roc_curves_10_function.Rdata')
load('~/Dropbox (Princeton)/SINATRA_Data/spheres10_5_2020-05-15/scrambled_roc_curves_10_function.Rdata')



#write.csv(g3,'~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_elastic_net_max.csv')
g3 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_elastic_net_max2.csv')[,-1][,-4]
roc_curve_frame3.1 = as.data.frame(g3)
library(ggplot2)
roc_curve_frame3.1 = roc_curve_frame3.1[as.numeric(as.character(roc_curve_frame3.1[,3])) == 1,]

roc_curve_frame3.1[,3] = rep('Baseline (EN Max)',dim(roc_curve_frame3.1)[1])
#g3.2 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
#                                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
#                                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
#                                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
#                                                   mode = 'sphere_baseline',num_cusps = cusps,
#                                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 0.5,reduce = mean)
#

#write.csv(g3.2,'~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_elastic_net_mean.csv')
g3.2 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_elastic_net_mean2.csv')[,-1][,-4]
roc_curve_frame3.2 = data.frame(g3.2)
library(ggplot2)
roc_curve_frame3.2 = roc_curve_frame3.2[roc_curve_frame3.2[,3] == 1,]
roc_curve_frame3.2[,3] = rep('Baseline (EN Mean)',dim(roc_curve_frame3.2)[1])
#g3.3 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
#                                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
#                                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
#                                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
#                                                   mode = 'sphere_baseline',num_cusps = cusps,
#                                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 1,reduce = max)
#

#write.csv(g3.3,'~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_ridge_max.csv')
g3.3 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_ridge_max2.csv')[,-1][,-4]
roc_curve_frame3.3 = data.frame(g3.3)
library(ggplot2)
roc_curve_frame3.3 = roc_curve_frame3.3[roc_curve_frame3.3[,3] == 1,]
roc_curve_frame3.3[,3] = rep('Baseline (RR Max)',dim(roc_curve_frame3.3)[1])
#g3.4 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
#                                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
#                                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
#                                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
#                                                   mode = 'sphere_baseline',num_cusps = cusps,
#                                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 1,reduce = mean)
#
#
#write.csv(g3.4,'~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_ridge_mean.csv')
g3.4 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_ridge_mean2.csv')[,-1][,-4]
roc_curve_frame3.4 = data.frame(g3.4)
library(ggplot2)
roc_curve_frame3.4 = roc_curve_frame3.4[roc_curve_frame3.4[,3] == 1,]
roc_curve_frame3.4[,3] = rep('Baseline (RR Mean)',dim(roc_curve_frame3.4)[1])
roc_curve_frame = rbind(roc_curve_frame3.1,roc_curve_frame3.2,roc_curve_frame3.3,roc_curve_frame3.4)
library(ggplot2)
load('~/Documents/SINATRA/Simulations/Results/new/df_ROC_causal5_shared10_10.RData')
rdfmeans = rdfmeans[rdfmeans$Class == 1,]
rdfmeans = rdfmeans[,-2]
rdfmeans = rdfmeans[,-2]
rdfmeans = rdfmeans[rdfmeans$Num_Cones==60,]
rdfmeans = rdfmeans[c('FPR','TPR')]
rdfmeans[,3] = rep('SINATRA',dim(rdfmeans)[1])
names(rdfmeans) = names(roc_curve_frame)
roc_curve_frame = rbind(roc_curve_frame,rdfmeans)
roc_curve_frame_limit3 = roc_curve_frame_limit3[roc_curve_frame_limit3[,3] == 'Limit Shapes Class 1',]
roc_curve_frame_limit3[,3] = 'Limit Shapes'
roc_curve_frame_limit_scrambled3 = roc_curve_frame_limit_scrambled3[roc_curve_frame_limit_scrambled3[,3] == 'Limit Shapes Class 1',]
roc_curve_frame_limit_scrambled3[,3] = 'Limit Shapes (Misspecified)'
roc_curve_frame = rbind(roc_curve_frame,roc_curve_frame_limit3)
roc_curve_frame = rbind(roc_curve_frame,roc_curve_frame_limit_scrambled3)
ggplot(roc_curve_frame, aes(x = V1,y = V2,group = V3)) + 
  geom_line(data = subset(roc_curve_frame, V3 == "SINATRA"), aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "Limit"), aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5, linetype = 4) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "RR"), aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "EN"), aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5, linetype = 4) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Method") +
  ggtitle(sprintf("10 Causal Regions, 5 Shared Regions, Size 10")) +
  coord_cartesian(xlim=c(0, 1.0)) + 
  #geom_abline(intercept = 0, slope = 1, alpha =0.5) + 
  #coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 14), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)
ggsave('~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared2.pdf')
write.csv(roc_curve_frame,file = '~/Documents/SINATRA/Scripts/Data/sphere_roc_10causal_5shared2.csv',row.names = FALSE)
