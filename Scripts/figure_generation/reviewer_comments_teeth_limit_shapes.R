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
roc_curve_teeth_baseline = function(data_dir,var_selection='elastic_net',ec_type = 'baseline',reduce = max,alpha = 0.5, truncated = 500){
  roc_curves1 = list()
  roc_curves2 = list()
  data_dirs = list.dirs(data_dir,recursive = FALSE)
  for (k in 1:length(data_dirs)){
      dir = data_dirs[k]
      dir1 = paste(dir,'/mesh/gp1',sep='')
      dir2 = paste(dir,'/mesh/gp2',sep='')
      data_files1 = list.files(dir1, full.names = TRUE)
      data_files1 = data_files1[str_detect(data_files1,'off')]
      data_files2 = list.files(dir2, full.names = TRUE)
      data_files2 = data_files1[str_detect(data_files2,'off')]
      class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE)
      class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)
      roc_curves = list()
      dirs = 0
      pset = list(num_cones = 0, cap_radius = 0.15, len = 100, directions_per_cone = 5,
                    dirs = 0)
      data_summary=real_data_summary(dir1=dir1,dir2 = dir2,direction=pset$dirs,
                                      len = pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type, reduce = reduce, alpha = alpha, mode = var_selection)
      roc_curve = compute_roc_curve_teeth(data_dir1 = dir1, data_dir2 = dir2, gamma = 0.25,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                            rate_values = data_summary$Rate2,directions_per_cone = pset$directions_per_cone,curve_length = pset$len,
                                            directions = pset$dirs,truncated = truncated,ball_radius = ball_radius, ball = ball, radius = 1, two_curves = TRUE,mode = 'baseline')
        #print(roc_curve)
#        roc_curve[,3] = dirs
        roc_frame = data.frame(roc_curve)
        roc_curves1[[k]] = roc_frame[1:truncated,]
        roc_curves2[[k]] = roc_frame[truncated:(2*truncated),]
  }
  roc_curve1 = summarize_list(roc_curves1)
  roc_curve2 = summarize_list(roc_curves2)
  return(list(roc_curves1,roc_curves2, roc_curve1, roc_curve2))
}

#### Fig 2.a. ####

base_dir = '~/Documents/new_aligned_shapesv3/'

#g1 = roc_curve_teeth_baseline(data_dir = base_dir,var_selection = 'elastic_net',ec_type = 'baseline',reduce = max, alpha = 0.5)
#write.csv(g1[[3]],'~/Documents/SINATRA/Scripts/Data/cariature_v3_elastic_net_max.csv')
g1 = read.csv('~/Documents/SINATRA/Scripts/Data/cariature_v3_elastic_net_max.csv')[,-1][,-4]
roc_curve_frame1.1 = data.frame(g1)
library(ggplot2)
roc_curve_frame1.1 = roc_curve_frame1.1[roc_curve_frame1.1[,3] == 1,]
roc_curve_frame1.1[,3] = rep('Baseline (EN Max)',dim(roc_curve_frame1.1)[1])
#g1.2 = roc_curve_teeth_baseline(data_dir = base_dir,var_selection = 'elastic_net',ec_type = 'baseline',reduce = mean, alpha = 0.5)
#write.csv(g1.2[[3]],'~/Documents/SINATRA/Scripts/Data/cariature_v3_elastic_net_mean.csv')
g1.2 = read.csv('~/Documents/SINATRA/Scripts/Data/cariature_v3_elastic_net_mean.csv')[,-1][,-4]
roc_curve_frame1.2 = data.frame(g1.2)
library(ggplot2)
roc_curve_frame1.2 = roc_curve_frame1.2[roc_curve_frame1.2[,3] == 1,]
roc_curve_frame1.2[,3] = rep('Baseline (EN Mean)',dim(roc_curve_frame1.2)[1])
#g1.3 = roc_curve_teeth_baseline(data_dir = base_dir,var_selection = 'elastic_net',ec_type = 'baseline',reduce = max, alpha = 1)
#write.csv(g1.3[[3]],'~/Documents/SINATRA/Scripts/Data/cariature_v3_ridge_net_max.csv')
g1.3 = read.csv('~/Documents/SINATRA/Scripts/Data/cariature_v3_ridge_net_max.csv')[,-1][,-4]
roc_curve_frame1.3 = data.frame(g1.3)
library(ggplot2)
roc_curve_frame1.3 = roc_curve_frame1.3[roc_curve_frame1.3[,3] == 1,]
roc_curve_frame1.3[,3] = rep('Baseline (RR Max)',dim(roc_curve_frame1.3)[1])
#g1.4 = roc_curve_teeth_baseline(data_dir = base_dir,var_selection = 'elastic_net',ec_type = 'baseline',reduce = mean, alpha = 1)
#write.csv(g1.4[[3]],'~/Documents/SINATRA/Scripts/Data/cariature_v3_ridge_net_mean.csv')
g1.4 = read.csv('~/Documents/SINATRA/Scripts/Data/cariature_v3_ridge_net_mean.csv')[,-1][,-4]
roc_curve_frame1.4 = data.frame(g1.4)
library(ggplot2)
roc_curve_frame1.4 = roc_curve_frame1.4[roc_curve_frame1.4[,3] == 1,]
roc_curve_frame1.4[,3] = rep('Baseline (RR Mean)',dim(roc_curve_frame1.4)[1])
library(ggplot2)
#### Load in the SINATRA Curves ####
#data_dirs = list.dirs(base_dir,recursive = FALSE)
#rdfmeans = read.csv(paste(data_dirs[1],'/roc_dirs1.csv', sep = ''))[,-4]
#for (k in 2:length(data_dirs)){
#  dir = data_dirs[k]
#  temp = read.csv(paste(dir,'/roc_dirs1.csv', sep = ''))[,-4]
#  rdfmeans = rdfmeans + temp
#}
#rdfmeans = rdfmeans/length(data_dirs)
#rdfmeans = data.frame(rdfmeans)
#rdfmeans = rdfmeans[rdfmeans$X3==35,]
#rdfmeans[,3] = rep('SINATRA',dim(rdfmeans)[1])
#names(rdfmeans) = names(roc_curve_frame)
#roc_curve_frame = rbind(roc_curve_frame,rdfmeans)


### Limit Shapes ####
path = '~/Documents/new_aligned_shapesv3'
dirs = list.dirs(path, full.names = TRUE, FALSE)
roc_curves1.1 = list()
roc_curves1.2 = list()
truncated = 500
functions_to_summarize = 10
for (i in 1:length(dirs)){
  dir = dirs[i]
  print(dir)
  total_rate_roc1 = matrix(0, nrow = truncated,ncol = 2)
  total_rate_roc2 = matrix(0, nrow = truncated,ncol = 2)
  dir1 = paste(dir,'/mesh/gp1/', sep = '')
  dir2 = paste(dir,'/mesh/gp2/', sep = '')
  causal_points1_vec = read.csv(paste(dir,'/gp1_spt.csv',sep = ''))
  causal_points2_vec = read.csv(paste(dir,'/gp2_spt.csv',sep = ''))
  
  class_1_true_vertices = which(causal_points1_vec > 0.25)
  class_2_true_vertices = which(causal_points2_vec > 0.25)
  
  class_1_false_vertices = setdiff(1:num_vertices,class_1_true_vertices)
  class_2_false_vertices = setdiff(1:num_vertices,class_2_true_vertices)
  
  limit_shapes_output_files1 = list.files(path = dir1, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs1 = limit_shapes_output_files1[str_detect(limit_shapes_output_files1,'csv')]
  limit_shapes_csvs1 = limit_shapes_csvs1[str_detect(limit_shapes_csvs1,'scrambled')== FALSE]
  
  limit_shapes_output_files2 = list.files(path = dir2, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs2 = limit_shapes_output_files2[str_detect(limit_shapes_output_files2,'csv')]
  limit_shapes_csvs2 = limit_shapes_csvs2[str_detect(limit_shapes_csvs2,'scrambled')== FALSE]
  
  #ROC curve for class 1
  
  counter = 0 
  for (file in limit_shapes_csvs1){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = truncated))){
      
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
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = truncated)) ){
      
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
total_roc1 = rbind(c(0,0), total_roc1)
total_roc2 = rbind(c(0,0), total_roc2)
total_roc1 =  cbind(total_roc1,rep('Limit Shapes Class 1', dim(total_roc1)[1]))
total_roc2 =  cbind(total_roc2,rep('Limit Shapes Class 2', dim(total_roc2)[1]))
roc_curve_frame_limit = as.data.frame(rbind(total_roc1, total_roc2))
roc_curve_frame_limit$V1 = as.numeric(as.character(roc_curve_frame_limit$V1))
roc_curve_frame_limit$V2 = as.numeric(as.character(roc_curve_frame_limit$V2))
save.image('~/Dropbox (Princeton)//new_aligned_shapesv3/limitshapes.Rdata')
load('~/Dropbox (Princeton)/SINATRA_Data/new_aligned_shapesv3/limitshapes.Rdata')

### Limit Shapes ####
path = '~/Documents/new_aligned_shapesv3'
dirs = list.dirs(path, full.names = TRUE, FALSE)
roc_curves1.1 = list()
roc_curves1.2 = list()
truncated = 500
functions_to_summarize = 10
for (i in 1:length(dirs)){
  dir = dirs[i]
  print(dir)
  total_rate_roc1 = matrix(0, nrow = truncated,ncol = 2)
  total_rate_roc2 = matrix(0, nrow = truncated,ncol = 2)
  dir1 = paste(dir,'/mesh/gp1/', sep = '')
  dir2 = paste(dir,'/mesh/gp2/', sep = '')
  causal_points1_vec = read.csv(paste(dir,'/gp1_spt.csv',sep = ''))
  causal_points2_vec = read.csv(paste(dir,'/gp2_spt.csv',sep = ''))
  
  class_1_true_vertices = which(causal_points1_vec > 0.25)
  class_2_true_vertices = which(causal_points2_vec > 0.25)
  
  class_1_false_vertices = setdiff(1:num_vertices,class_1_true_vertices)
  class_2_false_vertices = setdiff(1:num_vertices,class_2_true_vertices)
  
  limit_shapes_output_files1 = list.files(path = dir1, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs1 = limit_shapes_output_files1[str_detect(limit_shapes_output_files1,'csv')]
  limit_shapes_csvs1 = limit_shapes_csvs1[str_detect(limit_shapes_csvs1,'scrambled')== TRUE]
  
  limit_shapes_output_files2 = list.files(path = dir2, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs2 = limit_shapes_output_files2[str_detect(limit_shapes_output_files2,'csv')]
  limit_shapes_csvs2 = limit_shapes_csvs2[str_detect(limit_shapes_csvs2,'scrambled')== TRUE]
  
  #ROC curve for class 1
  
  counter = 0 
  for (file in limit_shapes_csvs1){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = truncated))){
      
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
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = truncated)) ){
      
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
total_roc1 = rbind(c(0,0), total_roc1)
total_roc2 = rbind(c(0,0), total_roc2)
total_roc1 =  cbind(total_roc1,rep('Limit Shapes Class 1', dim(total_roc1)[1]))
total_roc2 =  cbind(total_roc2,rep('Limit Shapes Class 2', dim(total_roc2)[1]))
roc_curve_frame_limit_scrambled = as.data.frame(rbind(total_roc1, total_roc2))
roc_curve_frame_limit_scrambled$V1 = as.numeric(as.character(roc_curve_frame_limit_scrambled$V1))
roc_curve_frame_limit_scrambled$V2 = as.numeric(as.character(roc_curve_frame_limit_scrambled$V2))
save.image('~/Documents/new_aligned_shapesv3/limitshapes_scrambled.Rdata')
load('~/Dropbox (Princeton)/SINATRA_Data/new_aligned_shapesv3/limitshapes_scrambled.Rdata')
roc_curve_frame = rbind(roc_curve_frame1.1,roc_curve_frame1.2,roc_curve_frame1.3,roc_curve_frame1.4)
roc_curve_frame = rbind(roc_curve_frame,rdfmeans)
roc_curve_frame_limit2 = roc_curve_frame_limit[roc_curve_frame_limit[,3] == 'Limit Shapes Class 1',]
roc_curve_frame_limit2[,3] = 'Limit Shapes'
roc_curve_frame_limit_scrambled2 = roc_curve_frame_limit_scrambled[roc_curve_frame_limit_scrambled[,3] == 'Limit Shapes Class 1',]
roc_curve_frame_limit_scrambled2[,3] = 'Limit Shapes (Misspecified)'
colnames(roc_curve_frame_limit2) = colnames(roc_curve_frame)
colnames(roc_curve_frame_limit_scrambled2) = colnames(roc_curve_frame)
roc_curve_frame = rbind(roc_curve_frame,roc_curve_frame_limit2)
roc_curve_frame = rbind(roc_curve_frame,roc_curve_frame_limit_scrambled2)

ggplot() + 
  geom_line(data = subset(roc_curve_frame, X3 == "SINATRA"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5) +
  geom_line(data = subset(roc_curve_frame, X3 %like% "Limit"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5, linetype = 4) +
  geom_line(data = subset(roc_curve_frame, X3 %like% "EN"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
#  geom_line(data = subset(roc_curve_frame, X3 %like% "Max"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
  # geom_line(stat = "identity",aes(color = factor(X3)), linetype = 'dotted') +
  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", color = "Method") +
  ggtitle(sprintf("3 Caricatured Peaks")) +
  coord_cartesian(xlim= c(0,1.0))+
#  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
#  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)
ggsave('~/Documents/SINATRA/Scripts/Data/cariature_v3_limit.pdf')
write.csv(roc_curve_frame,file = '~/Documents/SINATRA/Scripts/Data/3peaks_caricature_roc.csv',row.names = FALSE)
#### 5 Peaks ####
base_dir = '~/Documents/new_aligned_shapesv4/'

#g2 = roc_curve_teeth_baseline(data_dir = base_dir,var_selection = 'elastic_net',ec_type = 'baseline',reduce = max, alpha = 0.5, truncated = 500)
#write.csv(g2[[3]],'~/Documents/SINATRA/Scripts/Data/cariature_v4_elastic_net_max.csv')
g2 = read.csv('~/Documents/SINATRA/Scripts/Data/cariature_v4_elastic_net_max.csv')[,-1][,-4]
roc_curve_frame2.1 = data.frame(g2)
library(ggplot2)
roc_curve_frame2.1 = roc_curve_frame2.1[roc_curve_frame2.1[,3] == 1,]
roc_curve_frame2.1[,3] = rep('Baseline (EN Max)',dim(roc_curve_frame2.1)[1])
#g2.2 = roc_curve_teeth_baseline(data_dir = base_dir,var_selection = 'elastic_net',ec_type = 'baseline',reduce = mean, alpha = 0.5, truncated = 500)
#write.csv(g2.2[[3]],'~/Documents/SINATRA/Scripts/Data/cariature_v4_elastic_net_mean.csv')
g2.2 = read.csv('~/Documents/SINATRA/Scripts/Data/cariature_v4_elastic_net_mean.csv')[,-1][,-4]
roc_curve_frame2.2 = data.frame(g2.2)
library(ggplot2)
roc_curve_frame2.2 = roc_curve_frame2.2[roc_curve_frame2.2[,3] == 1,]
roc_curve_frame2.2[,3] = rep('Baseline (EN Mean)',dim(roc_curve_frame2.2)[1])
#g2.3 = roc_curve_teeth_baseline(data_dir = base_dir,var_selection = 'elastic_net',ec_type = 'baseline',reduce = max, alpha = 1, truncated = 500)
#write.csv(g2.3[[3]],'~/Documents/SINATRA/Scripts/Data/cariature_v4_ridge_net_max.csv')
g2.3 = read.csv('~/Documents/SINATRA/Scripts/Data/cariature_v4_ridge_net_max.csv')[,-1][,-4]
roc_curve_frame2.3 = data.frame(g2.3)
library(ggplot2)
roc_curve_frame2.3 = roc_curve_frame2.3[roc_curve_frame2.3[,3] == 1,]
roc_curve_frame2.3[,3] = rep('Baseline (RR Max)',dim(roc_curve_frame2.3)[1])
#g2.4 = roc_curve_teeth_baseline(data_dir = base_dir,var_selection = 'elastic_net',ec_type = 'baseline',reduce = mean, alpha = 1, truncated = 500)
#write.csv(g2.4[[3]],'~/Documents/SINATRA/Scripts/Data/cariature_v4_ridge_net_mean.csv')
g2.4 = read.csv('~/Documents/SINATRA/Scripts/Data/cariature_v4_ridge_net_mean.csv')[,-1][,-4]
roc_curve_frame2.4 = data.frame(g2.4)
library(ggplot2)
roc_curve_frame2.4 = roc_curve_frame2.4[roc_curve_frame2.4[,3] == 1,]
roc_curve_frame2.4[,3] = rep('Baseline (RR Mean)',dim(roc_curve_frame2.4)[1])
library(ggplot2)
#### Load in the SINATRA Curves ####
#data_dirs = list.dirs(base_dir,recursive = FALSE)
#roc_curve_frame = rbind(roc_curve_frame2.1,roc_curve_frame2.2,roc_curve_frame2.3,roc_curve_frame2.4)
#rdfmeans = read.csv(paste(data_dirs[1],'/roc_dirs1.csv', sep = ''))[,-4]
#for (k in 2:length(data_dirs)){
#  dir = data_dirs[k]
#  temp = read.csv(paste(dir,'/roc_dirs1.csv', sep = ''))[,-4]
#  rdfmeans = rdfmeans + temp
#}
#rdfmeans = rdfmeans/length(data_dirs)
#rdfmeans = data.frame(rdfmeans)
#rdfmeans = rdfmeans[rdfmeans$X3==35,]
#rdfmeans[,3] = rep('SINATRA',dim(rdfmeans)[1])
#names(rdfmeans) = names(roc_curve_frame)
#roc_curve_frame = rbind(roc_curve_frame,rdfmeans)
#ggplot() + 
#  geom_line(data = subset(roc_curve_frame, X3 == "SINATRA"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5) +
#  geom_line(data = subset(roc_curve_frame, X3 != "SINATRA"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5, linetype = 4) +
#  # geom_line(stat = "identity",aes(color = factor(X3)), linetype = 'dotted') +
#  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", color = "# Cones") +
#  ggtitle(sprintf("10 Causal Regions, 5 Shared Regions, Size 10")) +
#  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
#  coord_equal(ratio=1) +
#  theme_bw() +
#  theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=12),
#        axis.title=element_text(size=16,face="bold")) +
#  scale_colour_hue(l=40)
#ggsave('~/Documents/SINATRA/Scripts/Data/cariature_v4.pdf')


### Limit Shapes ####
path = '~/Documents/new_aligned_shapesv4'
dirs = list.dirs(path, full.names = TRUE, FALSE)
roc_curves2.1 = list()
roc_curves2.2 = list()
truncated = 1000
functions_to_summarize = 10
for (i in 1:length(dirs)){
  dir = dirs[i]
  print(dir)
  total_rate_roc1 = matrix(0, nrow = truncated,ncol = 2)
  total_rate_roc2 = matrix(0, nrow = truncated,ncol = 2)
  dir1 = paste(dir,'/mesh/gp1/', sep = '')
  dir2 = paste(dir,'/mesh/gp2/', sep = '')
  causal_points1_vec = read.csv(paste(dir,'/gp1_spt.csv',sep = ''))
  causal_points2_vec = read.csv(paste(dir,'/gp2_spt.csv',sep = ''))
  
  class_1_true_vertices = which(causal_points1_vec > 0.25)
  class_2_true_vertices = which(causal_points2_vec > 0.25)
  
  class_1_false_vertices = setdiff(1:num_vertices,class_1_true_vertices)
  class_2_false_vertices = setdiff(1:num_vertices,class_2_true_vertices)
  
  limit_shapes_output_files1 = list.files(path = dir1, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs1 = limit_shapes_output_files1[str_detect(limit_shapes_output_files1,'csv')]
  limit_shapes_csvs1 = limit_shapes_csvs1[str_detect(limit_shapes_csvs1,'scrambled')== FALSE]
  
  limit_shapes_output_files2 = list.files(path = dir2, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs2 = limit_shapes_output_files2[str_detect(limit_shapes_output_files2,'csv')]
  limit_shapes_csvs2 = limit_shapes_csvs2[str_detect(limit_shapes_csvs2,'scrambled')== FALSE]
  
  #ROC curve for class 1
  
  counter = 0 
  for (file in limit_shapes_csvs1){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = truncated))){
      
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
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = truncated)) ){
      
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
total_roc1 = rbind(c(0,0), total_roc1)
total_roc2 = rbind(c(0,0), total_roc2)
total_roc1 =  cbind(total_roc1,rep('Limit Shapes Class 1', dim(total_roc1)[1]))
total_roc2 =  cbind(total_roc2,rep('Limit Shapes Class 2', dim(total_roc2)[1]))
roc_curve_frame_limit = as.data.frame(rbind(total_roc1, total_roc2))
roc_curve_frame_limit$V1 = as.numeric(as.character(roc_curve_frame_limit$V1))
roc_curve_frame_limit$V2 = as.numeric(as.character(roc_curve_frame_limit$V2))
write.csv(roc_curve_frame_limit,'~/Documents/SINATRA/Scripts/Data/cariature_v4_limit_shape.csv')
roc_curve_frame_limit = read.csv('~/Documents/SINATRA/Scripts/Data/cariature_v4_limit_shape.csv')[,-1][,-4]
roc_curve_frame_limit = data.frame(roc_curve_frame_limit)
library(ggplot2)
roc_curve_frame_limit = roc_curve_frame_limit[roc_curve_frame_limit[,3] == 'Limit Shapes Class 1',]

save.image('~/Documents/new_aligned_shapesv4/limitshapes.Rdata')
load('~/Dropbox (Princeton)/SINATRA_Data/new_aligned_shapesv4/limitshapes.Rdata')

### Limit Shapes ####
path = '~/Documents/new_aligned_shapesv4'
dirs = list.dirs(path, full.names = TRUE, FALSE)
roc_curves2.1 = list()
roc_curves2.2 = list()
truncated = 500
functions_to_summarize = 10
for (i in 1:length(dirs)){
  dir = dirs[i]
  print(dir)
  total_rate_roc1 = matrix(0, nrow = truncated,ncol = 2)
  total_rate_roc2 = matrix(0, nrow = truncated,ncol = 2)
  dir1 = paste(dir,'/mesh/gp1/', sep = '')
  dir2 = paste(dir,'/mesh/gp2/', sep = '')
  causal_points1_vec = read.csv(paste(dir,'/gp1_spt.csv',sep = ''))
  causal_points2_vec = read.csv(paste(dir,'/gp2_spt.csv',sep = ''))
  
  class_1_true_vertices = which(causal_points1_vec > 0.25)
  class_2_true_vertices = which(causal_points2_vec > 0.25)
  
  class_1_false_vertices = setdiff(1:num_vertices,class_1_true_vertices)
  class_2_false_vertices = setdiff(1:num_vertices,class_2_true_vertices)
  
  limit_shapes_output_files1 = list.files(path = dir1, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs1 = limit_shapes_output_files1[str_detect(limit_shapes_output_files1,'csv')]
  limit_shapes_csvs1 = limit_shapes_csvs1[str_detect(limit_shapes_csvs1,'scrambled')== TRUE]
  
  limit_shapes_output_files2 = list.files(path = dir2, full.names = TRUE,recursive = FALSE)
  limit_shapes_csvs2 = limit_shapes_output_files2[str_detect(limit_shapes_output_files2,'csv')]
  limit_shapes_csvs2 = limit_shapes_csvs2[str_detect(limit_shapes_csvs2,'scrambled')== TRUE]
  
  #ROC curve for class 1
  
  counter = 0 
  for (file in limit_shapes_csvs1){
    rate_ROC <- matrix(0,nrow = 0,ncol = 2)
    g = read.csv(file, header = FALSE)
    g = g[,1:functions_to_summarize]
    g = apply(g, MARGIN = 2, FUN = scale_and_normalize)
    rate_values = apply(g, MARGIN = 1, FUN = max)
    counter = counter + 1
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = truncated))){
      
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
    
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = truncated)) ){
      
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
total_roc1 = rbind(c(0,0), total_roc1)
total_roc2 = rbind(c(0,0), total_roc2)
total_roc1 =  cbind(total_roc1,rep('Limit Shapes Class 1', dim(total_roc1)[1]))
total_roc2 =  cbind(total_roc2,rep('Limit Shapes Class 2', dim(total_roc2)[1]))
roc_curve_frame_limit_scrambled = as.data.frame(rbind(total_roc1, total_roc2))
roc_curve_frame_limit_scrambled$V1 = as.numeric(as.character(roc_curve_frame_limit_scrambled$V1))
roc_curve_frame_limit_scrambled$V2 = as.numeric(as.character(roc_curve_frame_limit_scrambled$V2))
save.image('~/Documents/new_aligned_shapesv4/limitshapes_scrambled.Rdata')
load('~/Dropbox (Princeton)/SINATRA_Data/new_aligned_shapesv4/limitshapes_scrambled.Rdata')
roc_curve_frame = rbind(roc_curve_frame2.1,roc_curve_frame2.2,roc_curve_frame2.3,roc_curve_frame2.4)
roc_curve_frame = rbind(roc_curve_frame,rdfmeans)
roc_curve_frame_limit2 = roc_curve_frame_limit[roc_curve_frame_limit[,3] == 'Limit Shapes Class 1',]
roc_curve_frame_limit2[,3] = 'Limit Shapes'
roc_curve_frame_limit_scrambled2 = roc_curve_frame_limit_scrambled[roc_curve_frame_limit_scrambled[,3] == 'Limit Shapes Class 1',]
roc_curve_frame_limit_scrambled2[,3] = 'Limit Shapes (Misspecified)'
colnames(roc_curve_frame_limit2) = colnames(roc_curve_frame)
colnames(roc_curve_frame_limit_scrambled2) = colnames(roc_curve_frame)
roc_curve_frame = rbind(roc_curve_frame,roc_curve_frame_limit2)
roc_curve_frame = rbind(roc_curve_frame,roc_curve_frame_limit_scrambled2)

ggplot() + 
  geom_line(data = subset(roc_curve_frame, X3 == "SINATRA"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5) +
  geom_line(data = subset(roc_curve_frame, X3 %like% "Limit"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5, linetype = 4) +
#  geom_line(data = subset(roc_curve_frame, X3 %like% "Mean"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
  geom_line(data = subset(roc_curve_frame, X3 %like% "EN"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
  # geom_line(stat = "identity",aes(color = factor(X3)), linetype = 'dotted') +
  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", color = "Method") +
  ggtitle(sprintf("5 Caricatured Peaks")) +
  coord_cartesian(xlim= c(0,1.0))+
#  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
#  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)

ggsave('~/Documents/SINATRA/Scripts/Data/cariature_v4_limit.pdf')
write.csv(roc_curve_frame,file = '~/Documents/SINATRA/Scripts/Data/5peaks_caricature_roc.csv',row.names = FALSE)

#### Test New Function ####
library(glmnet)
data2 = generate_data_sphere_simulation_new(nsim = nsim,dir = dirs, curve_length = len,noise_points = shared_points,
                                       causal_points = causal_points,ball_radius = ball_radius, subdivision = subdivision,
                                       cusps = cusps, causal_regions_1 = causal_regions_1, causal_regions_2 = causal_regions_2,
                                       shared_regions = shared_regions, ec_type = ec_type)
groups = rep(1:(dim(data2$data)[2]/3),each = 3)
#lasso = cv.gglasso(x = data2$data[,-1], y = data2$data[,1], group = groups, loss = "logit")
lasso = cv.glmnet(x = data2$data[,-1], y = data2$data[,1], alpha = 0.0,  family = "binomial")
coefs  = coef(lasso, s = 'lambda.min')

#Binary Option
nonzero = which(abs(coefs[-1])>0.8)
results = unique(groups[nonzero])


# Heatmap Option


g = cbind(groups, abs(coefs[-1]))

library(dplyr)
df = as.data.table(g)

new_df = aggregate(df[,2],list(df$groups),mean)

abs_coefs = new_df$V2
results = groups[which(abs_coefs > 0.5)]
cuts = cut(abs_coefs,50,labels = FALSE)


color1='blue'
color2='lightgreen'
color3='red'
#col_pal = c("#00007F", "blue", "#007FFF", "cyan",  "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
col_pal=c(color1,color1,color2,color2,color2,color3)
col_pal=c(color1,color2,color3)

colfunc <- colorRampPalette(col_pal)
heat_colors=colfunc(max(cuts) - min(cuts))[cuts - min(cuts)]

mesh = vcgSphere(subdivision = 3)
mesh$vb[1:3,] = t(data2$complex_points[[1]])
cols = rep('white', dim(mesh$vb)[2])
cols[data$causal_points1] = 'red'
cols[data$causal_points2] = 'red'
cols[data$noise] = 'blue'
cols[results] = 'green'

mfrow3d(1,2)
plot3d(mesh, col = cols)
plot3d(mesh, col = heat_colors)
#### Sanity Checks ####

j = cor(t(data_summary$data[,-1]))
heatmap(j)


generate_ROC_with_coned_directions(nsim = 10, curve_length = 50, grid_size = 25, distance_to_causal_point = 0.1, 
                                   causal_points = causal_points,shared_points = shared_points,  eta = 0.1, 
                                   truncated = -1, two_curves = TRUE, ball = TRUE, ball_radius = 2.5, type = 'vertex',
                                   min_points = 3,directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                   mode = 'sphere_baseline', 
                                   subdivision = 3,num_causal_region = 2, num_shared_region = 1)


#### Comp 2 ####
indices = c(261,295,833,224,939,278,252,593,177,357,4845,4215,4619,162,244,4618,4225,31,4705,4524,219,4301,4898,262,4860,69,64,4478,4452,4224,102,179,4176)

color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)

col_pal2=c(color1,color1,color2,color2,color3)
colfunc2 <- colorRampPalette(col_pal2)

