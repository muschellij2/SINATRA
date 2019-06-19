setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
library(R.utils)
sourceDirectory('~/Documents/SINATRA/SINATRA_Pipeline_Branch/')
sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')
#load('Caricature_ROCV2.Rdata')
#library(svd)
#set.seed(15)
#set.seed(35)
set.seed(15)

#Parameters for the Analysis
ball = TRUE
ball_radius = 0.5
ec_type = 'ECT'
rotation_matrix=matrix(c(0.99972576,0.02127766,0.00978078,0,0.01589701,-0.92330807,0.38373107,0,0.01719557,-0.38347024,-0.92339307,0,0,0,0,1),ncol=4,byrow=TRUE)
#Specifying Directories and the Associated Files

#### 3 v 3 ####
base_dir = '~/Documents/new_aligned_shapesv3/'
data_dirs = list.dirs(base_dir,recursive = FALSE)
#data_dirs = data_dirs[-(1:47)]
#data_dirs = data_dirs[1]
for (dir in data_dirs){
  old_data_dir = paste(dir,'/mesh/gp1',sep='')
  new_data_dir = paste(dir,'/mesh/gp2',sep='')
  old_data_files = list.files(old_data_dir, full.names = TRUE)
  new_data_files = list.files(new_data_dir, full.names = TRUE)
  class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE)
  class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)
  roc_curves = list()
  total_dirs = c(5,15,25,35)
  #total_dirs = c(1)
  for (k in 1:length(total_dirs)){
    dirs = total_dirs[k]
    print(paste("on dir", dirs))
    pset = list(num_cones = dirs, cap_radius = 0.15, len = 100, directions_per_cone = 5,
                      dirs = generate_equidistributed_cones(num_directions = dirs, cap_radius =  0.15, directions_per_cone = 5))
    data_summary=real_data_summary(dir1=new_data_dir,dir2 = old_data_dir,direction=pset$dirs,class1='Vegetable', class2='Vegetable',
                                      radius=0,accuracy=FALSE,len = pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
    roc_curve = compute_roc_curve_teeth(data_dir1 = old_data_dir, data_dir2 = new_data_dir, gamma = 0.25,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                 rate_values = data_summary$Rate2[,2],directions_per_cone = pset$directions_per_cone,curve_length = pset$len,
                                 directions = pset$dirs,truncated = 500,ball_radius = ball_radius, ball = ball, radius = 1, two_curves = TRUE)
    #print(roc_curve)
    roc_curve[,3] = dirs
    roc_frame = data.frame(roc_curve)
    roc_curves[[k]] = roc_frame
  }
  roc_curves_total1 = roc_curves[[1]][1:500,]
  for (k in 2:length(total_dirs)){
    roc_curves_total1 = rbind(roc_curves_total1,roc_curves[[k]][1:500,])
  }
  write.csv(roc_curves_total1,file = paste(dir,'/roc_dirs1.csv',sep=''),row.names = FALSE)
  roc_curves_total2 = roc_curves[[1]][501:1000,]
  for (k in 2:length(total_dirs)){
    roc_curves_total2 = rbind(roc_curves_total2,roc_curves[[k]][501:1000,])
  }
  write.csv(roc_curves_total2,file = paste(dir,'/roc_dirs2.csv',sep=''),row.names = FALSE)
  save.image(paste(dir,'/rocs_and_rate.R',sep=''))
}


# Making the ROC Curves
base_dir = '~/Documents/new_aligned_shapesv3/'
data_dirs = list.dirs(base_dir,recursive = FALSE)

roc_curve1 = read.csv('~/Documents/new_aligned_shapesv3/V13_3peak_2gp_v1_bw/roc_dirs1.csv')
roc_curve2 = read.csv('~/Documents/new_aligned_shapesv3/V13_3peak_2gp_v1_bw/roc_dirs2.csv')

data_dirs2 = data_dirs[-1]

for (dir in data_dirs2){
  roc_curve1.5 = read.csv(paste(dir,'/roc_dirs1.csv',sep=''))
  roc_curve2.5 = read.csv(paste(dir,'/roc_dirs2.csv',sep=''))
  roc_curve1 = roc_curve1 + roc_curve1.5
  roc_curve2 = roc_curve2 + roc_curve2.5
}
roc_curve1 = roc_curve1/50
roc_curve2 = roc_curve2/50

roc_curve1_frame = data.frame(roc_curve1)
roc_curve2_frame = data.frame(roc_curve2)

ROC_curve_plt <- ggplot(data <- roc_curve1_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Cones") +
  ggtitle(sprintf("ROC Curve for Class 1 Teeth")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt)
ggsave('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/3_peaks/Class_1_roc_dirs.pdf')
ROC_curve_plt2 <- ggplot(data <- roc_curve2_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Cones") +
  ggtitle(sprintf("ROC Curve for Class 2 Teeth")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt2)
ggsave('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/3_peaks/Class_2_roc_dirs.pdf')

# Plotting Now
dir = data_dirs[10]
load(paste(dir,'/rocs_and_rate.R',sep=''))
old_data_dir = paste(dir,'/mesh/gp1',sep='')
new_data_dir = paste(dir,'/mesh/gp2',sep='')
old_data_files = list.files(old_data_dir, full.names = TRUE)
new_data_files = list.files(new_data_dir, full.names = TRUE)
class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE)
class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)
ind = 3
old_veg1 = vcgImport(old_data_files[ind])
new_veg1 = vcgImport(new_data_files[ind])
og_mesh = vcgImport(file = 'Data/new_aligned_shapesv3/Old_Follivore/clean_V13_sas.off')
new_veg_1 = process_off_file_v3(new_data_files[ind])
old_veg_1 = process_off_file_v3(old_data_files[ind])
v1 = t(new_veg1$vb[-4,])
v2 = t(old_veg1$vb[-4,])


new_heat =  reconstruct_vertices_on_shape(dir = pset$dirs,complex = new_veg_1,rate_vals = data_summary$Rate2[,2],
                                              len = pset$len,cuts = 500,cone_size = pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
old_heat =  reconstruct_vertices_on_shape(dir = pset$dirs,complex = old_veg_1,rate_vals = data_summary$Rate2[,2],
                                              len = pset$len,cuts = 500,cone_size = pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)
col_pal = c("#00007F", "blue", "#007FFF", "cyan",  "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")

colfunc <- colorRampPalette(col_pal)
new_heat_colors=colfunc(max(new_heat[,1]) - min(new_heat[,1]))[new_heat[,1] - min(new_heat[,1])]
old_heat_colors=colfunc(max(old_heat[,1]) - min(old_heat[,1]))[old_heat[,1] - min(old_heat[,1])]
mfrow3d(2,3)
new_veg_colors = rep('white', dim(new_veg1$vb)[2])
plot3d(og_mesh,col = new_veg_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)
new_veg_colors[(which(class_2_probs > 0.25))] = 'red'
plot3d(new_veg1,col = new_veg_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(new_veg1,col = new_heat_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)

old_veg_colors = rep('white', dim(old_veg1$vb)[2])
plot3d(og_mesh,col = old_veg_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)
old_veg_colors[(which(class_1_probs > 0.25))] = 'red'

plot3d(old_veg1,col = old_veg_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)
plot3d(old_veg1,col = old_heat_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)

col.table=colfunc(1 + max(new_heat[,1]) - min(new_heat[,1]))
color.bar(col.table,min=10000 * round(min(new_heat[,2]) - min(new_heat[,2]),4),max=10000 * round(max(new_heat[,2])-min(new_heat[,2]),4))
# Varying Length Up

base_dir = '~/Documents/new_aligned_shapesv3/'
data_dirs = list.dirs(base_dir,recursive = FALSE)
data_dirs = data_dirs[-(1:47)]
#omit 6,47
for (dir in data_dirs){
  old_data_dir = paste(dir,'/mesh/gp1',sep='')
  new_data_dir = paste(dir,'/mesh/gp2',sep='')
  old_data_files = list.files(old_data_dir, full.names = TRUE)
  new_data_files = list.files(new_data_dir, full.names = TRUE)
  class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE)
  class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)
  roc_curves = list()
  dirs = 15
  total_lens = c(10,25,50,75)
  #total_dirs = c(1)
  for (k in 1:length(total_lens)){
    len = total_lens[k]
    print(paste("on length", len, 'dir', dirs))
    pset = list(num_cones = dirs, cap_radius = 0.15, len = len, directions_per_cone = 5,
                      dirs = generate_equidistributed_cones(num_directions = dirs, cap_radius =  0.15, directions_per_cone = 5))
    data_summary=real_data_summary(dir1=new_data_dir,dir2 = old_data_dir,direction=pset$dirs,class1='Vegetable', class2='Vegetable',
                                      radius=0,accuracy=FALSE,len = pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
    roc_curve = compute_roc_curve_teeth(data_dir1 = old_data_dir, data_dir2 = new_data_dir, gamma = 0.25,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                 rate_values = data_summary$Rate2[,2],directions_per_cone = pset$directions_per_cone,curve_length = pset$len,
                                 directions = pset$dirs,truncated = 500,ball_radius = ball_radius, ball = ball, radius = 1, two_curves = TRUE)
    #print(roc_curve)
    roc_curve[,3] = len
    roc_frame = data.frame(roc_curve)
    roc_curves[[k]] = roc_frame
  }
  roc_curves_total1 = roc_curves[[1]][1:500,]
  for (k in 2:length(total_lens)){
    roc_curves_total1 = rbind(roc_curves_total1,roc_curves[[k]][1:500,])
  }
  write.csv(roc_curves_total1,file = paste(dir,'/roc_lens1.csv',sep=''),row.names = FALSE)
  roc_curves_total2 = roc_curves[[1]][501:1000,]
  for (k in 2:length(total_lens)){
    roc_curves_total2 = rbind(roc_curves_total2,roc_curves[[k]][501:1000,])
  }
  write.csv(roc_curves_total2,file = paste(dir,'/roc_lens2.csv',sep=''),row.names = FALSE)
  save.image(paste(dir,'/rocs_and_rate_lens.R',sep=''))
}


# Making the ROC Curves
base_dir = '~/Documents/new_aligned_shapesv3/'
data_dirs = list.dirs(base_dir,recursive = FALSE)

roc_curve_len_1 = read.csv('~/Documents/new_aligned_shapesv3/V13_3peak_2gp_v1_bw/roc_lens1.csv')
roc_curve_len_2 = read.csv('~/Documents/new_aligned_shapesv3/V13_3peak_2gp_v1_bw/roc_lens2.csv')

data_dirs2 = data_dirs[-c(1,6,47)]

for (dir in data_dirs2){
  roc_curve_len_1.5 = read.csv(paste(dir,'/roc_lens1.csv',sep=''))
  roc_curve_len_2.5 = read.csv(paste(dir,'/roc_lens2.csv',sep=''))
  roc_curve_len_1 = roc_curve_len_1 + roc_curve_len_1.5
  roc_curve_len_2 = roc_curve_len_2 + roc_curve_len_2.5
}
roc_curve_len_1 = roc_curve_len_1/48
roc_curve_len_2 = roc_curve_len_2/48


roc_curve_1_100 = roc_curve1[roc_curve1[,3] == 15,]
roc_curve_2_100 = roc_curve2[roc_curve2[,3] == 15,]
roc_curve_1_100[,3] = 100
roc_curve_2_100[,3] = 100


roc_curve_len_1 = rbind(roc_curve_len_1,roc_curve_1_100)
roc_curve_len_2 = rbind(roc_curve_len_2,roc_curve_2_100)

roc_curve_len_1_frame = data.frame(roc_curve_len_1)
roc_curve_len_2_frame = data.frame(roc_curve_len_2)

ROC_curve_plt <- ggplot(data <- roc_curve_len_1_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Lengths") +
  ggtitle(sprintf("ROC Curve for Class 1 Teeth - 3 Peaks - Varying Length")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt)
ggsave('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/3_peaks/Class_1_roc_lengths.pdf')
ROC_curve_plt2 <- ggplot(data <- roc_curve_len_2_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Lengths") +
  ggtitle(sprintf("ROC Curve for Class 2 Teeth - 3 Peaks - Varying Length")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt2)
ggsave('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/3_peaks/Class_2_roc_lengths.pdf')
# Varying Angle Up

base_dir = '~/Documents/new_aligned_shapesv3/'
data_dirs = list.dirs(base_dir,recursive = FALSE)
data_dirs = data_dirs[-(1:32)]
#omit 6,47
for (dir in data_dirs){
  old_data_dir = paste(dir,'/mesh/gp1',sep='')
  new_data_dir = paste(dir,'/mesh/gp2',sep='')
  old_data_files = list.files(old_data_dir, full.names = TRUE)
  new_data_files = list.files(new_data_dir, full.names = TRUE)
  class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE)
  class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)
  roc_curves = list()
  dirs = 15
  len = 50
  total_angles = c(0.05,0.15,0.25,0.35,0.45)
  #total_dirs = c(1)
  for (k in 1:length(total_angles)){
    angle = total_angles[k]
    print(paste("on length", len, 'dir', dirs, 'angle', angle))
    pset = list(num_cones = dirs, cap_radius = angle, len = len, directions_per_cone = 5,
                      dirs = generate_equidistributed_cones(num_directions = dirs, cap_radius =  angle, directions_per_cone = 5))
    data_summary=real_data_summary(dir1=new_data_dir,dir2 = old_data_dir,direction=pset$dirs,class1='Vegetable', class2='Vegetable',
                                      radius=0,accuracy=FALSE,len = pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
    roc_curve = compute_roc_curve_teeth(data_dir1 = old_data_dir, data_dir2 = new_data_dir, gamma = 0.25,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                 rate_values = data_summary$Rate2[,2],directions_per_cone = pset$directions_per_cone,curve_length = pset$len,
                                 directions = pset$dirs,truncated = 500,ball_radius = ball_radius, ball = ball, radius = 1, two_curves = TRUE)
    #print(roc_curve)
    roc_curve[,3] = angle
    roc_frame = data.frame(roc_curve)
    roc_curves[[k]] = roc_frame
  }
  roc_curves_total1 = roc_curves[[1]][1:500,]
  for (k in 2:length(total_angles)){
    roc_curves_total1 = rbind(roc_curves_total1,roc_curves[[k]][1:500,])
  }
  write.csv(roc_curves_total1,file = paste(dir,'/roc_angles1.csv',sep=''),row.names = FALSE)
  roc_curves_total2 = roc_curves[[1]][501:1000,]
  for (k in 2:length(total_angles)){
    roc_curves_total2 = rbind(roc_curves_total2,roc_curves[[k]][501:1000,])
  }
  write.csv(roc_curves_total2,file = paste(dir,'/roc_angles2.csv',sep=''),row.names = FALSE)
  save.image(paste(dir,'/rocs_and_rate_angles.R',sep=''))
}
### ROC Curve for angle ###
base_dir = '~/Documents/new_aligned_shapesv3/'
data_dirs = list.dirs(base_dir,recursive = FALSE)

roc_curve1 = read.csv('~/Documents/new_aligned_shapesv3/V13_3peak_2gp_v1_bw/roc_angles1.csv')
roc_curve2 = read.csv('~/Documents/new_aligned_shapesv3/V13_3peak_2gp_v1_bw/roc_angles2.csv')

data_dirs2 = data_dirs[-1]

for (dir in data_dirs2){
  roc_curve1.5 = read.csv(paste(dir,'/roc_angles1.csv',sep=''))
  roc_curve2.5 = read.csv(paste(dir,'/roc_angles2.csv',sep=''))
  roc_curve1 = roc_curve1 + roc_curve1.5
  roc_curve2 = roc_curve2 + roc_curve2.5
}
roc_curve1 = roc_curve1/50
roc_curve2 = roc_curve2/50

roc_curve1_frame = data.frame(roc_curve1)
roc_curve2_frame = data.frame(roc_curve2)

ROC_curve_plt <- ggplot(data <- roc_curve1_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Cone Angle") +
  ggtitle(sprintf("ROC Curve for Class 1 Teeth - 3 Peaks")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt)
ggsave('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/3_peaks/Class_1_roc_angles.pdf')
ROC_curve_plt2 <- ggplot(data <- roc_curve2_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Cone Angle") +
  ggtitle(sprintf("ROC Curve for Class 2 Teeth - 3 Peaks")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt2)
ggsave('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/3_peaks/Class_2_roc_angles.pdf')



#dirpercone
base_dir = '~/Documents/new_aligned_shapesv3/'
data_dirs = list.dirs(base_dir,recursive = FALSE)
data_dirs = data_dirs[-(1:42)]
for (dir in data_dirs){
  old_data_dir = paste(dir,'/mesh/gp1',sep='')
  new_data_dir = paste(dir,'/mesh/gp2',sep='')
  old_data_files = list.files(old_data_dir, full.names = TRUE)
  new_data_files = list.files(new_data_dir, full.names = TRUE)
  class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE)
  class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)
  roc_curves = list()
  dirs = 15
  angle = 0.15
  len = 50
  dpc = c(1,3,5,7,9)
  #total_dirs = c(1)
  for (k in 1:length(dpc)){
    dirpercone = dpc[k]
    print(paste("on length", len, 'dir', dirs, 'angle', angle,'dirpercone',dirpercone))
    pset = list(num_cones = dirs, cap_radius = angle, len = len, directions_per_cone = dirpercone,
                      dirs = generate_equidistributed_cones(num_directions = dirs, cap_radius =  angle, directions_per_cone = dirpercone))
    data_summary=real_data_summary(dir1=new_data_dir,dir2 = old_data_dir,direction=pset$dirs,class1='Vegetable', class2='Vegetable',
                                      radius=0,accuracy=FALSE,len = pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
    roc_curve = compute_roc_curve_teeth(data_dir1 = old_data_dir, data_dir2 = new_data_dir, gamma = 0.25,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                 rate_values = data_summary$Rate2[,2],directions_per_cone = pset$directions_per_cone,curve_length = pset$len,
                                 directions = pset$dirs,truncated = 500,ball_radius = ball_radius, ball = ball, radius = 1, two_curves = TRUE)
    #print(roc_curve)
    roc_curve[,3] = dirpercone
    roc_frame = data.frame(roc_curve)
    roc_curves[[k]] = roc_frame
  }
  roc_curves_total1 = roc_curves[[1]][1:500,]
  for (k in 2:length(dpc)){
    roc_curves_total1 = rbind(roc_curves_total1,roc_curves[[k]][1:500,])
  }
  write.csv(roc_curves_total1,file = paste(dir,'/roc_dpc1.csv',sep=''),row.names = FALSE)
  roc_curves_total2 = roc_curves[[1]][501:1000,]
  for (k in 2:length(dpc)){
    roc_curves_total2 = rbind(roc_curves_total2,roc_curves[[k]][501:1000,])
  }
  write.csv(roc_curves_total2,file = paste(dir,'/roc_dpc2.csv',sep=''),row.names = FALSE)
  save.image(paste(dir,'/rocs_and_rate_dpc.R',sep=''))
}

#### Making ROC curves ####

base_dir = '~/Documents/new_aligned_shapesv3/'
data_dirs = list.dirs(base_dir,recursive = FALSE)
data_dirs = data_dirs[-c(42,33,34,35,25,20)]

roc_curve1 = read.csv('~/Documents/new_aligned_shapesv3/V13_3peak_2gp_v1_bw/roc_dpc1.csv')
roc_curve2 = read.csv('~/Documents/new_aligned_shapesv3/V13_3peak_2gp_v1_bw/roc_dpc2.csv')

data_dirs2 = data_dirs[-1]

for (dir in data_dirs2){
  roc_curve1.5 = read.csv(paste(dir,'/roc_dpc1.csv',sep=''))
  roc_curve2.5 = read.csv(paste(dir,'/roc_dpc2.csv',sep=''))
  roc_curve1 = roc_curve1 + roc_curve1.5
  roc_curve2 = roc_curve2 + roc_curve2.5
}
roc_curve1 = roc_curve1/44
roc_curve2 = roc_curve2/44

roc_curve1_frame = data.frame(roc_curve1)
roc_curve2_frame = data.frame(roc_curve2)

ROC_curve_plt <- ggplot(data <- roc_curve1_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Directions") +
  ggtitle(sprintf("ROC Curve for Class 1 Teeth - 3 Peaks")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt)
ggsave('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/3_peaks/Class_1_roc_dir_per_cone.pdf')
ROC_curve_plt2 <- ggplot(data <- roc_curve2_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Directions") +
  ggtitle(sprintf("ROC Curve for Class 2 Teeth - 3 Peaks")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt2)
ggsave('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/3_peaks/Class_2_roc_dir_per_cone.pdf')


#### 5 v 5 ####
base_dir = '~/Documents/new_aligned_shapesv4/'
data_dirs = list.dirs(base_dir,recursive = FALSE)
#data_dirs = data_dirs[(12:21)]
for (dir in data_dirs){
  old_data_dir = paste(dir,'/mesh/gp1',sep='')
  new_data_dir = paste(dir,'/mesh/gp2',sep='')
  old_data_files = list.files(old_data_dir, full.names = TRUE)
  new_data_files = list.files(new_data_dir, full.names = TRUE)
  class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE)
  class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)
  roc_curves = list()
  total_dirs = c(5,15,25,35)
  #total_dirs = c(1)
  for (k in 1:length(total_dirs)){
    dirs = total_dirs[k]
    print(paste("on dir", dirs))
    pset = list(num_cones = dirs, cap_radius = 0.15, len = 50, directions_per_cone = 5,
                      dirs = generate_equidistributed_cones(num_directions = dirs, cap_radius =  0.15, directions_per_cone = 5))
    data_summary=real_data_summary(dir1=new_data_dir,dir2 = old_data_dir,direction=pset$dirs,class1='Vegetable', class2='Vegetable',
                                      radius=0,accuracy=FALSE,len = pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
    roc_curve = compute_roc_curve_teeth(data_dir1 = old_data_dir, data_dir2 = new_data_dir, gamma = 0.25,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                 rate_values = data_summary$Rate2[,2],directions_per_cone = pset$directions_per_cone,curve_length = pset$len,
                                 directions = pset$dirs,truncated = 500,ball_radius = ball_radius, ball = ball, radius = 1, two_curves = TRUE)
    #print(roc_curve)
    roc_curve[,3] = dirs
    roc_frame = data.frame(roc_curve)
    roc_curves[[k]] = roc_frame
  }
  roc_curves_total1 = roc_curves[[1]][1:500,]
  for (k in 2:length(total_dirs)){
    roc_curves_total1 = rbind(roc_curves_total1,roc_curves[[k]][1:500,])
  }
  write.csv(roc_curves_total1,file = paste(dir,'/roc_dirs1.csv',sep=''),row.names = FALSE)
  roc_curves_total2 = roc_curves[[1]][501:1000,]
  for (k in 2:length(total_dirs)){
    roc_curves_total2 = rbind(roc_curves_total2,roc_curves[[k]][501:1000,])
  }
  write.csv(roc_curves_total2,file = paste(dir,'/roc_dirs2.csv',sep=''),row.names = FALSE)
  save.image(paste(dir,'/rocs_and_rate.R',sep=''))
}
# Making the ROC Curves
base_dir = '~/Documents/new_aligned_shapesv4/'
data_dirs = list.dirs(base_dir,recursive = FALSE)

roc_curve1 = read.csv('~/Documents/new_aligned_shapesv4/V13_5peak_2gp_v1_bw/roc_dirs1.csv')
roc_curve2 = read.csv('~/Documents/new_aligned_shapesv4/V13_5peak_2gp_v1_bw/roc_dirs2.csv')

data_dirs2 = data_dirs[-1]

for (dir in data_dirs2){
  roc_curve1.5 = read.csv(paste(dir,'/roc_dirs1.csv',sep=''))
  roc_curve2.5 = read.csv(paste(dir,'/roc_dirs2.csv',sep=''))
  roc_curve1 = roc_curve1 + roc_curve1.5
  roc_curve2 = roc_curve2 + roc_curve2.5
}
roc_curve1 = roc_curve1/50
roc_curve2 = roc_curve2/50

roc_curve1_frame = data.frame(roc_curve1)
roc_curve2_frame = data.frame(roc_curve2)

ROC_curve_plt <- ggplot(data <- roc_curve1_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Cones") +
  ggtitle(sprintf("ROC Curve for Class 1 Teeth - 5 Peaks")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt)
ggsave('~/Dropbox/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/5_peaks/Class_1_roc.pdf')
ROC_curve_plt2 <- ggplot(data <- roc_curve2_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Cones") +
  ggtitle(sprintf("ROC Curve for Class 2 Teeth - 5 Peaks")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt2)
ggsave('~/Dropbox/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/5_peaks/Class_2_roc.pdf')

# Varying Length Up

base_dir = '~/Documents/new_aligned_shapesv4/'
data_dirs = list.dirs(base_dir,recursive = FALSE)
data_dirs = data_dirs[-(1:25)]
for (dir in data_dirs){
  old_data_dir = paste(dir,'/mesh/gp1',sep='')
  new_data_dir = paste(dir,'/mesh/gp2',sep='')
  old_data_files = list.files(old_data_dir, full.names = TRUE)
  new_data_files = list.files(new_data_dir, full.names = TRUE)
  class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE)
  class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)
  roc_curves = list()
  dirs = 15
  total_lens = c(10,25,50,75)
  #total_dirs = c(1)
  for (k in 1:length(total_lens)){
    len = total_lens[k]
    print(paste("on length", len, 'dir', dirs))
    pset = list(num_cones = dirs, cap_radius = 0.15, len = len, directions_per_cone = 5,
                      dirs = generate_equidistributed_cones(num_directions = dirs, cap_radius =  0.15, directions_per_cone = 5))
    data_summary=real_data_summary(dir1=new_data_dir,dir2 = old_data_dir,direction=pset$dirs,class1='Vegetable', class2='Vegetable',
                                      radius=0,accuracy=FALSE,len = pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
    roc_curve = compute_roc_curve_teeth(data_dir1 = old_data_dir, data_dir2 = new_data_dir, gamma = 0.25,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                 rate_values = data_summary$Rate2[,2],directions_per_cone = pset$directions_per_cone,curve_length = pset$len,
                                 directions = pset$dirs,truncated = 500,ball_radius = ball_radius, ball = ball, radius = 1, two_curves = TRUE)
    #print(roc_curve)
    roc_curve[,3] = len
    roc_frame = data.frame(roc_curve)
    roc_curves[[k]] = roc_frame
  }
  roc_curves_total1 = roc_curves[[1]][1:500,]
  for (k in 2:length(total_lens)){
    roc_curves_total1 = rbind(roc_curves_total1,roc_curves[[k]][1:500,])
  }
  write.csv(roc_curves_total1,file = paste(dir,'/roc_lens1.csv',sep=''),row.names = FALSE)
  roc_curves_total2 = roc_curves[[1]][501:1000,]
  for (k in 2:length(total_lens)){
    roc_curves_total2 = rbind(roc_curves_total2,roc_curves[[k]][501:1000,])
  }
  write.csv(roc_curves_total2,file = paste(dir,'/roc_lens2.csv',sep=''),row.names = FALSE)
  save.image(paste(dir,'/rocs_and_rate_lens.R',sep=''))
}


# Making the ROC Curves
base_dir = '~/Documents/new_aligned_shapesv4/'
data_dirs = list.dirs(base_dir,recursive = FALSE)

roc_curve_len_1 = read.csv('~/Documents/new_aligned_shapesv4/V13_5peak_2gp_v1_bw/roc_lens1.csv')
roc_curve_len_2 = read.csv('~/Documents/new_aligned_shapesv4/V13_5peak_2gp_v1_bw/roc_lens2.csv')

data_dirs = data_dirs[-1]

for (dir in data_dirs){
  roc_curve_len_1.5 = read.csv(paste(dir,'/roc_lens1.csv',sep=''))
  roc_curve_len_2.5 = read.csv(paste(dir,'/roc_lens2.csv',sep=''))
  roc_curve_len_1 = roc_curve_len_1 + roc_curve_len_1.5
  roc_curve_len_2 = roc_curve_len_2 + roc_curve_len_2.5
}
roc_curve_len_1 = roc_curve_len_1/50
roc_curve_len_2 = roc_curve_len_2/50


roc_curve_1_100 = roc_curve1[roc_curve1[,3] == 15,]
roc_curve_2_100 = roc_curve2[roc_curve2[,3] == 15,]
roc_curve_1_100[,3] = 100
roc_curve_2_100[,3] = 100


roc_curve_len_1 = rbind(roc_curve_len_1,roc_curve_1_100)
roc_curve_len_2 = rbind(roc_curve_len_2,roc_curve_2_100)

roc_curve_len_1_frame = data.frame(roc_curve_len_1)
roc_curve_len_2_frame = data.frame(roc_curve_len_2)

ROC_curve_plt <- ggplot(data <- roc_curve_len_1_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Lengths") +
  ggtitle(sprintf("ROC Curve for Class 1 Teeth - 5 Peaks - Varying Length")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt)
ggsave('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/5_peaks/Class_1_roc_lengths.pdf')
ROC_curve_plt2 <- ggplot(data <- roc_curve_len_2_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Lengths") +
  ggtitle(sprintf("ROC Curve for Class 2 Teeth - 5 Peaks - Varying Length")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt2)
ggsave('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/5_peaks/Class_2_roc_lengths.pdf')

#Angles
base_dir = '~/Documents/new_aligned_shapesv4/'
data_dirs = list.dirs(base_dir,recursive = FALSE)
for (dir in data_dirs){
  old_data_dir = paste(dir,'/mesh/gp1',sep='')
  new_data_dir = paste(dir,'/mesh/gp2',sep='')
  old_data_files = list.files(old_data_dir, full.names = TRUE)
  new_data_files = list.files(new_data_dir, full.names = TRUE)
  class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE)
  class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)
  roc_curves = list()
  dirs = 15
  len = 50
  total_angles = c(0.05,0.15,0.25,0.35,0.45)
  #total_dirs = c(1)
  for (k in 1:length(total_angles)){
    angle = total_angles[k]
    print(paste("on length", len, 'dir', dirs, 'angle', angle))
    pset = list(num_cones = dirs, cap_radius = angle, len = len, directions_per_cone = 5,
                      dirs = generate_equidistributed_cones(num_directions = dirs, cap_radius =  angle, directions_per_cone = 5))
    data_summary=real_data_summary(dir1=new_data_dir,dir2 = old_data_dir,direction=pset$dirs,class1='Vegetable', class2='Vegetable',
                                      radius=0,accuracy=FALSE,len = pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
    roc_curve = compute_roc_curve_teeth(data_dir1 = old_data_dir, data_dir2 = new_data_dir, gamma = 0.25,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                 rate_values = data_summary$Rate2[,2],directions_per_cone = pset$directions_per_cone,curve_length = pset$len,
                                 directions = pset$dirs,truncated = 500,ball_radius = ball_radius, ball = ball, radius = 1, two_curves = TRUE)
    #print(roc_curve)
    roc_curve[,3] = angle
    roc_frame = data.frame(roc_curve)
    roc_curves[[k]] = roc_frame
  }
  roc_curves_total1 = roc_curves[[1]][1:500,]
  for (k in 2:length(total_angles)){
    roc_curves_total1 = rbind(roc_curves_total1,roc_curves[[k]][1:500,])
  }
  write.csv(roc_curves_total1,file = paste(dir,'/roc_angles1.csv',sep=''),row.names = FALSE)
  roc_curves_total2 = roc_curves[[1]][501:1000,]
  for (k in 2:length(total_angles)){
    roc_curves_total2 = rbind(roc_curves_total2,roc_curves[[k]][501:1000,])
  }
  write.csv(roc_curves_total2,file = paste(dir,'/roc_angles2.csv',sep=''),row.names = FALSE)
  save.image(paste(dir,'/rocs_and_rate_angles.R',sep=''))
}

### ROC Curve for angle ###
base_dir = '~/Documents/new_aligned_shapesv4/'
data_dirs = list.dirs(base_dir,recursive = FALSE)

roc_curve1 = read.csv('~/Documents/new_aligned_shapesv4/V13_5peak_2gp_v1_bw/roc_angles1.csv')
roc_curve2 = read.csv('~/Documents/new_aligned_shapesv4/V13_5peak_2gp_v1_bw/roc_angles2.csv')

data_dirs2 = data_dirs[-1]

for (dir in data_dirs2){
  roc_curve1.5 = read.csv(paste(dir,'/roc_angles1.csv',sep=''))
  roc_curve2.5 = read.csv(paste(dir,'/roc_angles2.csv',sep=''))
  roc_curve1 = roc_curve1 + roc_curve1.5
  roc_curve2 = roc_curve2 + roc_curve2.5
}
roc_curve1 = roc_curve1/50
roc_curve2 = roc_curve2/50

roc_curve1_frame = data.frame(roc_curve1)
roc_curve2_frame = data.frame(roc_curve2)

ROC_curve_plt <- ggplot(data <- roc_curve1_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Cone Angle") +
  ggtitle(sprintf("ROC Curve for Class 1 Teeth - 5 Peaks")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt)
ggsave('~/Dropbox/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/5_peaks/Class_1_roc_angles.pdf')
ROC_curve_plt2 <- ggplot(data <- roc_curve2_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Cone Angle") +
  ggtitle(sprintf("ROC Curve for Class 2 Teeth - 5 Peaks")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt2)
ggsave('~/Dropbox/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/5_peaks/Class_2_roc_angles.pdf')


#dirpercone
base_dir = '~/Documents/new_aligned_shapesv4/'
data_dirs = list.dirs(base_dir,recursive = FALSE)
for (dir in data_dirs){
  old_data_dir = paste(dir,'/mesh/gp1',sep='')
  new_data_dir = paste(dir,'/mesh/gp2',sep='')
  old_data_files = list.files(old_data_dir, full.names = TRUE)
  new_data_files = list.files(new_data_dir, full.names = TRUE)
  class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE)
  class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)
  roc_curves = list()
  dirs = 15
  angle = 0.15
  len = 50
  dpc = c(1,3,5,7,9)
  #total_dirs = c(1)
  for (k in 1:length(dpc)){
    dirpercone = dpc[k]
    print(paste("on length", len, 'dir', dirs, 'angle', angle,'dirpercone',dirpercone))
    pset = list(num_cones = dirs, cap_radius = angle, len = len, directions_per_cone = dirpercone,
                      dirs = generate_equidistributed_cones(num_directions = dirs, cap_radius =  angle, directions_per_cone = dirpercone))
    data_summary=real_data_summary(dir1=new_data_dir,dir2 = old_data_dir,direction=pset$dirs,class1='Vegetable', class2='Vegetable',
                                      radius=0,accuracy=FALSE,len = pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
    roc_curve = compute_roc_curve_teeth(data_dir1 = old_data_dir, data_dir2 = new_data_dir, gamma = 0.25,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                 rate_values = data_summary$Rate2[,2],directions_per_cone = pset$directions_per_cone,curve_length = pset$len,
                                 directions = pset$dirs,truncated = 500,ball_radius = ball_radius, ball = ball, radius = 1, two_curves = TRUE)
    #print(roc_curve)
    roc_curve[,3] = dirpercone
    roc_frame = data.frame(roc_curve)
    roc_curves[[k]] = roc_frame
  }
  roc_curves_total1 = roc_curves[[1]][1:500,]
  for (k in 2:length(dpc)){
    roc_curves_total1 = rbind(roc_curves_total1,roc_curves[[k]][1:500,])
  }
  write.csv(roc_curves_total1,file = paste(dir,'/roc_dpc1.csv',sep=''),row.names = FALSE)
  roc_curves_total2 = roc_curves[[1]][501:1000,]
  for (k in 2:length(dpc)){
    roc_curves_total2 = rbind(roc_curves_total2,roc_curves[[k]][501:1000,])
  }
  write.csv(roc_curves_total2,file = paste(dir,'/roc_dpc2.csv',sep=''),row.names = FALSE)
  save.image(paste(dir,'/rocs_and_rate_dpc.R',sep=''))
}
### ROC Curve for dpc ###
base_dir = '~/Documents/new_aligned_shapesv4/'
data_dirs = list.dirs(base_dir,recursive = FALSE)

roc_curve1 = read.csv('~/Documents/new_aligned_shapesv4/V13_5peak_2gp_v1_bw/roc_dpc1.csv')
roc_curve2 = read.csv('~/Documents/new_aligned_shapesv4/V13_5peak_2gp_v1_bw/roc_dpc2.csv')

data_dirs2 = data_dirs[-1]

for (dir in data_dirs2){
  roc_curve1.5 = read.csv(paste(dir,'/roc_dpc1.csv',sep=''))
  roc_curve2.5 = read.csv(paste(dir,'/roc_dpc2.csv',sep=''))
  roc_curve1 = roc_curve1 + roc_curve1.5
  roc_curve2 = roc_curve2 + roc_curve2.5
}
roc_curve1 = roc_curve1/50
roc_curve2 = roc_curve2/50

roc_curve1_frame = data.frame(roc_curve1)
roc_curve2_frame = data.frame(roc_curve2)

ROC_curve_plt <- ggplot(data <- roc_curve1_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Directions") +
  ggtitle(sprintf("ROC Curve for Class 1 Teeth - 5 Peaks")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt)
ggsave('~/Dropbox/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/5_peaks/Class_1_roc_dir_per_cone.pdf')
ROC_curve_plt2 <- ggplot(data <- roc_curve2_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Directions") +
  ggtitle(sprintf("ROC Curve for Class 2 Teeth - 5 Peaks")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt2)
ggsave('~/Dropbox/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/5_peaks/Class_2_roc_dir_per_cone.pdf')



# Plotting now
color1='blue'
color2 = 'blue'
color3 = 'lightblue'
color4='palegreen'
color5='orangered'
col_pal=c(color1,color1,color2,color2,color3,color3,color4,color5)
colfunc <- colorRampPalette(col_pal)
rotation_matrix=matrix(c(0.99972576,0.02127766,0.00978078,0,0.01589701,-0.92330807,0.38373107,0,0.01719557,-0.38347024,-0.92339307,0,0,0,0,1),ncol=4,byrow=TRUE)

dir = data_dirs[4]
load(paste(dir,'/rocs_and_rate.R',sep=''))
old_data_dir = paste(dir,'/mesh/gp1',sep='')
new_data_dir = paste(dir,'/mesh/gp2',sep='')
old_data_files = list.files(old_data_dir, full.names = TRUE)
new_data_files = list.files(new_data_dir, full.names = TRUE)
class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE)
class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)
ind = 3
old_veg1 = vcgImport(old_data_files[ind])
new_veg1 = vcgImport(new_data_files[ind])
og_mesh = vcgImport(file = 'Data/new_aligned_shapesv3/Old_Follivore/clean_V13_sas.off')
new_veg_1 = process_off_file_v3(new_data_files[ind])
old_veg_1 = process_off_file_v3(old_data_files[ind])
v1 = t(new_veg1$vb[-4,])
v2 = t(old_veg1$vb[-4,])


new_heat =  reconstruct_vertices_on_shape(dir = pset$dirs,complex = new_veg_1,rate_vals = data_summary$Rate2[,2],
                                              len = pset$len,cuts = 300,cone_size = pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
old_heat =  reconstruct_vertices_on_shape(dir = pset$dirs,complex = old_veg_1,rate_vals = data_summary$Rate2[,2],
                                              len = pset$len,cuts = 300,cone_size = pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)

colfunc <- colorRampPalette(col_pal)
new_heat_colors=colfunc(max(new_heat[,1]) - min(new_heat[,1]))[new_heat[,1] - min(new_heat[,1])]
old_heat_colors=colfunc(max(old_heat[,1]) - min(old_heat[,1]))[old_heat[,1] - min(old_heat[,1])]
mfrow3d(2,3)
new_veg_colors = rep('white', dim(new_veg1$vb)[2])
plot3d(og_mesh,col = new_veg_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)
new_veg_colors[(which(class_2_probs > 0.25))] = 'red'
plot3d(new_veg1,col = new_veg_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(new_veg1,col = new_heat_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)

old_veg_colors = rep('white', dim(old_veg1$vb)[2])
plot3d(og_mesh,col = old_veg_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)
old_veg_colors[(which(class_1_probs > 0.25))] = 'red'

plot3d(old_veg1,col = old_veg_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)
plot3d(old_veg1,col = old_heat_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)

#### 1 v 1 ####
base_dir = '~/Documents/new_aligned_shapesv5/'
data_dirs = list.dirs(base_dir,recursive = FALSE)
data_dirs = data_dirs[-(1:9)]
for (dir in data_dirs){
  old_data_dir = paste(dir,'/mesh/gp1',sep='')
  new_data_dir = paste(dir,'/mesh/gp2',sep='')
  old_data_files = list.files(old_data_dir, full.names = TRUE)
  new_data_files = list.files(new_data_dir, full.names = TRUE)
  class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE)
  class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)
  roc_curves = list()
  total_dirs = c(5,15,25,35)
  #total_dirs = c(1)
  for (k in 1:length(total_dirs)){
    dirs = total_dirs[k]
    print(paste("on dir", dirs))
    pset = list(num_cones = dirs, cap_radius = 0.15, len = 100, directions_per_cone = 5,
                      dirs = generate_equidistributed_cones(num_directions = dirs, cap_radius =  0.15, directions_per_cone = 5))
    data_summary=real_data_summary(dir1=new_data_dir,dir2 = old_data_dir,direction=pset$dirs,class1='Vegetable', class2='Vegetable',
                                      radius=0,accuracy=FALSE,len = pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
    roc_curve = compute_roc_curve_teeth(data_dir1 = old_data_dir, data_dir2 = new_data_dir, gamma = 0.25,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                 rate_values = data_summary$Rate2[,2],directions_per_cone = pset$directions_per_cone,curve_length = pset$len,
                                 directions = pset$dirs,truncated = 500,ball_radius = ball_radius, ball = ball, radius = 1, two_curves = TRUE)
    #print(roc_curve)
    roc_curve[,3] = dirs
    roc_frame = data.frame(roc_curve)
    roc_curves[[k]] = roc_frame
  }
  roc_curves_total1 = roc_curves[[1]][1:500,]
  for (k in 2:length(total_dirs)){
    roc_curves_total1 = rbind(roc_curves_total1,roc_curves[[k]][1:500,])
  }
  write.csv(roc_curves_total1,file = paste(dir,'/roc_dirs1.csv',sep=''),row.names = FALSE)
  roc_curves_total2 = roc_curves[[1]][501:1000,]
  for (k in 2:length(total_dirs)){
    roc_curves_total2 = rbind(roc_curves_total2,roc_curves[[k]][501:1000,])
  }
  write.csv(roc_curves_total2,file = paste(dir,'/roc_dirs2.csv',sep=''),row.names = FALSE)
  save.image(paste(dir,'/rocs_and_rate.R',sep=''))
}
# Making the ROC Curves
base_dir = '~/Documents/new_aligned_shapesv4/'
data_dirs = list.dirs(base_dir,recursive = FALSE)

roc_curve1 = read.csv('~/Documents/new_aligned_shapesv4/V13_5peak_2gp_v1_bw/roc_dirs1.csv')
roc_curve2 = read.csv('~/Documents/new_aligned_shapesv4/V13_5peak_2gp_v1_bw/roc_dirs2.csv')

data_dirs2 = data_dirs[-1]

for (dir in data_dirs2){
  roc_curve1.5 = read.csv(paste(dir,'/roc_dirs1.csv',sep=''))
  roc_curve2.5 = read.csv(paste(dir,'/roc_dirs2.csv',sep=''))
  roc_curve1 = roc_curve1 + roc_curve1.5
  roc_curve2 = roc_curve2 + roc_curve2.5
}
roc_curve1 = roc_curve1/50
roc_curve2 = roc_curve2/50

roc_curve1_frame = data.frame(roc_curve1)
roc_curve2_frame = data.frame(roc_curve2)

ROC_curve_plt <- ggplot(data <- roc_curve1_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Cones") +
  ggtitle(sprintf("ROC Curve for Class 1 Teeth - 5 Peaks")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt)
ggsave('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/5_peaks/Class_1_roc.pdf')
ROC_curve_plt2 <- ggplot(data <- roc_curve2_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Cones") +
  ggtitle(sprintf("ROC Curve for Class 2 Teeth - 5 Peaks")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt2)
ggsave('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/5_peaks/Class_2_roc.pdf')

# Varying Length Up

base_dir = '~/Documents/new_aligned_shapesv4/'
data_dirs = list.dirs(base_dir,recursive = FALSE)
data_dirs = data_dirs[-(1:25)]
for (dir in data_dirs){
  old_data_dir = paste(dir,'/mesh/gp1',sep='')
  new_data_dir = paste(dir,'/mesh/gp2',sep='')
  old_data_files = list.files(old_data_dir, full.names = TRUE)
  new_data_files = list.files(new_data_dir, full.names = TRUE)
  class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE)
  class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)
  roc_curves = list()
  dirs = 15
  total_lens = c(10,25,50,75)
  #total_dirs = c(1)
  for (k in 1:length(total_lens)){
    len = total_lens[k]
    print(paste("on length", len, 'dir', dirs))
    pset = list(num_cones = dirs, cap_radius = 0.15, len = len, directions_per_cone = 5,
                      dirs = generate_equidistributed_cones(num_directions = dirs, cap_radius =  0.15, directions_per_cone = 5))
    data_summary=real_data_summary(dir1=new_data_dir,dir2 = old_data_dir,direction=pset$dirs,class1='Vegetable', class2='Vegetable',
                                      radius=0,accuracy=FALSE,len = pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
    roc_curve = compute_roc_curve_teeth(data_dir1 = old_data_dir, data_dir2 = new_data_dir, gamma = 0.25,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                 rate_values = data_summary$Rate2[,2],directions_per_cone = pset$directions_per_cone,curve_length = pset$len,
                                 directions = pset$dirs,truncated = 500,ball_radius = ball_radius, ball = ball, radius = 1, two_curves = TRUE)
    #print(roc_curve)
    roc_curve[,3] = len
    roc_frame = data.frame(roc_curve)
    roc_curves[[k]] = roc_frame
  }
  roc_curves_total1 = roc_curves[[1]][1:500,]
  for (k in 2:length(total_lens)){
    roc_curves_total1 = rbind(roc_curves_total1,roc_curves[[k]][1:500,])
  }
  write.csv(roc_curves_total1,file = paste(dir,'/roc_lens1.csv',sep=''),row.names = FALSE)
  roc_curves_total2 = roc_curves[[1]][501:1000,]
  for (k in 2:length(total_lens)){
    roc_curves_total2 = rbind(roc_curves_total2,roc_curves[[k]][501:1000,])
  }
  write.csv(roc_curves_total2,file = paste(dir,'/roc_lens2.csv',sep=''),row.names = FALSE)
  save.image(paste(dir,'/rocs_and_rate_lens.R',sep=''))
}


# Making the ROC Curves
base_dir = '~/Documents/new_aligned_shapesv4/'
data_dirs = list.dirs(base_dir,recursive = FALSE)

roc_curve_len_1 = read.csv('~/Documents/new_aligned_shapesv4/V13_5peak_2gp_v1_bw/roc_lens1.csv')
roc_curve_len_2 = read.csv('~/Documents/new_aligned_shapesv4/V13_5peak_2gp_v1_bw/roc_lens2.csv')

data_dirs = data_dirs[-1]

for (dir in data_dirs){
  roc_curve_len_1.5 = read.csv(paste(dir,'/roc_lens1.csv',sep=''))
  roc_curve_len_2.5 = read.csv(paste(dir,'/roc_lens2.csv',sep=''))
  roc_curve_len_1 = roc_curve_len_1 + roc_curve_len_1.5
  roc_curve_len_2 = roc_curve_len_2 + roc_curve_len_2.5
}
roc_curve_len_1 = roc_curve_len_1/50
roc_curve_len_2 = roc_curve_len_2/50


roc_curve_1_100 = roc_curve1[roc_curve1[,3] == 15,]
roc_curve_2_100 = roc_curve2[roc_curve2[,3] == 15,]
roc_curve_1_100[,3] = 100
roc_curve_2_100[,3] = 100


roc_curve_len_1 = rbind(roc_curve_len_1,roc_curve_1_100)
roc_curve_len_2 = rbind(roc_curve_len_2,roc_curve_2_100)

roc_curve_len_1_frame = data.frame(roc_curve_len_1)
roc_curve_len_2_frame = data.frame(roc_curve_len_2)

ROC_curve_plt <- ggplot(data <- roc_curve_len_1_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Lengths") +
  ggtitle(sprintf("ROC Curve for Class 1 Teeth - 5 Peaks - Varying Length")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt)
ggsave('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/5_peaks/Class_1_roc_lengths.pdf')
ROC_curve_plt2 <- ggplot(data <- roc_curve_len_2_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Lengths") +
  ggtitle(sprintf("ROC Curve for Class 2 Teeth - 5 Peaks - Varying Length")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt2)
ggsave('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Old Drafts/Draft 2/Figures/CaricaturedTeeth/5_peaks/Class_2_roc_lengths.pdf')

# Plotting now
color1='blue'
color2 = 'blue'
color3 = 'lightblue'
color4='palegreen'
color5='orangered'
col_pal=c(color1,color1,color2,color2,color3,color3,color4,color5)
colfunc <- colorRampPalette(col_pal)
rotation_matrix=matrix(c(0.99972576,0.02127766,0.00978078,0,0.01589701,-0.92330807,0.38373107,0,0.01719557,-0.38347024,-0.92339307,0,0,0,0,1),ncol=4,byrow=TRUE)

dir = data_dirs[2]
load(paste(dir,'/rocs_and_rate.R',sep=''))
old_data_dir = paste(dir,'/mesh/gp1',sep='')
new_data_dir = paste(dir,'/mesh/gp2',sep='')
old_data_files = list.files(old_data_dir, full.names = TRUE)
new_data_files = list.files(new_data_dir, full.names = TRUE)
class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE)
class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)
ind = 3
old_veg1 = vcgImport(old_data_files[ind])
new_veg1 = vcgImport(new_data_files[ind])
og_mesh = vcgImport(file = 'Data/new_aligned_shapesv3/Old_Follivore/clean_V13_sas.off')
new_veg_1 = process_off_file_v3(new_data_files[ind])
old_veg_1 = process_off_file_v3(old_data_files[ind])
v1 = t(new_veg1$vb[-4,])
v2 = t(old_veg1$vb[-4,])


new_heat =  reconstruct_vertices_on_shape(dir = pset$dirs,complex = new_veg_1,rate_vals = data_summary$Rate2[,2],
                                              len = pset$len,cuts = 300,cone_size = pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
old_heat =  reconstruct_vertices_on_shape(dir = pset$dirs,complex = old_veg_1,rate_vals = data_summary$Rate2[,2],
                                              len = pset$len,cuts = 300,cone_size = pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)

colfunc <- colorRampPalette(col_pal)
new_heat_colors=colfunc(max(new_heat[,1]) - min(new_heat[,1]))[new_heat[,1] - min(new_heat[,1])]
old_heat_colors=colfunc(max(old_heat[,1]) - min(old_heat[,1]))[old_heat[,1] - min(old_heat[,1])]
mfrow3d(2,3)
new_veg_colors = rep('white', dim(new_veg1$vb)[2])
plot3d(og_mesh,col = new_veg_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)
new_veg_colors[(which(class_2_probs > 0.25))] = 'red'
plot3d(new_veg1,col = new_veg_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(new_veg1,col = new_heat_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)

old_veg_colors = rep('white', dim(old_veg1$vb)[2])
plot3d(og_mesh,col = old_veg_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)
old_veg_colors[(which(class_1_probs > 0.25))] = 'red'

plot3d(old_veg1,col = old_veg_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)
plot3d(old_veg1,col = old_heat_colors, axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
rgl.viewpoint(userMatrix = rotation_matrix)
