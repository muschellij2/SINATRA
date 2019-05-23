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
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Number of Cones") +
  ggtitle(sprintf("ROC Curve for Class 1 Teeth")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt)
ROC_curve_plt2 <- ggplot(data <- roc_curve2_frame,aes(x = X1, y = X2, group = X3)) +
  geom_line(stat = "identity",aes(color = factor(X3) )) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Number of Cones") +
  ggtitle(sprintf("ROC Curve for Class 2 Teeth")) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt2)




#### 5 v 5 ####
base_dir = '~/Documents/new_aligned_shapesv4/'
data_dirs = list.dirs(base_dir,recursive = FALSE)
data_dirs = data_dirs[-(1:11)]
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
#### Supplemental stuff####

roc_curves_len = list()
total_lens = c(10,25,50,75,100)
class_1_probs = read.csv('Data/new_aligned_shapesv3/V13_v2/gp1_spt.csv', header = FALSE)
class_2_probs = read.csv('Data/new_aligned_shapesv3/V13_v2/gp2_spt.csv', header = FALSE)
for (k in 1:length(total_lens)){
  lens = total_lens[k]
  print(paste("on len", lens))
  pset = list(num_cones = 35, cap_radius = 0.15, len = lens, directions_per_cone = 5,
                  dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))
  data_summary=real_data_summary(dir1=new_veg_dir,dir2 = old_veg_dir,direction=pset$dirs,class1='Vegetable', class2='Vegetable',
                                radius=0,accuracy=FALSE,len = pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
  c9 = compute_roc_curve_teeth(data_dir1 = old_veg_dir, data_dir2 = new_veg_dir, gamma = 0.5,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                             rate_values = data_summary$Rate2[,2],directions_per_cone = pset$directions_per_cone,curve_length = pset$len,
                             directions = pset$dirs,truncated = 300,ball_radius = ball_radius, ball = ball, radius = 1)
  c9[,3] = c9[,3] +k
  c9_frame = data.frame(c9)
  roc_curves_len[[k]] = c9_frame
}
roc_curves_len_total = roc_curves_len[[1]][301:600,]
roc_curves_len_total[,3] = 10
for (k in 2:length(total_lens)){
  roc_curves_len_total = rbind(roc_curves_len_total,roc_curves_len[[k]][301:600,])
  roc_curves_len_total[((k-1)*300+1): ((k)*300) ,3] = total_lens[k]
}
ggplot(roc_curves_len_total, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 1,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Sub Level Sets') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "ROC Curve Varying Sub Level Sets")

save.image('Caricature_ROCV3.Rdata')

color1='blue'
color2 = 'blue'
color3 = 'lightblue'
color4='palegreen'
color5='orangered'
col_pal=c(color1,color1,color2,color2,color3,color3,color4,color5)
colfunc <- colorRampPalette(col_pal)
rotation_matrix=matrix(c(0.99972576,0.02127766,0.00978078,0,0.01589701,-0.92330807,0.38373107,0,0.01719557,-0.38347024,-0.92339307,0,0,0,0,1),ncol=4,byrow=TRUE)

#### Old Veg New Veg ####
dir = data_dirs[50]
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


open3d()
mfrow3d(1,2)
#plot_selected_landmarks(old_fruit1,old_fruit1_vert,num_landmarks = 20)
new_veg1_vert = compute_selected_vertices_cones(dir = pset$dirs, complex =new_veg_1, rate_vals = data_summary$Rate2[,2], len = pset$len, threshold = (pset$directions_per_cone/(3*dim(data_summary$Rate2)[1])),
                                                cone_size = pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

old_veg1_vert = compute_selected_vertices_cones(dir = pset$dirs, complex = old_veg_1, rate_vals = data_summary$Rate2[,2], len = pset$len, threshold = (pset$directions_per_cone/(3*dim(data_summary$Rate2)[1])),
                                                cone_size = pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
new_veg_colors = rep('white', dim(new_veg1$vb)[2])
new_veg_colors[new_veg1_vert] = 'red'
plot3d(new_veg1,col = new_veg_colors,specular="black", axes = FALSE)
#plot_selected_landmarks(new_veg1,new_veg1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

old_veg_colors = rep('white', dim(old_veg1$vb)[2])
old_veg_colors[old_veg1_vert] = 'blue'
plot3d(old_veg1,col = old_veg_colors,specular="black", axes = FALSE)
#plot_selected_landmarks(old_veg1,old_veg1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

veg_new_heat =  reconstruct_vertices_on_shape(dir = pset$dirs,complex = new_veg_1,rate_vals = data_summary$Rate2[,2],
                                              len = pset$len,cuts = 300,cone_size = pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
veg_old_heat =  reconstruct_vertices_on_shape(dir = pset$dirs,complex = old_veg_1,rate_vals = data_summary$Rate2[,2],
                                              len = pset$len,cuts = 300,cone_size = pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

#veg_new_colors=colfunc(max(veg_new_heat[,1]))[veg_new_heat[,1]]
#veg_old_colors=colfunc(max(veg_old_heat[,1]))[veg_old_heat[,1]]

col_pal=c(color1,'blue','lightblue','lightgreen','lightgreen','red')
colfunc <- colorRampPalette(col_pal)
veg_new_colors=colfunc(max(veg_new_heat[,1]) - min(veg_new_heat[,1]))[veg_new_heat[,1] - min(veg_new_heat[,1])]
veg_old_colors=colfunc(max(veg_old_heat[,1]) - min(veg_old_heat[,1]))[veg_old_heat[,1] - min(veg_old_heat[,1])]
mfrow3d(1,2)
plot3d(new_veg1,col = veg_new_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(new_veg1,new_veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(old_veg1,col = veg_old_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(old_veg1,old_veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)


mfrow3d(2,3)
og_mesh = vcgImport(file = 'Data/new_aligned_shapesv3/Old_Follivore/clean_V13_sas.off')
plot3d(og_mesh, col = 'white', axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
new_veg_colors = rep('white', dim(new_veg1$vb)[2])
new_veg_colors[(which(class_2_probs > 0.05))] = 'green'
plot3d(new_veg1,col = new_veg_colors,specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
plot3d(new_veg1,col = veg_new_colors, specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(og_mesh, col = 'white',specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
old_veg_colors = rep('white', dim(old_veg1$vb)[2])
old_veg_colors[(which(class_1_probs > 0.05))] = 'green'

plot3d(old_veg1,col = old_veg_colors, axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
#plot_selected_landmarks(new_veg1,new_veg1_vert,num_landmarks = 20)

plot3d(old_veg1,col = veg_old_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(old_veg1,old_veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)


