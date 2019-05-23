setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
library(R.utils)
sourceDirectory('~/Documents/SINATRA/SINATRA_Pipeline_Branch/')
sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')
load('Caricature_ROCV2.Rdata')
#library(svd)
#set.seed(15)
#set.seed(35)
set.seed(15)

#Specifying Directories and the Associated Files
old_five_dir='Data/new_aligned_shapesv3/V13_10peak_2gp_mix/mesh/gp1'
new_five_dir='Data/new_aligned_shapesv3/V13_10peak_2gp_mix/mesh/gp2'
old_five_files = list.files(old_five_dir, full.names = TRUE)
new_five_files = list.files(new_five_dir, full.names = TRUE)
class_1_probs5 = read.csv('Data/new_aligned_shapesv3/V13_10peak_2gp_mix/gp1_spt.csv', header = FALSE)
class_2_probs5 = read.csv('Data/new_aligned_shapesv3/V13_10peak_2gp_mix/gp2_spt.csv', header = FALSE)


#Parameters for the Analysis
ball = TRUE
ball_radius = 0.5
ec_type = 'ECT'
#### Start Comparison #### 


#### 5 Causal Landmarks ####
roc_curves5 = list()
#total_dirs = c(2,5,10,15,20,25,30,35)
total_dirs5 = c(35)
for (k in 1:length(total_dirs5)){
  dirs = total_dirs5[k]
  print(paste("on dir", dirs))
  #five_pset = list(num_cones = dirs, cap_radius = 0.15, len = 100, directions_per_cone = 5,
  #                 dirs = generate_equidistributed_cones(num_directions = dirs, cap_radius =  0.15, directions_per_cone = 5))
  #five_five_old=real_data_summary(dir1=new_five_dir,dir2 = old_five_dir,direction=five_pset$dirs,class1='Vegetable', class2='Vegetable',
  #                                radius=0,accuracy=FALSE,len = five_pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
  c5 = compute_roc_curve_teeth(data_dir1 = old_five_dir, data_dir2 = new_five_dir, gamma = 0.25,class_1_probs = class_1_probs5,class_2_probs = class_2_probs5,
                               rate_values = five_five_old$Rate2[,2],directions_per_cone = five_pset$directions_per_cone,curve_length = five_pset$len,
                               directions = five_pset$dirs,truncated = 300,ball_radius = ball_radius, ball = ball, radius = 1)
  c5[,3] = c5[,3] +k
  c5_frame = data.frame(c5)
  roc_curves5[[k]] = c5_frame
}
#save.image('Caricature_ROC.Rdata')
ggplot(roc_curve_total_10_peakss, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 1,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Detection Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Number of Cones') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "ROC Curve for Caricatured Teeth")

#### 3 Causal Landmarks ####
old_three_dir='Data/new_aligned_shapesv3/V13_v2/gp1'
new_three_dir='Data/new_aligned_shapesv3/V13_v2/gp2'
old_three_files = list.files(old_three_dir, full.names = TRUE)
new_three_files = list.files(new_three_dir, full.names = TRUE)
class_1_probs3 = read.csv('Data/new_aligned_shapesv3/V13_v2/gp1_spt.csv', header = FALSE)
class_2_probs3 = read.csv('Data/new_aligned_shapesv3/V13_v2/gp2_spt.csv', header = FALSE)


#Parameters for the Analysis
ball = TRUE
ball_radius = 0.5
ec_type = 'ECT'

# loop over directions
roc_curves3 = list()
#total_dirs = c(2,5,10,15,20,25,30,35)
total_dirs3 = c(35)
for (k in 1:length(total_dirs3)){
  dirs = total_dirs3[k]
  print(paste("on dir", dirs))
  #three_pset = list(num_cones = dirs, cap_radius = 0.15, len = 100, directions_per_cone = 5,
  #                  dirs = generate_equidistributed_cones(num_directions = dirs, cap_radius =  0.15, directions_per_cone = 5))
  #three_three_old=real_data_summary(dir1=new_three_dir,dir2 = old_three_dir,direction=three_pset$dirs,class1='Vegetable', class2='Vegetable',
  #                                  radius=0,accuracy=FALSE,len = three_pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
  c3 = compute_roc_curve_teeth(data_dir1 = old_three_dir, data_dir2 = new_three_dir, gamma = 0.25,class_1_probs = class_1_probs3,class_2_probs = class_2_probs3,
                               rate_values = three_three_old$Rate2[,2],directions_per_cone = three_pset$directions_per_cone,curve_length = three_pset$len,
                               directions = three_pset$dirs,truncated = 300,ball_radius = ball_radius, ball = ball, radius = 1)
  c3[,3] = c3[,3] +k
  c3_frame = data.frame(c3)
  roc_curves3[[k]] = c3_frame
}

#### 1 Causal Landmarks ####

old_one_dir='Data/new_aligned_shapesv3/V13_2peak_2gp/mesh/gp1'
new_one_dir='Data/new_aligned_shapesv3/V13_2peak_2gp/mesh/gp2'
old_one_files = list.files(old_one_dir, full.names = TRUE)
new_one_files = list.files(new_one_dir, full.names = TRUE)
class_1_probs1 = read.csv('Data/new_aligned_shapesv3/V13_2peak_2gp/gp1_spt.csv', header = FALSE)
class_2_probs1 = read.csv('Data/new_aligned_shapesv3/V13_2peak_2gp/gp2_spt.csv', header = FALSE)


#Parameters for the Analysis
ball = TRUE
ball_radius = 0.5
ec_type = 'ECT'

# loop over directions
roc_curves1 = list()
#total_dirs = c(2,5,10,15,20,25,30,35)
total_dirs1 = c(35)
for (k in 1:length(total_dirs1)){
  dirs = total_dirs1[k]
  print(paste("on dir", dirs))
  #one_pset = list(num_cones = dirs, cap_radius = 0.15, len = 100, directions_per_cone = 5,
  #                dirs = generate_equidistributed_cones(num_directions = dirs, cap_radius =  0.15, directions_per_cone = 5))
  #one_one_old=real_data_summary(dir1=new_one_dir,dir2 = old_one_dir,direction=one_pset$dirs,class1='Vegetable', class2='Vegetable',
  #                              radius=0,accuracy=FALSE,len = one_pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
  c1 = compute_roc_curve_teeth(data_dir1 = old_one_dir, data_dir2 = new_one_dir, gamma = 0.25,class_1_probs = class_1_probs1,class_2_probs = class_2_probs1,
                               rate_values = one_one_old$Rate2[,2],directions_per_cone = one_pset$directions_per_cone,curve_length = one_pset$len,
                               directions = one_pset$dirs,truncated = 300,ball_radius = ball_radius, ball = ball, radius = 1)
  c1[,3] = c1[,3] +k
  c1_frame = data.frame(c1)
  roc_curves1[[k]] = c1_frame
}
#### Combine the ROC curves ####
roc_curve_total_3_peaks = roc_curves3[[1]][301:600,]
roc_curve_total_3_peaks = (roc_curves3[[1]][301:600,] + roc_curves3[[1]][1:300,])/2
roc_curve_total_3_peaks[,3] = 3
roc_curve_total_10_peakss = roc_curves5[[1]][1:300,]
roc_curve_total_10_peakss = (roc_curves5[[1]][301:600,] + roc_curves5[[1]][1:300,])/2
roc_curve_total_10_peakss[,3] = 5
roc_curve_total_1_peakss = roc_curves1[[1]][1:300,]
roc_curve_total_1_peakss = (roc_curves1[[1]][301:600,] + roc_curves1[[1]][1:300,])/2
roc_curve_total_1_peakss[,3] = 1

roc_curve_total = rbind(roc_curve_total_3_peaks,roc_curve_total_10_peakss,roc_curve_total_1_peakss)
ggplot(roc_curve_total, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 1,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Detection Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Number of Causal Landmarks') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "ROC Curve for Caricatured Teeth")

### Supplemental stuff####

roc_curves_len = list()
total_lens = c(10,25,50,75,100)
class_1_probs = read.csv('Data/new_aligned_shapesv3/V13_v2/gp1_spt.csv', header = FALSE)
class_2_probs = read.csv('Data/new_aligned_shapesv3/V13_v2/gp2_spt.csv', header = FALSE)
for (k in 1:length(total_lens)){
  lens = total_lens[k]
  print(paste("on len", lens))
  veg_pset = list(num_cones = 35, cap_radius = 0.15, len = lens, directions_per_cone = 5,
                  dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))
  veg_veg_old=real_data_summary(dir1=new_veg_dir,dir2 = old_veg_dir,direction=veg_pset$dirs,class1='Vegetable', class2='Vegetable',
                                radius=0,accuracy=FALSE,len = veg_pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
  c9 = compute_roc_curve_teeth(data_dir1 = old_veg_dir, data_dir2 = new_veg_dir, gamma = 0.5,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                             rate_values = veg_veg_old$Rate2[,2],directions_per_cone = veg_pset$directions_per_cone,curve_length = veg_pset$len,
                             directions = veg_pset$dirs,truncated = 300,ball_radius = ball_radius, ball = ball, radius = 1)
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
old_veg_dir='Data/new_aligned_shapesv3/V13_10peak_2gp/mesh/gp1'
new_veg_dir='Data/new_aligned_shapesv3/V13_10peak_2gp/mesh/gp2'
old_veg_files = list.files(old_veg_dir, full.names = TRUE)
new_veg_files = list.files(new_veg_dir, full.names = TRUE)
class_1_probs = read.csv('Data/new_aligned_shapesv3/V13_10peak_2gp/gp1_spt.csv', header = FALSE)
class_2_probs = read.csv('Data/new_aligned_shapesv3/V13_10peak_2gp/gp2_spt.csv', header = FALSE)
ind = 3
old_veg1 = vcgImport(old_veg_files[ind])
new_veg1 = vcgImport(new_veg_files[ind])
new_veg_1 = process_off_file_v3(new_veg_files[ind])
old_veg_1 = process_off_file_v3(old_veg_files[ind])
v1 = t(new_veg1$vb[-4,])
v2 = t(old_veg1$vb[-4,])

mfrow3d(2,2)
new_veg_colors = rep('white', dim(new_veg1$vb)[2])
plot3d(new_veg1,col = new_veg_colors, axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
new_veg_colors[(which(class_2_probs > 0.25))] = 'red'
plot3d(new_veg1,col = new_veg_colors, axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)

old_veg_colors = rep('white', dim(old_veg1$vb)[2])
plot3d(old_veg1,col = old_veg_colors, axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
old_veg_colors[(which(class_1_probs > 0.25))] = 'red'

plot3d(old_veg1,col = old_veg_colors, axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
#plot3d(new_veg1,col = veg_new_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(new_veg1,new_veg1_vert,num_landmarks = 20)
#rgl.viewpoint(userMatrix = rotation_matrix)
#
#plot3d(old_veg1,col = veg_old_colors, specular="black", axes = FALSE)
##plot_selected_landmarks(old_veg1,old_veg1_vert,num_landmarks = 20)
#rgl.viewpoint(userMatrix = rotation_matrix)

open3d()
mfrow3d(1,2)
#plot_selected_landmarks(old_fruit1,old_fruit1_vert,num_landmarks = 20)
new_veg1_vert = compute_selected_vertices_cones(dir = veg_pset$dirs, complex =new_veg_1, rate_vals = veg_veg_old$Rate2[,2], len = veg_pset$len, threshold = (veg_pset$directions_per_cone/(3*dim(veg_veg_old$Rate2)[1])),
                                                cone_size = veg_pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

old_veg1_vert = compute_selected_vertices_cones(dir = veg_pset$dirs, complex = old_veg_1, rate_vals = veg_veg_old$Rate2[,2], len = veg_pset$len, threshold = (veg_pset$directions_per_cone/(3*dim(veg_veg_old$Rate2)[1])),
                                                cone_size = veg_pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
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

veg_new_heat =  reconstruct_vertices_on_shape(dir = veg_pset$dirs,complex = new_veg_1,rate_vals = veg_veg_old$Rate2[,2],
                                              len = veg_pset$len,cuts = 300,cone_size = veg_pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
veg_old_heat =  reconstruct_vertices_on_shape(dir = veg_pset$dirs,complex = old_veg_1,rate_vals = veg_veg_old$Rate2[,2],
                                              len = veg_pset$len,cuts = 300,cone_size = veg_pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

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


