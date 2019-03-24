setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
library(R.utils)
sourceDirectory('~/Documents/SINATRA/SINATRA_Pipeline_Branch/')
sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')

#library(svd)
set.seed(15)

#Specifying Directories and the Associated Files

#old_fruit_dir = 'Data/HDM/by_diet_scaled_aligned/Old_Frugivore/'
#old_veg_dir='Data/HDM/by_diet_scaled_aligned/Old_Follivore/'
old_bug_dir='Data/HDM//by_diet_scaled_aligned/Old_Insectivore/'
#new_fruit_dir = 'Data/HDM/by_diet_scaled_aligned/New_Frugivore/'
#new_veg_dir='Data/HDM/by_diet_scaled_aligned/New_Follivore/'
#new_bug_dir='Data/HDM//by_diet_scaled_aligned/New_Insectivore/'


old_bug_dir='Data/new_aligned_shapesv3/V13_v4/gp1'
new_bug_dir='Data/new_aligned_shapesv3/V13_v4/gp2'
old_fruit_dir = 'Data/new_aligned_shapesv3/V13_v5/gp1'
new_fruit_dir = 'Data/new_aligned_shapesv3/V13_v5/gp2'
old_veg_dir='Data/new_aligned_shapesv3/V13_v2/gp1'
new_veg_dir='Data/new_aligned_shapesv3/V13_v2/gp2'
old_fruit_dir = 'Data/new_aligned_shapesv3/h16/gp1'
new_fruit_dir = 'Data/new_aligned_shapesv3/h16/gp3'

old_fruit_files=list.files(path=old_fruit_dir,full.names = TRUE)
old_veg_files=list.files(path=old_veg_dir,full.names = TRUE)
old_bug_files=list.files(path=old_bug_dir,full.names = TRUE)

new_fruit_files=list.files(path=new_fruit_dir,full.names = TRUE)
new_veg_files=list.files(path=new_veg_dir,full.names = TRUE)
new_bug_files=list.files(path=new_bug_dir,full.names = TRUE)


files = c(old_fruit_files, old_veg_files, old_bug_files, new_fruit_files, new_veg_files, new_bug_files)
#Parameters for the Analysis
cap_radius = 0.15
num_cones = 25
directions_per_cone = 5
len = 100
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  cap_radius, directions_per_cone = directions_per_cone)
rot_mat = generate_random_rotation(n = 3)
new_dirs = dirs %*% rot_mat
#plot3d(dirs, size = 3, col = 'blue')
rgl.points(new_dirs, size = 3, col = 'red')
plot3d(new_dirs, size = 3, col = 'red')
ball = FALSE
ball_radius = 0.5
ec_type = 'ECT'
#### Start Comparison #### 

fruit_fruit_old=real_data_summary(dir1=new_fruit_dir,dir2 = old_fruit_dir,direction=dirs,class1='Fruit', class2='Fruit',
                                              radius=0,accuracy=FALSE,len = len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)

veg_veg_old=real_data_summary(dir1=new_veg_dir,dir2 = old_veg_dir,direction=dirs,class1='Vegetable', class2='Vegetable',
                                          radius=0,accuracy=FALSE,len = len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
bug_bug_old=real_data_summary(dir1=new_bug_dir,dir2 = old_bug_dir,direction=dirs,class1='Bug', class2='Bug',
                                          radius=0,accuracy=FALSE,len = len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)


color1='blue'
color2 = 'blue'
color3 = 'lightblue'
color4='palegreen'
color5='orangered'
col_pal=c(color1,color1,color2,color2,color3,color3,color4,color4,color5)
#col_pal =c(color1,color2,color3)
colfunc <- colorRampPalette(col_pal)
rotation_matrix=matrix(c(0.99972576,0.02127766,0.00978078,0,0.01589701,-0.92330807,0.38373107,0,0.01719557,-0.38347024,-0.92339307,0,0,0,0,1),ncol=4,byrow=TRUE)

#### Old Fruit Fruit ####
class_1_probs = read.csv('Data/new_aligned_shapesv3/h16/gp1_spt.csv',header = TRUE)
class_2_probs = read.csv('Data/new_aligned_shapesv3/h16/gp3_spt.csv', header = TRUE)
ind = 10
new_fruit1 = vcgImport(new_fruit_files[ind])
new_fruit_1 = process_off_file_v3(new_fruit_files[ind])

ind = 10
old_fruit1 = vcgImport(old_fruit_files[ind])
old_fruit_1 = process_off_file_v3(old_fruit_files[ind])

mfrow3d(1,2)
new_fruit1_vert = compute_selected_vertices_cones(dir = dirs, complex =new_fruit_1, rate_vals = fruit_fruit_old$Rate2[,2], len = len, threshold = (directions_per_cone/(1*dim(fruit_fruit_old$Rate2)[1])),
                                                  cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

old_fruit1_vert = compute_selected_vertices_cones(dir = dirs, complex = old_fruit_1, rate_vals = fruit_fruit_old$Rate2[,2], len = len, threshold = (directions_per_cone/(1*dim(fruit_fruit_old$Rate2)[1])),
                                                  cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
new_fruit_colors = rep('white', dim(new_fruit1$vb)[2])
new_fruit_colors[new_fruit1_vert] = 'red'
v1 = t(new_fruit1$vb[-4,])
v2 = t(old_fruit1$vb[-4,])
plot3d(new_fruit1,col = new_fruit_colors)
rgl.points(v1[(1 + which(class_2_probs[,1] > 0.05)),], size = 3, col = 'green')
#plot_selected_landmarks(new_fruit1,new_fruit1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

old_fruit_colors = rep('white', dim(old_fruit1$vb)[2])
old_fruit_colors[old_fruit1_vert] = 'blue'
plot3d(old_fruit1,col = old_fruit_colors)
rgl.points(v2[which(class_1_probs[,1] > 0.0),], size = 3, col = 'green')
#plot_selected_landmarks(old_fruit1,old_fruit1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)


fruit_new_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = new_fruit_1,rate_vals = fruit_fruit_old$Rate2[,2],
                                                len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
fruit_old_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = old_fruit_1,rate_vals = fruit_fruit_old$Rate2[,2],
                                                len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

#col_pal=c(color1,'lightgreen','orangered',color5)
#col_pal =c(color1,color2,color3)
colfunc <- colorRampPalette(col_pal)
fruit_new_colors=colfunc(max(fruit_new_heat[,1]))[(fruit_new_heat[,1] - min(fruit_new_heat[,1]))]
fruit_old_colors=colfunc(max(fruit_old_heat[,1]))[(fruit_old_heat[,1] -min(fruit_old_heat[,1]))]

fruit_new_colors=colfunc(max(fruit_new_heat[,1]))[fruit_new_heat[,1]]
fruit_old_colors=colfunc(max(fruit_old_heat[,1]))[fruit_old_heat[,1]]

mfrow3d(1,2)
plot3d(new_fruit1,col = fruit_new_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(new_fruit1,new_fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(old_fruit1,col = fruit_old_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(old_fruit1,old_fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

#class_1_probs = read.csv('Data/new_aligned_shapesv3/V13_v3/gp1_spt.csv')
#class_2_probs = read.csv('Data/new_aligned_shapesv3/V13_v3/gp2_spt.csv')

c1 = compute_roc_curve_teeth(data_dir1 = old_fruit_dir, data_dir2 = new_fruit_dir, gamma = 0.7,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                   rate_values = fruit_fruit_old$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = len,
                                   directions = dirs,truncated = 50,class = 1,ball_radius = ball_radius, ball = ball)
c1 = compute_roc_curve_teeth(data_dir1 = new_fruit_dir, data_dir2 = old_fruit_dir, gamma = 0.8,class_1_probs = class_2_probs,class_2_probs = class_1_probs,
                                   rate_values = fruit_fruit_old$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = len,
                                   directions = dirs,truncated = 50,class = 1,ball_radius = ball_radius, ball = ball)

c1_frame = data.frame(c1)
ggplot(c1_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "V3 RATE ROC with Gamma = 0.8")

c2 = compute_roc_curve_teeth(data_dir1 = old_fruit_dir, data_dir2 = new_fruit_dir, gamma = 0.00001,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                   rate_values = fruit_fruit_old$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = len,
                                   directions = dirs,truncated = 50,class = 1,ball_radius = ball_radius, ball = ball)
c2 = compute_roc_curve_teeth(data_dir1 = old_fruit_dir, data_dir2 = new_fruit_dir, gamma = 0.7,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                             rate_values = fruit_fruit_old$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = len,
                             directions = dirs,truncated = 50,class = 1,ball_radius = ball_radius, ball = ball)

c2_frame = data.frame(c2)
ggplot(c2_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "V3 RATE ROC with Gamma = 0.7")
c3 = compute_roc_curve_teeth(data_dir1 = old_fruit_dir, data_dir2 = new_fruit_dir, gamma = 0.000001,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                   rate_values = fruit_fruit_old$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = len,
                                   directions = dirs,truncated = 50,class = 1,ball_radius = ball_radius, ball = ball)

c3_frame = data.frame(c3)
ggplot(c3_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "V3 RATE ROC with Gamma = 0.000001")

#### Old Veg New Veg ####
ind = 3
old_veg1 = vcgImport(old_veg_files[ind])
new_veg1 = vcgImport(new_veg_files[ind])
new_veg_1 = process_off_file_v3(new_veg_files[ind])
old_veg_1 = process_off_file_v3(old_veg_files[ind])

mfrow3d(1,2)
new_veg1_vert = compute_selected_vertices_cones(dir = dirs, complex =new_veg_1, rate_vals = veg_veg_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(veg_veg_old$Rate2)[1])),
                                                cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

old_veg1_vert = compute_selected_vertices_cones(dir = dirs, complex = old_veg_1, rate_vals = veg_veg_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(veg_veg_old$Rate2)[1])),
                                                cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
new_veg_colors = rep('white', dim(new_veg1$vb)[2])
new_veg_colors[new_veg1_vert] = 'red'
plot3d(new_veg1,col = new_veg_colors)
#plot_selected_landmarks(new_veg1,new_veg1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

old_veg_colors = rep('white', dim(old_veg1$vb)[2])
old_veg_colors[old_veg1_vert] = 'blue'
plot3d(old_veg1,col = old_veg_colors)
#plot_selected_landmarks(old_veg1,old_veg1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

veg_new_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = new_veg_1,rate_vals = veg_veg_old$Rate2[,2],
                                              len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
veg_old_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = old_veg_1,rate_vals = veg_veg_old$Rate2[,2],
                                              len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

veg_new_colors=colfunc(max(veg_new_heat[,1]))[veg_new_heat[,1]]
veg_old_colors=colfunc(max(veg_old_heat[,1]))[veg_old_heat[,1]]

mfrow3d(1,2)
plot3d(new_veg1,col = veg_new_colors, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(new_veg1,new_veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(old_veg1,col = veg_old_colors, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(old_veg1,old_veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)


class_1_probs = read.csv('Data/new_aligned_shapesv3/V13_v2/gp1_spt.csv')
class_2_probs = read.csv('Data/new_aligned_shapesv3/V13_v2/gp2_spt.csv')

c4 = compute_roc_curve_teeth(data_dir1 = old_veg_dir, data_dir2 = new_veg_dir, gamma = 0.5,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                             rate_values = veg_veg_old$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = len,
                             directions = dirs,truncated = 50,class = 1,ball_radius = ball_radius, ball = ball)

c4_frame = data.frame(c4)
ggplot(c4_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "V2 RATE ROC with Gamma = 0.5")

c5 = compute_roc_curve_teeth(data_dir1 = old_veg_dir, data_dir2 = new_veg_dir, gamma = 0.6,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                             rate_values = veg_veg_old$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = len,
                             directions = dirs,truncated = 50,class = 1,ball_radius = ball_radius, ball = ball)

c5_frame = data.frame(c5)
ggplot(c5_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "V2 RATE ROC with Gamma = 0.6")
c6 = compute_roc_curve_teeth(data_dir1 = old_veg_dir, data_dir2 = new_veg_dir, gamma = 0.7,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                             rate_values = veg_veg_old$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = len,
                             directions = dirs,truncated = 50,class = 1,ball_radius = ball_radius, ball = ball)

c6_frame = data.frame(c6)
ggplot(c6_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "V2 RATE ROC with Gamma = 0.7")




#### Old Bug New Bug ####
ind = 3
old_bug1 = vcgImport(old_bug_files[ind])
new_bug1 = vcgImport(new_bug_files[ind])
new_bug_1 = process_off_file_v3(new_bug_files[ind])
old_bug_1 = process_off_file_v3(old_bug_files[ind])

mfrow3d(1,2)
new_bug1_vert = compute_selected_vertices_cones(dir = dirs, complex =new_bug_1, rate_vals = bug_bug_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(bug_bug_old$Rate2)[1])),
                                                cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

old_bug1_vert = compute_selected_vertices_cones(dir = dirs, complex = old_bug_1, rate_vals = bug_bug_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(bug_bug_old$Rate2)[1])),
                                                cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
new_bug_colors = rep('white', dim(new_bug1$vb)[2])
new_bug_colors[new_bug1_vert] = 'red'
plot3d(new_bug1,col = new_bug_colors)
#plot_selected_landmarks(new_bug1,new_bug1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

old_bug_colors = rep('white', dim(old_bug1$vb)[2])
old_bug_colors[old_bug1_vert] = 'blue'
plot3d(old_bug1,col = old_bug_colors)
#plot_selected_landmarks(old_bug1,old_bug1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

class_1_probs = read.csv('Data/new_aligned_shapesv3/V13_v4/gp1_spt.csv')
class_2_probs = read.csv('Data/new_aligned_shapesv3/V13_v4/gp2_spt.csv')

c7 = compute_roc_curve_teeth(data_dir1 = old_bug_dir, data_dir2 = new_bug_dir, gamma = 0.5,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                             rate_values = bug_bug_old$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = len,
                             directions = dirs,truncated = 50,class = 1,ball_radius = ball_radius, ball = ball)

c7_frame = data.frame(c7)
ggplot(c7_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "RATE ROC with Gamma = 0.5")

c8 = compute_roc_curve_teeth(data_dir1 = old_bug_dir, data_dir2 = new_bug_dir, gamma = 0.6,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                             rate_values = bug_bug_old$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = len,
                             directions = dirs,truncated = 50,class = 1,ball_radius = ball_radius, ball = ball)

c8_frame = data.frame(c8)
ggplot(c8_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "RATE ROC with Gamma = 0.6")
c9 = compute_roc_curve_teeth(data_dir1 = old_bug_dir, data_dir2 = new_bug_dir, gamma = 0.7,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                             rate_values = bug_bug_old$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = len,
                             directions = dirs,truncated = 50,class = 1,ball_radius = ball_radius, ball = ball)

c9_frame = data.frame(c9)
ggplot(c9_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "RATE ROC with Gamma = 0.7")


#### Heatmap for New Old ####



new_fruit1_vert = compute_selected_vertices_cones(dir = dirs, complex =new_fruit_1, rate_vals = fruit_fruit_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(fruit_fruit_old$Rate2)[1])),
                                                  cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

old_fruit1_vert = compute_selected_vertices_cones(dir = dirs, complex = old_fruit_1, rate_vals = fruit_fruit_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(fruit_fruit_old$Rate2)[1])),
                                                  cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

fruit_new_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = new_fruit_1,rate_vals = fruit_fruit_old$Rate2[,2],
                                                len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
fruit_old_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = old_fruit_1,rate_vals = fruit_fruit_old$Rate2[,2],
                                                len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

fruit_new_colors=colfunc(max(fruit_new_heat[,1]))[fruit_new_heat[,1]]
fruit_old_colors=colfunc(max(fruit_old_heat[,1]))[fruit_old_heat[,1]]



new_veg1_vert = compute_selected_vertices_cones(dir = dirs, complex =new_veg_1, rate_vals = veg_veg_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(veg_veg_old$Rate2)[1])),
                                                cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

old_veg1_vert = compute_selected_vertices_cones(dir = dirs, complex = old_veg_1, rate_vals = veg_veg_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(veg_veg_old$Rate2)[1])),
                                                cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

veg_new_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = new_veg_1,rate_vals = veg_veg_old$Rate2[,2],
                                              len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
veg_old_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = old_veg_1,rate_vals = veg_veg_old$Rate2[,2],
                                              len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

veg_new_colors=colfunc(max(veg_new_heat[,1]))[veg_new_heat[,1]]
veg_old_colors=colfunc(max(veg_old_heat[,1]))[veg_old_heat[,1]]



new_bug1_vert = compute_selected_vertices_cones(dir = dirs, complex =new_bug_1, rate_vals = bug_bug_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(bug_bug_old$Rate2)[1])),
                                                cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

old_bug1_vert = compute_selected_vertices_cones(dir = dirs, complex = old_bug_1, rate_vals = bug_bug_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(bug_bug_old$Rate2)[1])),
                                                cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

bug_new_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = new_bug_1,rate_vals = bug_bug_old$Rate2[,2],
                                              len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
bug_old_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = old_bug_1,rate_vals = bug_bug_old$Rate2[,2],
                                              len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

bug_new_colors=colfunc(max(bug_new_heat[,1]))[bug_new_heat[,1]]
bug_old_colors=colfunc(max(bug_old_heat[,1]))[bug_old_heat[,1]]


#### Plotting ####

mfrow3d(2,3)

plot3d(new_fruit1,col = fruit_new_colors, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(new_fruit1,new_fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)


plot3d(new_veg1,col = veg_new_colors,  specular="black", axes = FALSE)
#plot_selected_landmarks(new_veg1,new_veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)


plot3d(new_bug1,col = bug_new_colors,  specular="black", axes = FALSE)
#plot_selected_landmarks(new_bug1,new_bug1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(old_fruit1,col = fruit_old_colors, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(old_fruit1,old_fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(old_veg1,col = veg_old_colors, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(old_veg1,old_veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(old_bug1,col = bug_old_colors, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(old_bug1,old_bug1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

mfrow3d(5,5)
for (k in 1:25){
  file = vcgImport(new_fruit_files[k])
  plot3d(file, col = 'white', axes = FALSE)
  rgl.viewpoint(userMatrix = rotation_matrix)
}

sphere1 = vcgImport(old_bug_files[3])

sphere1$vb[1:3,] = sphere1$vb[1:3,]  * rnorm(dim(sphere1$vb)[2], mean = 1, sd = 0.035)
plot3d(sphere1, col = 'white')
