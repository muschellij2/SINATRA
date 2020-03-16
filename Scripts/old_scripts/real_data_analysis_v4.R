setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
library(R.utils)
sourceDirectory('~/Documents/SINATRA/SINATRA_Pipeline_Branch/')
sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')

#library(svd)
#set.seed(15)
#set.seed(35)
set.seed(15)

#Specifying Directories and the Associated Files

#old_fruit_dir = 'Data/HDM/by_diet_scaled_aligned/Old_Frugivore/'
#old_veg_dir='Data/HDM/by_diet_scaled_aligned/Old_Follivore/'
#old_bug_dir='Data/HDM//by_diet_scaled_aligned/Old_Insectivore/'
#new_fruit_dir = 'Data/HDM/by_diet_scaled_aligned/New_Frugivore/'
#new_veg_dir='Data/HDM/by_diet_scaled_aligned/New_Follivore/'
#new_bug_dir='Data/HDM//by_diet_scaled_aligned/New_Insectivore/'


old_bug_dir='Data/new_aligned_shapesv3/V13_v2/gp1'
new_bug_dir='Data/new_aligned_shapesv3/V13_v2/gp2'
#old_fruit_dir = 'Data/new_aligned_shapesv3/V13_v5/gp1'
#new_fruit_dir = 'Data/new_aligned_shapesv3/V13_v5/gp2'
old_veg_dir='Data/new_aligned_shapesv3/V13_v2/gp1'
new_veg_dir='Data/new_aligned_shapesv3/V13_v2/gp2'
old_fruit_dir = 'Data/new_aligned_shapesv3/V13_v2/gp1'
new_fruit_dir = 'Data/new_aligned_shapesv3/V13_v2/gp2'

old_fruit_files=list.files(path=old_fruit_dir,full.names = TRUE)
old_veg_files=list.files(path=old_veg_dir,full.names = TRUE)
old_bug_files=list.files(path=old_bug_dir,full.names = TRUE)

new_fruit_files=list.files(path=new_fruit_dir,full.names = TRUE)
new_veg_files=list.files(path=new_veg_dir,full.names = TRUE)
new_bug_files=list.files(path=new_bug_dir,full.names = TRUE)

veg_pset = list(num_cones = 35, cap_radius = 0.15, len = 100, directions_per_cone = 5,
                dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))
fruit_pset = list(num_cones = 35, cap_radius = 0.25, len = 100, directions_per_cone = 5,
                dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))
bug_pset = list(num_cones = 10, cap_radius = 0.15, len = 75, directions_per_cone = 5,
                dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))
bug2_pset = list(num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
                dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))
test_pset = list(num_cones = 10, cap_radius = 0.25, len = 25, directions_per_cone = 5,
                dirs = generate_equidistributed_cones(num_directions = 10, cap_radius =  0.15, directions_per_cone = 5))
files = c(old_fruit_files, old_veg_files, old_bug_files, new_fruit_files, new_veg_files, new_bug_files)
#Parameters for the Analysis
ball = TRUE
ball_radius = 0.5
ec_type = 'ECT'
#### Start Comparison #### 

# 35, 0.25 cap radius
fruit_fruit_old=real_data_summary(dir1=new_veg_dir,dir2 = old_veg_dir,direction=fruit_pset$dirs,class1='Fruit', class2='Fruit',
                                              radius=0,accuracy=FALSE,len = fruit_pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)

veg_veg_old=real_data_summary(dir1=new_veg_dir,dir2 = old_veg_dir,direction=veg_pset$dirs,class1='Vegetable', class2='Vegetable',
                                          radius=0,accuracy=FALSE,len = veg_pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
#35, 75 length 0.25 cap radius
bug_bug_old=real_data_summary(dir1=new_veg_dir,dir2 = old_veg_dir,direction=bug_pset$dirs,class1='Bug', class2='Bug',
                                          radius=0,accuracy=FALSE,len = bug_pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
#0.15 , 75 length cap radius
bug_bug_old2=real_data_summary(dir1=new_veg_dir,dir2 = old_veg_dir,direction=bug2_pset$dirs,class1='Bug', class2='Bug',
                                          radius=0,accuracy=FALSE,len = bug2_pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)

test_rate=real_data_summary(dir1=new_veg_dir,dir2 = old_veg_dir,direction=bug2_pset$dirs,class1='Bug', class2='Bug',
                                          radius=0,accuracy=FALSE,len = bug2_pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)

color1='blue'
color2 = 'blue'
color3 = 'lightblue'
color4='palegreen'
color5='orangered'
col_pal=c(color1,color1,color2,color2,color3,color3,color4,color5)
colfunc <- colorRampPalette(col_pal)
rotation_matrix=matrix(c(0.99972576,0.02127766,0.00978078,0,0.01589701,-0.92330807,0.38373107,0,0.01719557,-0.38347024,-0.92339307,0,0,0,0,1),ncol=4,byrow=TRUE)

#### Old Fruit Fruit ####
class_1_probs = (read.csv('Data/new_aligned_shapesv3/V13_v2/gp1_spt.csv',header = FALSE))
class_2_probs = (read.csv('Data/new_aligned_shapesv3/V13_v2/gp2_spt.csv', header = FALSE))
ind = 3
new_fruit1 = vcgImport(new_fruit_files[ind])
new_fruit_1 = process_off_file_v3(new_fruit_files[ind])

ind = 3
old_fruit1 = vcgImport(old_fruit_files[ind])
old_fruit_1 = process_off_file_v3(old_fruit_files[ind])

mfrow3d(1,2)
new_fruit1_vert = compute_selected_vertices_cones(dir = fruit_pset$dirs, complex =new_fruit_1, rate_vals = fruit_fruit_old$Rate2[,2], len = len, threshold = (directions_per_cone/(1*dim(fruit_fruit_old$Rate2)[1])),
                                                  cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball )

old_fruit1_vert = compute_selected_vertices_cones(dir = fruit_pset$dirs, complex = old_fruit_1, rate_vals = fruit_fruit_old$Rate2[,2], len = len, threshold = (directions_per_cone/(1*dim(fruit_fruit_old$Rate2)[1])),
                                                  cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
new_fruit_colors = rep('white', dim(new_fruit1$vb)[2])
new_fruit_colors[new_fruit1_vert] = 'red'
v1 = t(new_fruit1$vb[-4,])
v2 = t(old_fruit1$vb[-4,])
plot3d(new_fruit1,col = new_fruit_colors)
rgl.points(v1[( which(class_2_probs > 0.05)),], size = 3, col = 'green')
#plot_selected_landmarks(new_fruit1,new_fruit1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

old_fruit_colors = rep('white', dim(old_fruit1$vb)[2])
old_fruit_colors[old_fruit1_vert] = 'blue'
plot3d(old_fruit1,col = old_fruit_colors)
rgl.points(v2[which(class_1_probs > 0.05),], size = 3, col = 'green')
#plot_selected_landmarks(old_fruit1,old_fruit1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)


fruit_new_heat =  reconstruct_vertices_on_shape(dir = fruit_pset$dirs,complex = new_fruit_1,rate_vals = fruit_fruit_old$Rate2[,2],
                                                len = fruit_pset$len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
fruit_old_heat =  reconstruct_vertices_on_shape(dir = fruit_pset$dirs,complex = old_fruit_1,rate_vals = fruit_fruit_old$Rate2[,2],
                                                len = fruit_pset$len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

col_pal=c(color1,'lightgreen',color5)
colfunc <- colorRampPalette(col_pal)
fruit_new_colors=colfunc(max(fruit_new_heat[,1])- min(fruit_new_heat[,1]))[(fruit_new_heat[,1] - min(fruit_new_heat[,1]))]
fruit_old_colors=colfunc(max(fruit_old_heat[,1])- min(fruit_new_heat[,1]))[(fruit_old_heat[,1] -min(fruit_old_heat[,1]))]

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
class_1_probs = (read.csv('data/new_aligned_shapesv3/v13_v2/gp1_spt.csv',header = FALSE))
class_2_probs = (read.csv('Data/new_aligned_shapesv3/V13_v2/gp2_spt.csv', header = FALSE))

c1 = compute_roc_curve_teeth(data_dir1 = old_fruit_dir, data_dir2 = new_fruit_dir, gamma = 0.5,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                   rate_values = fruit_fruit_old$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = len,
                                   directions = dirs,truncated = 30,ball_radius = ball_radius, ball = ball)

c1_frame = data.frame(c1)
ggplot(c1_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "V2 35 Cones, 0.25 Cap Radius RATE ROC with Gamma = 0.5")
c2 = compute_roc_curve_teeth(data_dir1 = old_fruit_dir, data_dir2 = new_fruit_dir, gamma = 0.05,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                   rate_values = fruit_fruit_old$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = fruit_pset$len,
                                   directions = fruit_pset$dirs,truncated = 30,ball_radius = ball_radius, ball = ball)

c2_frame = data.frame(c2)
ggplot(c2_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "V2 35 Cones, 0.25 Cap Radius RATE ROC with Gamma = 0.05")


mfrow3d(2,3)
og_mesh = vcgImport(file = 'Data/new_aligned_shapesv3/Old_Follivore/clean_V13_sas.off')
plot3d(og_mesh, col = 'white', axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
new_fruit_colors = rep('white', dim(new_fruit1$vb)[2])
new_fruit_colors[(which(class_2_probs > 0.05))] = 'green'
plot3d(new_fruit1,col = new_fruit_colors,specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
plot3d(new_fruit1,col = fruit_new_colors, specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(og_mesh, col = 'white',specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
old_fruit_colors = rep('white', dim(old_fruit1$vb)[2])
old_fruit_colors[(which(class_1_probs > 0.05))] = 'green'

plot3d(old_fruit1,col = old_fruit_colors, axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
#plot_selected_landmarks(new_fruit1,new_fruit1_vert,num_landmarks = 20)

plot3d(old_fruit1,col = fruit_old_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(old_fruit1,old_fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)


#### Old Veg New Veg ####
class_1_probs = read.csv('Data/new_aligned_shapesv3/V13_v2/gp1_spt.csv', header = FALSE)
class_2_probs = read.csv('Data/new_aligned_shapesv3/V13_v2/gp2_spt.csv', header = FALSE)
ind = 3
old_veg1 = vcgImport(old_veg_files[ind])
new_veg1 = vcgImport(new_veg_files[ind])
new_veg_1 = process_off_file_v3(new_veg_files[ind])
old_veg_1 = process_off_file_v3(old_veg_files[ind])
v1 = t(new_veg1$vb[-4,])
v2 = t(old_veg1$vb[-4,])

mfrow3d(2,2)
new_veg_colors = rep('white', dim(new_veg1$vb)[2])
new_veg_colors[(which(class_2_probs > 0.7))] = 'green'
plot3d(new_veg1,col = new_veg_colors, axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)

old_veg_colors = rep('white', dim(old_veg1$vb)[2])
old_veg_colors[(which(class_1_probs > 0.7))] = 'green'

plot3d(old_veg1,col = old_veg_colors, axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
plot3d(new_veg1,col = veg_new_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(new_veg1,new_veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(old_veg1,col = veg_old_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(old_veg1,old_veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

open3d()
mfrow3d(1,2)
#plot_selected_landmarks(old_fruit1,old_fruit1_vert,num_landmarks = 20)
new_veg1_vert = compute_selected_vertices_cones(dir = veg_pset$dirs, complex =new_veg_1, rate_vals = veg_veg_old$Rate2[,2], len = veg_pset$len, threshold = (veg_pset$directions_per_cone/(2*dim(veg_veg_old$Rate2)[1])),
                                                cone_size = veg_pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

old_veg1_vert = compute_selected_vertices_cones(dir = veg_pset$dirs, complex = old_veg_1, rate_vals = veg_veg_old$Rate2[,2], len = veg_pset$len, threshold = (veg_pset$directions_per_cone/(2*dim(veg_veg_old$Rate2)[1])),
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
                                              len = veg_pset$len,cuts = 100,cone_size = veg_pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
veg_old_heat =  reconstruct_vertices_on_shape(dir = veg_pset$dirs,complex = old_veg_1,rate_vals = veg_veg_old$Rate2[,2],
                                              len = veg_pset$len,cuts = 100,cone_size = veg_pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

#veg_new_colors=colfunc(max(veg_new_heat[,1]))[veg_new_heat[,1]]
#veg_old_colors=colfunc(max(veg_old_heat[,1]))[veg_old_heat[,1]]

col_pal=c(color1,'lightblue','lightgreen',color5)
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


class_1_probs = read.csv('Data/new_aligned_shapesv3/V13_v2/gp1_spt.csv', header = FALSE)
class_2_probs = read.csv('Data/new_aligned_shapesv3/V13_v2/gp2_spt.csv', header = FALSE)


c9 = compute_roc_curve_teeth(data_dir1 = old_veg_dir, data_dir2 = new_veg_dir, gamma = 0.5,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                             rate_values = veg_veg_old$Rate2[,2],directions_per_cone = veg_pset$directions_per_cone,curve_length = veg_pset$len,
                             directions = veg_pset$dirs,truncated = 100,ball_radius = ball_radius, ball = ball, radius = 2)

c9_frame = data.frame(c9)
ggplot(c9_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "V2 +/- 2 35 cones 0.15 cap radius 100 curve length RATE ROC with Gamma = 0.5")

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






#### Veg 40 dirs ####
class_1_probs = (read.csv('Data/new_aligned_shapesv3/V13_v2/gp1_spt.csv',header = FALSE))
class_2_probs = (read.csv('Data/new_aligned_shapesv3/V13_v2/gp2_spt.csv', header = FALSE))
ind = 3
new_bug1 = vcgImport(new_bug_files[ind])
new_bug_1 = process_off_file_v3(new_bug_files[ind])

ind = 3
old_bug1 = vcgImport(old_bug_files[ind])
old_bug_1 = process_off_file_v3(old_bug_files[ind])

mfrow3d(1,2)
new_bug1_vert = compute_selected_vertices_cones(dir = bug_pset$dirs, complex =new_bug_1, rate_vals = bug_bug_old$Rate2[,2], len = len, threshold = (directions_per_cone/(1*dim(bug_bug_old$Rate2)[1])),
                                                cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

old_bug1_vert = compute_selected_vertices_cones(dir = bug_pset$dirs, complex = old_bug_1, rate_vals = bug_bug_old$Rate2[,2], len = len, threshold = (directions_per_cone/(1*dim(bug_bug_old$Rate2)[1])),
                                                cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
new_bug_colors = rep('white', dim(new_bug1$vb)[2])
new_bug_colors[new_bug1_vert] = 'red'
v1 = t(new_bug1$vb[-4,])
v2 = t(old_bug1$vb[-4,])
plot3d(new_bug1,col = new_bug_colors)
rgl.points(v1[( which(class_2_probs > 0.05)),], size = 3, col = 'green')
#plot_selected_landmarks(new_bug1,new_bug1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

old_bug_colors = rep('white', dim(old_bug1$vb)[2])
old_bug_colors[old_bug1_vert] = 'blue'
plot3d(old_bug1,col = old_bug_colors)
rgl.points(v2[which(class_1_probs > 0.05),], size = 3, col = 'green')
#plot_selected_landmarks(old_bug1,old_bug1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)


bug_new_heat =  reconstruct_vertices_on_shape(dir = bug_pset$dirs,complex = new_bug_1,rate_vals = bug_bug_old$Rate2[,2],
                                              len = bug_pset$len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
bug_old_heat =  reconstruct_vertices_on_shape(dir = bug_pset$dirs,complex = old_bug_1,rate_vals = bug_bug_old$Rate2[,2],
                                              len = bug_pset$len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

#col_pal=c(color1,'lightgreen','orangered',color5)
colfunc <- colorRampPalette(col_pal)
bug_new_colors=colfunc(max(bug_new_heat[,1])- min(bug_new_heat[,1]))[(bug_new_heat[,1] - min(bug_new_heat[,1]))]
bug_old_colors=colfunc(max(bug_old_heat[,1])- min(bug_old_heat[,1]))[(bug_old_heat[,1] -min(bug_old_heat[,1]))]

#bug_new_colors=colfunc(max(bug_new_heat[,1]))[bug_new_heat[,1]]
#bug_old_colors=colfunc(max(bug_old_heat[,1]))[bug_old_heat[,1]]

mfrow3d(1,2)
plot3d(new_bug1,col = bug_new_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(new_bug1,new_bug1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(old_bug1,col = bug_old_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(old_bug1,old_bug1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

#class_1_probs = read.csv('Data/new_aligned_shapesv3/V13_v3/gp1_spt.csv')
#class_2_probs = read.csv('Data/new_aligned_shapesv3/V13_v3/gp2_spt.csv')
class_1_probs = (read.csv('data/new_aligned_shapesv3/v13_v2/gp1_spt.csv',header = FALSE))
class_2_probs = (read.csv('Data/new_aligned_shapesv3/V13_v2/gp2_spt.csv', header = FALSE))

c2 = compute_roc_curve_teeth(data_dir1 = old_bug_dir, data_dir2 = new_bug_dir, gamma = 0.5,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                             rate_values = bug_bug_old$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = bug_pset$len,
                             directions = bug_pset$dirs,truncated = 30,ball_radius = ball_radius, ball = ball)

c2_frame = data.frame(c2)
ggplot(c2_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "V2, 35, 75 length, 0.25 radius cones RATE ROC with Gamma = 0.5")

mfrow3d(2,3)
og_mesh = vcgImport(file = 'Data/new_aligned_shapesv3/Old_Follivore/clean_V13_sas.off')
plot3d(og_mesh, col = 'white', axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
new_bug_colors = rep('white', dim(new_bug1$vb)[2])
new_bug_colors[(which(class_2_probs > 0.05))] = 'green'
plot3d(new_bug1,col = new_bug_colors,specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
plot3d(new_bug1,col = bug_new_colors, specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(og_mesh, col = 'white',specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
old_bug_colors = rep('white', dim(old_bug1$vb)[2])
old_bug_colors[(which(class_1_probs > 0.05))] = 'green'

plot3d(old_bug1,col = old_bug_colors, axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
#plot_selected_landmarks(new_bug1,new_bug1_vert,num_landmarks = 20)

plot3d(old_bug1,col = bug_old_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(old_bug1,old_bug1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

#### 0.25 cap radius ####
class_1_probs = (read.csv('Data/new_aligned_shapesv3/V13_v2/gp1_spt.csv',header = FALSE))
class_2_probs = (read.csv('Data/new_aligned_shapesv3/V13_v2/gp2_spt.csv', header = FALSE))
ind = 3
new_test1 = vcgImport(new_bug_files[ind])
new_test_1 = process_off_file_v3(new_bug_files[ind])

ind = 3
old_test1 = vcgImport(old_bug_files[ind])
old_test_1 = process_off_file_v3(old_bug_files[ind])

mfrow3d(1,2)
new_test1_vert = compute_selected_vertices_cones(dir = bug2_pset$dirs, complex =new_test_1, rate_vals = bug_bug_old2$Rate2[,2], len = bug2_pset$len, threshold = (directions_per_cone/(1*dim(bug_bug_old2$Rate2)[1])),
                                                 cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

old_test1_vert = compute_selected_vertices_cones(dir = bug2_pset$dirs, complex = old_test_1, rate_vals = bug_bug_old2$Rate2[,2], len = bug2_pset$len, threshold = (directions_per_cone/(1*dim(bug_bug_old2$Rate2)[1])),
                                                 cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
new_test_colors = rep('white', dim(new_test1$vb)[2])
new_test_colors[new_test1_vert] = 'red'
v1 = t(new_test1$vb[-4,])
v2 = t(old_test1$vb[-4,])
plot3d(new_test1,col = new_test_colors)
rgl.points(v1[( which(class_2_probs > 0.05)),], size = 3, col = 'green')
#plot_selected_landmarks(new_test1,new_test1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

old_test_colors = rep('white', dim(old_test1$vb)[2])
old_test_colors[old_test1_vert] = 'blue'
plot3d(old_test1,col = old_test_colors)
rgl.points(v2[which(class_1_probs > 0.05),], size = 3, col = 'green')
#plot_selected_landmarks(old_test1,old_test1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)


test_new_heat =  reconstruct_vertices_on_shape(dir = bug2_pset$dirs,complex = new_test_1,rate_vals = bug_bug_old2$Rate2[,2],
                                               len = bug2_pset$len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
test_old_heat =  reconstruct_vertices_on_shape(dir = bug2_pset$dirs,complex = old_test_1,rate_vals = bug_bug_old2$Rate2[,2],
                                               len = bug2_pset$len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

#col_pal=c(color1,'lightgreen','orangered',color5)
colfunc <- colorRampPalette(col_pal)
test_new_colors=colfunc(max(test_new_heat[,1]- min(test_new_heat[,1])))[(test_new_heat[,1] - min(test_new_heat[,1]))]
test_old_colors=colfunc(max(test_old_heat[,1]- min(test_old_heat[,1])))[(test_old_heat[,1] -min(test_old_heat[,1]))]


mfrow3d(1,2)
plot3d(new_test1,col = test_new_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(new_test1,new_test1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(old_test1,col = test_old_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(old_test1,old_test1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

#class_1_probs = read.csv('Data/new_aligned_shapesv3/V13_v3/gp1_spt.csv')
#class_2_probs = read.csv('Data/new_aligned_shapesv3/V13_v3/gp2_spt.csv')
class_1_probs = (read.csv('data/new_aligned_shapesv3/v13_v2/gp1_spt.csv',header = FALSE))
class_2_probs = (read.csv('Data/new_aligned_shapesv3/V13_v2/gp2_spt.csv', header = FALSE))

c3 = compute_roc_curve_teeth(data_dir1 = old_bug_dir, data_dir2 = new_bug_dir, gamma = 0.5,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                             rate_values = bug_bug_old2$Rate2[,2],directions_per_cone = directions_per_cone,curve_length = bug2_pset$len,
                             directions = bug2_pset$dirs,truncated = 30,ball_radius = ball_radius, ball = ball)

c3_frame = data.frame(c3)
ggplot(c3_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "V2, 35,0.15, curve_length 75 cap radius cones RATE ROC with Gamma = 0.5")
mfrow3d(2,3)
og_mesh = vcgImport(file = 'Data/new_aligned_shapesv3/Old_Follivore/clean_V13_sas.off')
plot3d(og_mesh, col = 'white', axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
new_test_colors = rep('white', dim(new_test1$vb)[2])
new_test_colors[(which(class_2_probs > 0.5))] = 'green'
plot3d(new_test1,col = new_test_colors,specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
plot3d(new_test1,col = test_new_colors, specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(og_mesh, col = 'white',specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
old_test_colors = rep('white', dim(old_test1$vb)[2])
old_test_colors[(which(class_1_probs > 0.5))] = 'green'

plot3d(old_test1,col = old_test_colors, axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)
#plot_selected_landmarks(new_test1,new_test1_vert,num_landmarks = 20)

plot3d(old_test1,col = test_old_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(old_test1,old_test1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)


top_dirs = c(1,3)
find_dirs = matrix(bug_pset$dirs[1,],ncol = 3)

d =  reconstruct_vertices_on_shape(dir = find_dirs,complex = new_bug_1,rate_vals = bug_bug_old$Rate2[,2],
                                   len = bug_pset$len,cuts = 100,cone_size = 1,ball_radius = ball_radius, ball = TRUE)

test_colors=colfunc(max(d[,1])- min(d[,1]))[(d[,1] -min(d[,1]))]
plot3d(new_bug1,col = test_colors, specular="black")

sp = vcgSphere(4)
#sp$vb = sp$vb/2
sp = scalemesh(sp,0.5 )
new_sp = convert_off_file(sp)

hist(new_sp$Vertices[,1:3]%*%find_dirs[1,])

d1 =  reconstruct_vertices_on_shape(dir = find_dirs,complex = new_sp,rate_vals = bug_bug_old$Rate2[,2],
                                   len = 75,cuts = 10,cone_size = 1,ball_radius = ball_radius, ball = TRUE)


test_colors1=colfunc(max(d1[,1])- min(d1[,1]))[(d1[,1] -min(d1[,1]))]

plot3d(sp,col = test_colors1, specular="black")


#### Debugging ####
i=5
projections <- new_sp$Vertices[,1:3]%*%find_dirs[i,]
l = 20
buckets <- seq(-ball_radius,ball_radius,length.out = l+1)

#bucket these projections into curve_length number of groups; could have also solved this with the cut function
step_length <- (max(buckets) - min(buckets))/l
#Replace projections by buckets
projection_buckets <- apply((projections - min(buckets))/step_length,1, function(float) as.integer(float)) + (l)*(i-1)
#print(step_l)
projection_buckets=projection_buckets+1
print(paste(min(projection_buckets), max(projection_buckets)))

d = which(projection_buckets %in% (c(1,2,10,20)+(i - 1)*l))

sphere_colors=colfunc(max(projection_buckets))[projection_buckets]
sp_col = rep('white', dim(sp$vb)[2])
sp_col[d] = 'red'
plot3d(sp,col = sphere_colors, specular = 'black', axes = FALSE)
plot3d(sp,col = sp_col, specular = 'black', axes = FALSE)


#### Debug part 2 ####

buckets2 <- seq(-ball_radius,ball_radius,length.out = l+1)

# map vertex projection to the feature index
projection_bucket2 <- cut(projections, buckets2, labels = FALSE)

# update index to reflect rate values
projection_buckets2 <- projection_bucket2 + (i - 1)*l

print(paste(min(projection_buckets2), max(projection_buckets2)))

d2 = which(projection_buckets2 %in% (c(1,2,10,20)+(i - 1)*l))

sphere_colors2=colfunc(max(projection_buckets2))[projection_buckets2]
sp_col2 = rep('white', dim(sp$vb)[2])
sp_col2[d2] = 'red'
plot3d(sp,col = sphere_colors2, specular = 'black', axes = FALSE)
plot3d(sp,col = sp_col2, specular = 'black', axes = FALSE)


#### Debug ####

g = projection_buckets == projection_buckets2
b = which(g == FALSE)
