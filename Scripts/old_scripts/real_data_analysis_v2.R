setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
library(R.utils)
sourceDirectory('~/Documents/SINATRA/SINATRA_Pipeline_Branch/')
sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')

#library(svd)
set.seed(15)

#Specifying Directories and the Associated Files

#old_fruit_dir = 'Data/HDM/by_diet_scaled_aligned/Old_Frugivore/'
#old_veg_dir='Data/HDM/by_diet_scaled_aligned/Old_Follivore/'
#old_bug_dir='Data/HDM//by_diet_scaled_aligned/Old_Insectivore/'
#new_fruit_dir = 'Data/HDM/by_diet_scaled_aligned/New_Frugivore/'
#new_veg_dir='Data/HDM/by_diet_scaled_aligned/New_Follivore/'
#new_bug_dir='Data/HDM//by_diet_scaled_aligned/New_Insectivore/'

old_fruit_dir = 'Data/new_aligned_shapesv3/Old_Frugivore/'
old_veg_dir='Data/new_aligned_shapesv3/Old_Follivore/'
old_bug_dir='Data/new_aligned_shapesv3/Old_Insectivore/'
new_fruit_dir = 'Data/new_aligned_shapesv3/New_Frugivore/'
new_veg_dir='Data/new_aligned_shapesv3/New_Follivore/'
new_bug_dir='Data/new_aligned_shapesv3/New_Insectivore/'

old_fruit_files=list.files(path=old_fruit_dir,full.names = TRUE)
old_veg_files=list.files(path=old_veg_dir,full.names = TRUE)
old_bug_files=list.files(path=old_bug_dir,full.names = TRUE)

new_fruit_files=list.files(path=new_fruit_dir,full.names = TRUE)
new_veg_files=list.files(path=new_veg_dir,full.names = TRUE)
new_bug_files=list.files(path=new_bug_dir,full.names = TRUE)


files = c(old_fruit_files, old_veg_files, old_bug_files, new_fruit_files, new_veg_files, new_bug_files)
#Parameters for the Analysis
cap_radius = 0.25
num_cones = 25
directions_per_cone = 7
len = 100
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  cap_radius, directions_per_cone = directions_per_cone)
ball = FALSE
ball_radius = 0.5
ec_type = 'ECT'
#### Start Comparison #### 

fruit_veg_new=real_data_summary(dir1=new_fruit_dir,dir2 = new_veg_dir,direction=dirs,class1='Fruit', class2='Vegetable',
                                            radius=0,accuracy=FALSE,len = len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)

fruit_bug_new=real_data_summary(dir1 = new_fruit_dir,dir2 = new_bug_dir,class1='Fruit', direction=dirs,len=len,
                                            class2='Bug',accuracy=FALSE,radius=0, ball = ball, ball_radius = ball_radius, ec_type = ec_type)

veg_bug_new=real_data_summary(dir1 = new_veg_dir,dir2 = new_bug_dir,class1='Veg', direction=dirs,len=len,class2='Bug',
                                          radius=0,accuracy=FALSE, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
fruit_veg_old=real_data_summary(dir1=old_fruit_dir,dir2 = old_veg_dir,direction=dirs,class1='Fruit', class2='Vegetable',
                                            radius=0,accuracy=FALSE,len = len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)

fruit_bug_old=real_data_summary(dir1 = old_fruit_dir,dir2 = old_bug_dir,class1='Fruit', direction=dirs,len=len,
                                            class2='Bug',accuracy=FALSE,radius=0, ball = ball, ball_radius = ball_radius, ec_type = ec_type)

veg_bug_old=real_data_summary(dir1 = old_veg_dir,dir2 = old_bug_dir,class1='Veg', direction=dirs,len=len,class2='Bug',
                                          radius=0,accuracy=FALSE, ball = ball, ball_radius = ball_radius, ec_type = ec_type)

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

#### New Fruit Veg ####
ind = 4
rotation_matrix=matrix(c(0.99972576,0.02127766,0.00978078,0,0.01589701,-0.92330807,0.38373107,0,0.01719557,-0.38347024,-0.92339307,0,0,0,0,1),ncol=4,byrow=TRUE)
veg1 = vcgImport(new_veg_files[ind])
fruit1 = vcgImport(new_fruit_files[ind])
fruit_1 = process_off_file_v3(new_fruit_files[ind])
veg_1 = process_off_file_v3(new_veg_files[ind])

#Also finding the landmarks

mfrow3d(1,2)
fruit1_vert = compute_selected_vertices_cones(dir = dirs, complex =fruit_1, rate_vals = fruit_veg_new$Rate2[,2], len = len, threshold = ((directions_per_cone)/(dim(fruit_veg_new$Rate2)[1])),
                                              cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

veg1_vert = compute_selected_vertices_cones(dir = dirs, complex = veg_1, rate_vals = fruit_veg_new$Rate2[,2], len = len, threshold = ((directions_per_cone)/(dim(fruit_veg_new$Rate2)[1])),
                                            cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
fruit_colors = rep('white', dim(fruit1$vb)[2])
fruit_colors[fruit1_vert] = 'red'

plot3d(fruit1,col = fruit_colors)
#plot_selected_landmarks(fruit1,fruit1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

veg_colors = rep('white', dim(veg1$vb)[2])
veg_colors[veg1_vert] = 'blue'
plot3d(veg1,col = veg_colors)
#plot_selected_landmarks(veg1,veg1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

fruit_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = fruit_1,rate_vals = fruit_veg_new$Rate2[,2],
                                            len = len,cuts = 200,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
veg_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = veg_1,rate_vals = fruit_veg_new$Rate2[,2],
                                          len = len,cuts = 200,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

col_pal=c(color1,color2,color3, color4, color5)
#col_pal =c(color1,color2,color3)
colfunc <- colorRampPalette(col_pal)
fruit_colors=colfunc(max(fruit_heat[,1]))[fruit_heat[,1]]
veg_colors=colfunc(max(veg_heat[,1]))[veg_heat[,1]]

mfrow3d(1,2)
plot3d(fruit1,col = fruit_colors, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(fruit1,fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(veg1,col = veg_colors, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(veg1,veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)


#### New Fruit Bug ####

ind = 1
fruit1.5 = vcgImport(new_fruit_files[ind])
fruit_1.5 = process_off_file_v3(new_fruit_files[ind])
ind = 3
bug1.5 = vcgImport(new_bug_files[ind])
bug_1.5 = process_off_file_v3(new_bug_files[ind])

mfrow3d(1,2)
fruit1.5_vert = compute_selected_vertices_cones(dir = dirs, complex =fruit_1.5, rate_vals = fruit_bug_new$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(fruit_bug_new$Rate2)[1])),
                                                cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

bug1.5_vert = compute_selected_vertices_cones(dir = dirs, complex = bug_1.5, rate_vals = fruit_bug_new$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(fruit_bug_new$Rate2)[1])),
                                              cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
fruit_colors = rep('white', dim(fruit1.5$vb)[2])
fruit_colors[fruit1.5_vert] = 'red'
plot3d(fruit1.5,col = fruit_colors)
#plot_selected_landmarks(fruit1.5,fruit1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

bug_colors = rep('white', dim(bug1.5$vb)[2])
bug_colors[bug1.5_vert] = 'blue'
plot3d(bug1.5,col = bug_colors)
#plot_selected_landmarks(bug1.5,bug1.5_vert)
rgl.viewpoint(userMatrix = rotation_matrix)


fruit_heat1.5 =  reconstruct_vertices_on_shape(dir = dirs,complex = fruit_1.5,rate_vals = fruit_bug_new$Rate2[,2],
                                               len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
bug_heat1.5 =  reconstruct_vertices_on_shape(dir = dirs,complex = bug_1.5,rate_vals = fruit_bug_new$Rate2[,2],
                                             len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

fruit_colors1.5=colfunc(max(fruit_heat1.5[,1]))[fruit_heat1.5[,1]]
bug_colors1.5=colfunc(max(bug_heat1.5[,1]))[bug_heat1.5[,1]]

mfrow3d(1,2)
plot3d(fruit1.5,col = fruit_colors1.5, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(fruit1.5,fruit1.5_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(bug1.5,col = bug_colors1.5, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(bug1.5,bug1.5_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)


#### New Bug Vegs ####

ind = 3

bug2 = vcgImport(new_bug_files[ind])
veg2 = vcgImport(new_veg_files[ind])
veg_2 = process_off_file_v3(new_veg_files[ind])
bug_2 = process_off_file_v3(new_bug_files[ind])


mfrow3d(1,2)
veg2_vert = compute_selected_vertices_cones(dir = dirs, complex =veg_2, rate_vals = veg_bug_new$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(veg_bug_new$Rate2)[1])),
                                            cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

bug2_vert = compute_selected_vertices_cones(dir = dirs, complex = bug_2, rate_vals = veg_bug_new$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(veg_bug_new$Rate2)[1])),
                                            cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
veg_colors2 = rep('white', dim(veg2$vb)[2])
veg_colors2[veg2_vert] = 'red'
plot3d(veg2,col = veg_colors2)
#plot_selected_landmarks(veg2,veg2_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

bug_colors2 = rep('white', dim(bug2$vb)[2])
bug_colors2[bug2_vert] = 'blue'
plot3d(bug2,col = bug_colors2)
#plot_selected_landmarks(bug2,bug2_vert)
rgl.viewpoint(userMatrix = rotation_matrix)


veg_heat2 =  reconstruct_vertices_on_shape(dir = dirs,complex = veg_2,rate_vals = veg_bug_new$Rate2[,2],
                                           len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
bug_heat2 =  reconstruct_vertices_on_shape(dir = dirs,complex = bug_2,rate_vals = veg_bug_new$Rate2[,2],
                                           len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

veg_colors2=colfunc(max(veg_heat2[,1]))[veg_heat2[,1]]
bug_colors2=colfunc(max(bug_heat2[,1]))[bug_heat2[,1]]

mfrow3d(1,2)
plot3d(veg2,col = veg_colors2, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(veg2,veg2_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(bug2,col = bug_colors2, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(bug2,bug2_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

#### Comparisons Between Old World Monkeys ####

ind = 5
fruit3 = vcgImport(old_fruit_files[ind])
fruit_3 = process_off_file_v3(old_fruit_files[ind])

ind = 2
veg3 = vcgImport(old_veg_files[ind])
veg_3 = process_off_file_v3(old_veg_files[ind])

#Also finding the landmarks

mfrow3d(1,2)
fruit3_vert = compute_selected_vertices_cones(dir = dirs, complex =fruit_3, rate_vals = fruit_veg_old$Rate2[,2], len = len, threshold = ((directions_per_cone)/(dim(fruit_veg_old$Rate2)[1])),
                                              cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

veg3_vert = compute_selected_vertices_cones(dir = dirs, complex = veg_3, rate_vals = fruit_veg_old$Rate2[,2], len = len, threshold = ((directions_per_cone)/(dim(fruit_veg_old$Rate2)[1])),
                                            cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
fruit_colors3 = rep('white', dim(fruit3$vb)[2])
fruit_colors3[fruit3_vert] = 'red'

plot3d(fruit3,col = fruit_colors3)
#plot_selected_landmarks(fruit3,fruit3_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

veg_colors3 = rep('white', dim(veg3$vb)[2])
veg_colors3[veg3_vert] = 'blue'
plot3d(veg3,col = veg_colors3)
#plot_selected_landmarks(veg3,veg3_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

fruit_heat3 =  reconstruct_vertices_on_shape(dir = dirs,complex = fruit_3,rate_vals = fruit_veg_old$Rate2[,2],
                                             len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
veg_heat3 =  reconstruct_vertices_on_shape(dir = dirs,complex = veg_3,rate_vals = fruit_veg_old$Rate2[,2],
                                           len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

fruit_colors3=colfunc(max(fruit_heat3[,1]))[fruit_heat3[,1]]
veg_colors3=colfunc(max(veg_heat3[,1]))[veg_heat3[,1]]

mfrow3d(1,2)
plot3d(fruit3,col = fruit_colors3, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(fruit3,fruit3_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(veg3,col = veg_colors3, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(veg3,veg3_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

#### Old Fruit Bug ####
ind = 5
fruit4 = vcgImport(old_fruit_files[ind])
fruit_4 = process_off_file_v3(old_fruit_files[ind])
ind = 13
bug4 = vcgImport(old_bug_files[ind])
bug_4 = process_off_file_v3(old_bug_files[ind])

mfrow3d(1,2)
fruit4_vert = compute_selected_vertices_cones(dir = dirs, complex =fruit_4, rate_vals = fruit_bug_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(fruit_bug_old$Rate2)[1])),
                                              cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

bug4_vert = compute_selected_vertices_cones(dir = dirs, complex = bug_4, rate_vals = fruit_bug_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(fruit_bug_old$Rate2)[1])),
                                            cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
fruit_colors = rep('white', dim(fruit4$vb)[2])
fruit_colors[fruit4_vert] = 'red'
plot3d(fruit4,col = fruit_colors)
#plot_selected_landmarks(fruit4,fruit1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

bug_colors = rep('white', dim(bug4$vb)[2])
bug_colors[bug4_vert] = 'blue'
plot3d(bug4,col = bug_colors)
#plot_selected_landmarks(bug4,bug4_vert)
rgl.viewpoint(userMatrix = rotation_matrix)


fruit_heat4 =  reconstruct_vertices_on_shape(dir = dirs,complex = fruit_4,rate_vals = fruit_bug_old$Rate2[,2],
                                             len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
bug_heat4 =  reconstruct_vertices_on_shape(dir = dirs,complex = bug_4,rate_vals = fruit_bug_old$Rate2[,2],
                                           len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

fruit_colors4=colfunc(max(fruit_heat4[,1]))[fruit_heat4[,1]]
bug_colors4=colfunc(max(bug_heat4[,1]))[bug_heat4[,1]]

mfrow3d(1,2)
plot3d(fruit4,col = fruit_colors4, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(fruit4,fruit4_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(bug4,col = bug_colors4, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(bug4,bug4_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



#### New Bug Vegs ####

ind = 2

veg5 = vcgImport(old_veg_files[ind])
veg_5 = process_off_file_v3(old_veg_files[ind])
ind = 7
bug5 = vcgImport(old_bug_files[ind])
bug_5 = process_off_file_v3(old_bug_files[ind])


mfrow3d(1,2)
veg5_vert = compute_selected_vertices_cones(dir = dirs, complex =veg_5, rate_vals = veg_bug_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(veg_bug_old$Rate2)[1])),
                                            cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

bug5_vert = compute_selected_vertices_cones(dir = dirs, complex = bug_5, rate_vals = veg_bug_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(veg_bug_old$Rate2)[1])),
                                            cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
veg_colors5 = rep('white', dim(veg5$vb)[2])
veg_colors5[veg5_vert] = 'red'
plot3d(veg5,col = veg_colors5)
#plot_selected_landmarks(veg5,veg5_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

bug_colors5 = rep('white', dim(bug5$vb)[2])
bug_colors5[bug5_vert] = 'blue'
plot3d(bug5,col = bug_colors5)
#plot_selected_landmarks(bug5,bug5_vert)
rgl.viewpoint(userMatrix = rotation_matrix)


veg_heat5 =  reconstruct_vertices_on_shape(dir = dirs,complex = veg_5,rate_vals = veg_bug_old$Rate2[,2],
                                           len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
bug_heat5 =  reconstruct_vertices_on_shape(dir = dirs,complex = bug_5,rate_vals = veg_bug_old$Rate2[,2],
                                           len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

veg_colors5=colfunc(max(veg_heat5[,1]))[veg_heat5[,1]]
bug_colors5=colfunc(max(bug_heat5[,1]))[bug_heat5[,1]]

mfrow3d(1,2)
plot3d(veg5,col = veg_colors5, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(veg5,veg5_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(bug5,col = bug_colors5, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(bug5,bug5_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

#### Old Fruit Fruit ####
ind = 7
new_fruit1 = vcgImport(new_fruit_files[ind])
new_fruit_1 = process_off_file_v3(new_fruit_files[ind])

ind = 5
old_fruit1 = vcgImport(old_fruit_files[ind])
old_fruit_1 = process_off_file_v3(old_fruit_files[ind])

mfrow3d(1,2)
new_fruit1_vert = compute_selected_vertices_cones(dir = dirs, complex =new_fruit_1, rate_vals = fruit_fruit_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(fruit_fruit_old$Rate2)[1])),
                                                  cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

old_fruit1_vert = compute_selected_vertices_cones(dir = dirs, complex = old_fruit_1, rate_vals = fruit_fruit_old$Rate2[,2], len = len, threshold = (directions_per_cone/(dim(fruit_fruit_old$Rate2)[1])),
                                                  cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
new_fruit_colors = rep('white', dim(new_fruit1$vb)[2])
new_fruit_colors[new_fruit1_vert] = 'red'
plot3d(new_fruit1,col = new_fruit_colors)
#plot_selected_landmarks(new_fruit1,new_fruit1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

old_fruit_colors = rep('white', dim(old_fruit1$vb)[2])
old_fruit_colors[old_fruit1_vert] = 'blue'
plot3d(old_fruit1,col = old_fruit_colors)
#plot_selected_landmarks(old_fruit1,old_fruit1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)


fruit_new_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = new_fruit_1,rate_vals = fruit_fruit_old$Rate2[,2],
                                                len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
fruit_old_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = old_fruit_1,rate_vals = fruit_fruit_old$Rate2[,2],
                                                len = len,cuts = 100,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

fruit_new_colors=colfunc(max(fruit_new_heat[,1]))[fruit_new_heat[,1]]
fruit_old_colors=colfunc(max(fruit_old_heat[,1]))[fruit_old_heat[,1]]


mfrow3d(1,2)
plot3d(new_fruit1,col = fruit_new_colors, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(new_fruit1,new_fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(old_fruit1,col = fruit_old_colors, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(old_fruit1,old_fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

#### Old Veg New Veg ####
ind = 6
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
plot3d(old_veg1,col = old_veg_colors, axes = FALSE)
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



#### Old Bug New Bug ####
ind = 2
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


plot3d(new_veg1,col = veg_new_colors, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(new_veg1,new_veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)


plot3d(new_bug1,col = bug_new_colors, back="lines", specular="black", axes = FALSE)
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



#### Outputs for Tingran ###



process_tingran_files(new_fruit_dir,comp1 = fruit_veg_new, comp2 = fruit_bug_new, comp3 = fruit_fruit_old, 
                      name1 = '_new_fruit_veg',name2 = '_new_fruit_bug',name3 = '_new_old_fruit') 


process_tingran_files(new_veg_dir,comp1 = fruit_veg_new, comp2 = veg_bug_new, comp3 = veg_veg_old, 
                      name1 = '_new_fruit_veg',name2 = '_new_veg_bug',name3 = '_new_old_veg') 

process_tingran_files(new_bug_dir,comp1 = fruit_bug_new, comp2 = veg_bug_new, comp3 = bug_bug_old, 
                      name1 = '_new_fruit_bug',name2 = '_new_veg_bug',name3 = '_new_old_bug') 



process_tingran_files(old_fruit_dir,comp1 = fruit_veg_old, comp2 = fruit_bug_old, comp3 = fruit_fruit_old, 
                      name1 = '_old_fruit_veg',name2 = '_old_fruit_bug',name3 = '_new_old_fruit') 


process_tingran_files(old_veg_dir,comp1 = fruit_veg_old, comp2 = veg_bug_old, comp3 = veg_veg_old, 
                      name1 = '_old_fruit_veg',name2 = '_old_veg_bug',name3 = '_new_old_veg') 

process_tingran_files(old_bug_dir,comp1 = fruit_bug_old, comp2 = veg_bug_old, comp3 = bug_bug_old, 
                      name1 = '_old_fruit_bug',name2 = '_old_veg_bug',name3 = '_new_old_bug') 




