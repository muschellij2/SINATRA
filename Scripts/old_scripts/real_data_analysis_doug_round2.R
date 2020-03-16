setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
library(R.utils)
sourceDirectory('~/Documents/SINATRA/SINATRA_Pipeline_Branch/')
sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')
set.seed(55)
load('new_world_monkey_comparison.Rdata')
#Specifying Directories and the Associated Files

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


#Parameters for the Analysis
cap_radius = 0.15
num_cones = 35
directions_per_cone = 5
len = 100 
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  cap_radius, directions_per_cone = directions_per_cone)
ball = TRUE
ball_radius = 0.5
ec_type = 'ECT'
#### Start Comparison #### 

pset1 = list(dir1 = new_fruit_dir,dir2 = new_veg_dir,
             num_cones = 35, cap_radius = 0.15, len = 100, directions_per_cone = 5,
                dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))
pset2 =  list(dir1 = new_fruit_dir,dir2 = new_bug_dir,
              num_cones = 35, cap_radius = 0.15, len = 100, directions_per_cone = 5,
              dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))
pset3 = list(dir1 = new_veg_dir,dir2 = new_bug_dir,
             num_cones = 35, cap_radius = 0.15, len = 100, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))

comp1 = real_data_summary(dir1=pset1$dir1,dir2 = pset1$dir2,direction=pset1$dirs,class1=pset1$dir1, class2=pset1$dir1,
                                    radius=0,accuracy=FALSE,len = pset1$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
comp2 = real_data_summary(dir1=pset2$dir1,dir2 = pset2$dir2,direction=pset2$dirs,class1=pset2$dir1, class2=pset2$dir1,
                          radius=0,accuracy=FALSE,len = pset2$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
comp3 = real_data_summary(dir1=pset3$dir1,dir2 = pset3$dir2,direction=pset3$dirs,class1=pset3$dir1, class2=pset3$dir1,
                          radius=0,accuracy=FALSE,len = pset3$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)

color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)
#col_pal =c(color1,color2,color3)
colfunc <- colorRampPalette(col_pal)
rotation_matrix=matrix(c(0.8065218,0.5911149,0.01028626,0,0.5583186,-0.7558216,0.34207344,0,0.1944301,-0.2816325,-0.93961692,0,0,0,0,1),ncol=4,byrow=TRUE)

#### New Fruit Veg ####
ind = 13
c1_mesh1 = vcgImport(list.files(pset1$dir1,full.names = TRUE)[ind])
c1_mesh2 = vcgImport(list.files(pset1$dir2,full.names = TRUE)[ind])
c1_mesh_1 = process_off_file_v3(list.files(pset1$dir1,full.names = TRUE)[ind])
c1_mesh_2 = process_off_file_v3(list.files(pset1$dir2,full.names = TRUE)[ind])


mfrow3d(1,2)
c1_mesh1_vert = compute_selected_vertices_cones(dir = pset1$dirs, complex =c1_mesh_1, rate_vals = comp1$Rate2[,2], len = pset1$len, threshold = ((pset1$directions_per_cone)/(dim(comp1$Rate2)[1])),
                                              cone_size = pset1$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

c1_mesh2_vert = compute_selected_vertices_cones(dir = pset1$dirs, complex =c1_mesh_2, rate_vals = comp1$Rate2[,2], len = pset1$len, threshold = ((pset1$directions_per_cone)/(dim(comp1$Rate2)[1])),
                                              cone_size = pset1$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
c1_mesh1_colors = rep('white', dim(c1_mesh1$vb)[2])
c1_mesh1_colors[c1_mesh1_vert] = 'red'

plot3d(c1_mesh1,col = c1_mesh1_colors)
rgl.viewpoint(userMatrix = rotation_matrix)

c1_mesh2_colors = rep('white', dim(c1_mesh2$vb)[2])
c1_mesh2_colors[c1_mesh2_vert] = 'red'

plot3d(c1_mesh2,col = c1_mesh2_colors)
rgl.viewpoint(userMatrix = rotation_matrix)


c1_mesh1_heat =  reconstruct_vertices_on_shape(dir = pset1$dirs,complex = c1_mesh_1,rate_vals = comp1$Rate2[,2],
                                            len = pset1$len,cuts = 1000,cone_size = pset1$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)
c1_mesh2_heat =  reconstruct_vertices_on_shape(dir = pset1$dirs,complex = c1_mesh_2,rate_vals = comp1$Rate2[,2],
                                            len = pset1$len,cuts = 1000,cone_size = pset1$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)

c1_mesh1_colors_heat=colfunc(1 + max(c1_mesh1_heat[,1]) - min(c1_mesh1_heat[,1]))[1 + c1_mesh1_heat[,1] - min(c1_mesh1_heat[,1])]
c1_mesh2_colors_heat=colfunc(1 + max(c1_mesh2_heat[,1]) - min(c1_mesh2_heat[,1]))[1 + c1_mesh2_heat[,1] - min(c1_mesh2_heat[,1])]

c1_mesh1_colors_heat=colfunc(1 + max(c1_mesh1_heat[,1]) )[1 + c1_mesh1_heat[,1] ]
c1_mesh2_colors_heat=colfunc(1 + max(c1_mesh2_heat[,1]) )[1 + c1_mesh2_heat[,1]]

mfrow3d(1,2)
plot3d(c1_mesh1,col = c1_mesh1_colors_heat, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(fruit1,fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(c1_mesh2,col = c1_mesh2_colors_heat, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(veg1,veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)


c1_mesh1_colors_list = get_heat_colors(dir_name = pset1$dir1,cuts = 300,pset = pset1,comp = comp1, colfunc = colfunc)
c1_mesh2_colors_list = get_heat_colors(dir_name = pset1$dir2,cuts = 300,pset = pset1,comp = comp1, colfunc = colfunc)



#### New Fruit Bug ####

ind = 2
c2_mesh1 = vcgImport(list.files(pset2$dir1,full.names = TRUE)[ind])
c2_mesh2 = vcgImport(list.files(pset2$dir2,full.names = TRUE)[ind])
c2_mesh_1 = process_off_file_v3(list.files(pset2$dir1,full.names = TRUE)[ind])
c2_mesh_2 = process_off_file_v3(list.files(pset2$dir2,full.names = TRUE)[ind])


mfrow3d(1,2)
c2_mesh1_vert = compute_selected_vertices_cones(dir = pset2$dirs, complex =c2_mesh_1, rate_vals = comp2$Rate2[,2], len = pset2$len, threshold = ((pset2$directions_per_cone)/(dim(comp2$Rate2)[1])),
                                                cone_size = pset2$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

c2_mesh2_vert = compute_selected_vertices_cones(dir = pset2$dirs, complex =c2_mesh_2, rate_vals = comp2$Rate2[,2], len = pset2$len, threshold = ((pset2$directions_per_cone)/(dim(comp2$Rate2)[1])),
                                                cone_size = pset2$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
c2_mesh1_colors = rep('white', dim(c2_mesh1$vb)[2])
c2_mesh1_colors[c2_mesh1_vert] = 'red'

plot3d(c2_mesh1,col = c2_mesh1_colors)
rgl.viewpoint(userMatrix = rotation_matrix)

c2_mesh2_colors = rep('white', dim(c2_mesh2$vb)[2])
c2_mesh2_colors[c2_mesh2_vert] = 'red'

plot3d(c2_mesh2,col = c2_mesh2_colors)
rgl.viewpoint(userMatrix = rotation_matrix)


c2_mesh1_heat =  reconstruct_vertices_on_shape(dir = pset2$dirs,complex = c2_mesh_1,rate_vals = comp2$Rate2[,2],
                                               len = pset2$len,cuts = 100,cone_size = pset2$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)
c2_mesh2_heat =  reconstruct_vertices_on_shape(dir = pset2$dirs,complex = c2_mesh_2,rate_vals = comp2$Rate2[,2],
                                               len = pset2$len,cuts = 100,cone_size = pset2$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)

c2_mesh1_colors_heat=colfunc(1 + max(c2_mesh1_heat[,1]) - min(c2_mesh1_heat[,1]))[1 + c2_mesh1_heat[,1] - min(c2_mesh1_heat[,1])]
c2_mesh2_colors_heat=colfunc(1 + max(c2_mesh2_heat[,1]) - min(c2_mesh2_heat[,1]))[1 + c2_mesh2_heat[,1] - min(c2_mesh2_heat[,1])]

mfrow3d(1,2)
plot3d(c2_mesh1,col = c2_mesh1_colors_heat, back="lines", specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(c2_mesh2,col = c2_mesh2_colors_heat, back="lines", specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)

c2_mesh1_colors_list = get_heat_colors(dir_name = pset2$dir1,cuts = 300,pset = pset2,comp = comp2, colfunc = colfunc)
c2_mesh2_colors_list = get_heat_colors(dir_name = pset2$dir2,cuts = 300,pset = pset2,comp = comp2, colfunc = colfunc)

#### New Veg Bug ####
ind = 2
c3_mesh1 = vcgImport(list.files(pset3$dir1,full.names = TRUE)[ind])
c3_mesh2 = vcgImport(list.files(pset3$dir2,full.names = TRUE)[ind])
c3_mesh_1 = process_off_file_v3(list.files(pset3$dir1,full.names = TRUE)[ind])
c3_mesh_2 = process_off_file_v3(list.files(pset3$dir2,full.names = TRUE)[ind])

#Also finding the landmarks

mfrow3d(1,2)
c3_mesh1_vert = compute_selected_vertices_cones(dir = pset3$dirs, complex =c3_mesh_1, rate_vals = comp3$Rate2[,2], len = pset3$len, threshold = ((pset3$directions_per_cone)/(dim(comp3$Rate2)[1])),
                                                cone_size = pset3$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

c3_mesh2_vert = compute_selected_vertices_cones(dir = pset3$dirs, complex =c3_mesh_2, rate_vals = comp3$Rate2[,2], len = pset3$len, threshold = ((pset3$directions_per_cone)/(dim(comp3$Rate2)[1])),
                                                cone_size = pset3$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
c3_mesh1_colors = rep('white', dim(c3_mesh1$vb)[2])
c3_mesh1_colors[c3_mesh1_vert] = 'red'

plot3d(c3_mesh1,col = c3_mesh1_colors)
rgl.viewpoint(userMatrix = rotation_matrix)

c3_mesh2_colors = rep('white', dim(c3_mesh2$vb)[2])
c3_mesh2_colors[c3_mesh2_vert] = 'red'

plot3d(c3_mesh2,col = c3_mesh2_colors)
rgl.viewpoint(userMatrix = rotation_matrix)


c3_mesh1_heat =  reconstruct_vertices_on_shape(dir = pset3$dirs,complex = c3_mesh_1,rate_vals = comp3$Rate2[,2],
                                               len = pset3$len,cuts = 100,cone_size = pset3$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)
c3_mesh2_heat =  reconstruct_vertices_on_shape(dir = pset3$dirs,complex = c3_mesh_2,rate_vals = comp3$Rate2[,2],
                                               len = pset3$len,cuts = 100,cone_size = pset3$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)

c3_mesh1_colors_heat=colfunc(1+max(c3_mesh1_heat[,1]) - min(c3_mesh1_heat[,1]))[1+c3_mesh1_heat[,1] - min(c3_mesh1_heat[,1])]
c3_mesh2_colors_heat=colfunc(1+max(c3_mesh2_heat[,1]) - min(c3_mesh2_heat[,1]))[1+c3_mesh2_heat[,1] - min(c3_mesh2_heat[,1])]

mfrow3d(1,2)
plot3d(c3_mesh1,col = c3_mesh1_colors_heat, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(fruit1,fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(c3_mesh2,col = c3_mesh2_colors_heat, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(veg1,veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

c3_mesh1_colors_list = get_heat_colors(dir_name = pset3$dir1,cuts = 300,pset = pset3,comp = comp3, colfunc = colfunc)
c3_mesh2_colors_list = get_heat_colors(dir_name = pset3$dir2,cuts = 300,pset = pset3,comp = comp3, colfunc = colfunc)

mfrow3d(3,10)
for (k in 1:length(list.files(pset3$dir1,full.names = TRUE))){
  mesh = vcgImport(list.files(pset3$dir1,full.names = TRUE)[k])
  plot3d(mesh, col = c3_mesh1_colors_list[[k]],back="lines", specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
  rgl.viewpoint(userMatrix = rotation_matrix)
}
for (k in 1:length(list.files(pset3$dir2,full.names = TRUE))){
  mesh = vcgImport(list.files(pset3$dir2,full.names = TRUE)[k])
  plot3d(mesh, col = c3_mesh2_colors_list[[k]],back="lines", specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
  rgl.viewpoint(userMatrix = rotation_matrix)
}

mesh = vcgImport('mushroom.off')
plot3d(mesh, col = 'white')
rgl.viewpoint(userMatrix = rotation_matrix)
save.image('new_world_monkey_comparison.Rdata')
