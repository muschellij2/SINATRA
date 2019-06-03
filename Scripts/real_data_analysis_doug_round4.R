setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
library(R.utils)
sourceDirectory('~/Documents/SINATRA/SINATRA_Pipeline_Branch/')
sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')
set.seed(55)
#load('chris_doug_analysis.Rdata')
#Specifying Directories and the Associated Files

#old_fruit_dir = 'Data/HDM/by_diet_scaled_aligned/Old_Frugivore/'
#old_veg_dir='Data/HDM/by_diet_scaled_aligned/Old_Follivore/'
#old_bug_dir='Data/HDM//by_diet_scaled_aligned/Old_Insectivore/'
#new_fruit_dir = 'Data/HDM/by_diet_scaled_aligned/New_Frugivore/'
#new_veg_dir='Data/HDM/by_diet_scaled_aligned/New_Follivore/'
#new_bug_dir='Data/HDM//by_diet_scaled_aligned/New_Insectivore/'

#Parameters for the Analysis
pset1 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Microcebus/',dir2 = '~/Documents/doug_new_teeth_by_species/Mirza/',
             num_cones = 35, cap_radius = 0.15, len = 100, directions_per_cone = 5,
                dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))
pset2 =  list(dir1 = '~/Documents/doug_new_teeth_by_species/Microcebus/',dir2 = '~/Documents/doug_new_teeth_by_species/Tarsius/',
              num_cones = 35, cap_radius = 0.15, len = 100, directions_per_cone = 5,
              dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))
pset3 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Microcebus/',dir2 = '~/Documents/doug_new_teeth_by_species/Saimiri/',
             num_cones = 35, cap_radius = 0.15, len = 100, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))
pset4 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Tarsius/',dir2 = '~/Documents/doug_new_teeth_by_species/Saimiri/',
             num_cones = 35, cap_radius = 0.15, len = 100, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))
pset5 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Tarsius/',dir2 = '~/Documents/doug_new_teeth_by_species/Mirza/',
             num_cones = 35, cap_radius = 0.15, len = 100, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))
pset6 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Saimiri/',dir2 = '~/Documents/doug_new_teeth_by_species/Mirza/',
             num_cones = 35, cap_radius = 0.15, len = 100, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))

cap_radius = 0.15
num_cones = 35
directions_per_cone = 5
len = 100
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  cap_radius, directions_per_cone = directions_per_cone)
ball = TRUE
ball_radius = 0.5
ec_type = 'ECT'
#### Start Comparison ####

comp1 = real_data_summary(dir1=pset1$dir1,dir2 = pset1$dir2,direction=pset1$dirs,class1=pset1$dir1, class2=pset1$dir1,
                                    radius=0,accuracy=FALSE,len = pset1$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
comp2 = real_data_summary(dir1=pset2$dir1,dir2 = pset2$dir2,direction=pset2$dirs,class1=pset2$dir1, class2=pset2$dir1,
                          radius=0,accuracy=FALSE,len = pset2$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
comp3 = real_data_summary(dir1=pset3$dir1,dir2 = pset3$dir2,direction=pset3$dirs,class1=pset3$dir1, class2=pset3$dir1,
                          radius=0,accuracy=FALSE,len = pset3$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
comp4 = real_data_summary(dir1=pset4$dir1,dir2 = pset4$dir2,direction=pset4$dirs,class1=pset4$dir1, class2=pset4$dir1,
                          radius=0,accuracy=FALSE,len = pset4$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
comp5 = real_data_summary(dir1=pset5$dir1,dir2 = pset5$dir2,direction=pset5$dirs,class1=pset5$dir1, class2=pset5$dir1,
                          radius=0,accuracy=FALSE,len = pset5$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
comp6 = real_data_summary(dir1=pset6$dir1,dir2 = pset6$dir2,direction=pset6$dirs,class1=pset6$dir1, class2=pset6$dir1,
                          radius=0,accuracy=FALSE,len = pset6$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
#comp7 = real_data_summary(dir1=pset7$dir1,dir2 = pset7$dir2,direction=pset7$dirs,class1=pset7$dir1, class2=pset7$dir1,
#                          radius=0,accuracy=FALSE,len = pset7$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)





color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)

col_pal2=c(color1,color1,color2,color2,color3)
colfunc2 <- colorRampPalette(col_pal2)


#### Comp 1 ####
ind = 2
rotation_matrix=matrix(c(0.8065218,0.5911149,0.01028626,0,0.5583186,-0.7558216,0.34207344,0,0.1944301,-0.2816325,-0.93961692,0,0,0,0,1),ncol=4,byrow=TRUE)
c1_mesh1 = vcgImport(list.files(pset1$dir1,full.names = TRUE)[ind])
c1_mesh2 = vcgImport(list.files(pset1$dir2,full.names = TRUE)[ind])
c1_mesh_1 = process_off_file_v3(list.files(pset1$dir1,full.names = TRUE)[ind])
c1_mesh_2 = process_off_file_v3(list.files(pset1$dir2,full.names = TRUE)[ind])

#Also finding the landmarks

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
                                            len = pset1$len,cuts = 300,cone_size = pset1$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)
c1_mesh2_heat =  reconstruct_vertices_on_shape(dir = pset1$dirs,complex = c1_mesh_2,rate_vals = comp1$Rate2[,2],
                                            len = pset1$len,cuts = 300,cone_size = pset1$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)

c1_mesh1_colors_heat=colfunc(1 + max(c1_mesh1_heat[,1]) - min(c1_mesh1_heat[,1]))[1 + c1_mesh1_heat[,1] - min(c1_mesh1_heat[,1])]
c1_mesh2_colors_heat=colfunc(1 + max(c1_mesh2_heat[,1]) - min(c1_mesh2_heat[,1]))[1 + c1_mesh2_heat[,1] - min(c1_mesh2_heat[,1])]

mfrow3d(1,2)
plot3d(c1_mesh1,col = c1_mesh1_colors_heat, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(fruit1,fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(c1_mesh2,col = c1_mesh2_colors_heat, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(veg1,veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)


c1_mesh1_colors_list = get_heat_colors(dir_name = pset1$dir1,cuts = 300,pset = pset1,comp = comp1, colfunc = colfunc)
c1_mesh2_colors_list = get_heat_colors(dir_name = pset1$dir2,cuts = 300,pset = pset1,comp = comp1, colfunc = colfunc)
#### Comp 2 ####

ind = 2
c2_mesh1 = vcgImport(list.files(pset2$dir1,full.names = TRUE)[ind])
c2_mesh2 = vcgImport(list.files(pset2$dir2,full.names = TRUE)[ind])
c2_mesh_1 = process_off_file_v3(list.files(pset2$dir1,full.names = TRUE)[ind])
c2_mesh_2 = process_off_file_v3(list.files(pset2$dir2,full.names = TRUE)[ind])

#Also finding the landmarks

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

#### Comp 3 ####
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

mfrow3d(2,4)
for (k in 1:length(list.files(pset3$dir2,full.names = TRUE))){
  mesh = vcgImport(list.files(pset3$dir2,full.names = TRUE)[k])
  plot3d(mesh, col = c3_mesh2_colors_list[[k]],back="lines", specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
  rgl.viewpoint(userMatrix = rotation_matrix)
}




#### Comp 4 ####

ind = 2
c4_mesh1 = vcgImport(list.files(pset4$dir1,full.names = TRUE)[ind])
c4_mesh2 = vcgImport(list.files(pset4$dir2,full.names = TRUE)[ind])
c4_mesh_1 = process_off_file_v3(list.files(pset4$dir1,full.names = TRUE)[ind])
c4_mesh_2 = process_off_file_v3(list.files(pset4$dir2,full.names = TRUE)[ind])

#Also finding the landmarks

mfrow3d(1,2)
c4_mesh1_vert = compute_selected_vertices_cones(dir = pset4$dirs, complex =c4_mesh_1, rate_vals = comp4$Rate2[,2], len = pset4$len, threshold = ((pset4$directions_per_cone)/(dim(comp4$Rate2)[1])),
                                                cone_size = pset4$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

c4_mesh2_vert = compute_selected_vertices_cones(dir = pset4$dirs, complex =c4_mesh_2, rate_vals = comp4$Rate2[,2], len = pset4$len, threshold = ((pset4$directions_per_cone)/(dim(comp4$Rate2)[1])),
                                                cone_size = pset4$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
c4_mesh1_colors = rep('white', dim(c4_mesh1$vb)[2])
c4_mesh1_colors[c4_mesh1_vert] = 'red'
plot3d(c4_mesh1,col = c4_mesh1_colors)
rgl.viewpoint(userMatrix = rotation_matrix)

c4_mesh2_colors = rep('white', dim(c4_mesh2$vb)[2])
c4_mesh2_colors[c4_mesh2_vert] = 'red'

plot3d(c4_mesh2,col = c4_mesh2_colors)
rgl.viewpoint(userMatrix = rotation_matrix)


c4_mesh1_heat =  reconstruct_vertices_on_shape(dir = pset4$dirs,complex = c4_mesh_1,rate_vals = comp4$Rate2[,2],
                                               len = pset4$len,cuts = 100,cone_size = pset4$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)
c4_mesh2_heat =  reconstruct_vertices_on_shape(dir = pset4$dirs,complex = c4_mesh_2,rate_vals = comp4$Rate2[,2],
                                               len = pset4$len,cuts = 100,cone_size = pset4$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)

c4_mesh1_colors_heat=colfunc(max(c4_mesh1_heat[,1]+1) - min(c4_mesh1_heat[,1]))[c4_mesh1_heat[,1] - min(c4_mesh1_heat[,1])+1]
c4_mesh2_colors_heat=colfunc(max(c4_mesh2_heat[,1]+1) - min(c4_mesh2_heat[,1]))[c4_mesh2_heat[,1] - min(c4_mesh2_heat[,1])+1]

mfrow3d(1,2)
plot3d(c4_mesh1,col = c4_mesh1_colors_heat, back="lines", specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(c4_mesh2,col = c4_mesh2_colors_heat, back="lines", specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)

c4_mesh1_colors_list = get_heat_colors(dir_name = pset4$dir1,cuts = 300,pset = pset4,comp = comp4, colfunc = colfunc)
c4_mesh2_colors_list = get_heat_colors(dir_name = pset4$dir2,cuts = 300,pset = pset4,comp = comp4, colfunc = colfunc)


#### Comp 5 ####
ind = 2
c5_mesh1 = vcgImport(list.files(pset5$dir1,full.names = TRUE)[ind])
c5_mesh2 = vcgImport(list.files(pset5$dir2,full.names = TRUE)[ind])
c5_mesh_1 = process_off_file_v3(list.files(pset5$dir1,full.names = TRUE)[ind])
c5_mesh_2 = process_off_file_v3(list.files(pset5$dir2,full.names = TRUE)[ind])

#Also finding the landmarks

mfrow3d(1,2)
c5_mesh1_vert = compute_selected_vertices_cones(dir = pset5$dirs, complex =c5_mesh_1, rate_vals = comp5$Rate2[,2], len = pset5$len, threshold = ((pset5$directions_per_cone)/(dim(comp5$Rate2)[1])),
                                                cone_size = pset5$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

c5_mesh2_vert = compute_selected_vertices_cones(dir = pset5$dirs, complex =c5_mesh_2, rate_vals = comp5$Rate2[,2], len = pset5$len, threshold = ((pset5$directions_per_cone)/(dim(comp5$Rate2)[1])),
                                                cone_size = pset5$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
c5_mesh1_colors = rep('white', dim(c5_mesh1$vb)[2])
c5_mesh1_colors[c5_mesh1_vert] = 'red'

plot3d(c5_mesh1,col = c5_mesh1_colors)
rgl.viewpoint(userMatrix = rotation_matrix)

c5_mesh2_colors = rep('white', dim(c5_mesh2$vb)[2])
c5_mesh2_colors[c5_mesh2_vert] = 'red'

plot3d(c5_mesh2,col = c5_mesh2_colors)
rgl.viewpoint(userMatrix = rotation_matrix)


c5_mesh1_heat =  reconstruct_vertices_on_shape(dir = pset5$dirs,complex = c5_mesh_1,rate_vals = comp5$Rate2[,2],
                                               len = pset5$len,cuts = 100,cone_size = pset5$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)
c5_mesh2_heat =  reconstruct_vertices_on_shape(dir = pset5$dirs,complex = c5_mesh_2,rate_vals = comp5$Rate2[,2],
                                               len = pset5$len,cuts = 100,cone_size = pset5$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)

c5_mesh1_colors_heat=colfunc(1 + max(c5_mesh1_heat[,1]) - min(c5_mesh1_heat[,1]))[1 + c5_mesh1_heat[,1] - min(c5_mesh1_heat[,1])]
c5_mesh2_colors_heat=colfunc(1 + max(c5_mesh2_heat[,1]) - min(c5_mesh2_heat[,1]))[1 + c5_mesh2_heat[,1] - min(c5_mesh2_heat[,1])]

mfrow3d(1,2)
plot3d(c5_mesh1,col = c5_mesh1_colors_heat, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(fruit1,fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(c5_mesh2,col = c5_mesh2_colors_heat, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(veg1,veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

c5_mesh1_colors_list = get_heat_colors(dir_name = pset5$dir1,cuts = 300,pset = pset5,comp = comp5, colfunc = colfunc)
c5_mesh2_colors_list = get_heat_colors(dir_name = pset5$dir2,cuts = 300,pset = pset5,comp = comp5, colfunc = colfunc)


#### Comp 6 ####

ind = 2
c6_mesh1 = vcgImport(list.files(pset6$dir1,full.names = TRUE)[ind])
c6_mesh2 = vcgImport(list.files(pset6$dir2,full.names = TRUE)[ind])
c6_mesh_1 = process_off_file_v3(list.files(pset6$dir1,full.names = TRUE)[ind])
c6_mesh_2 = process_off_file_v3(list.files(pset6$dir2,full.names = TRUE)[ind])

#Also finding the landmarks

mfrow3d(1,2)
c6_mesh1_vert = compute_selected_vertices_cones(dir = pset6$dirs, complex =c6_mesh_1, rate_vals = comp6$Rate2[,2], len = pset6$len, threshold = ((pset6$directions_per_cone)/(dim(comp6$Rate2)[1])),
                                                cone_size = pset6$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

c6_mesh2_vert = compute_selected_vertices_cones(dir = pset6$dirs, complex =c6_mesh_2, rate_vals = comp6$Rate2[,2], len = pset6$len, threshold = ((pset6$directions_per_cone)/(dim(comp6$Rate2)[1])),
                                                cone_size = pset6$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
c6_mesh1_colors = rep('white', dim(c6_mesh1$vb)[2])
c6_mesh1_colors[c6_mesh1_vert] = 'red'

plot3d(c6_mesh1,col = c6_mesh1_colors)
rgl.viewpoint(userMatrix = rotation_matrix)

c6_mesh2_colors = rep('white', dim(c6_mesh2$vb)[2])
c6_mesh2_colors[c6_mesh2_vert] = 'red'

plot3d(c6_mesh2,col = c6_mesh2_colors)
rgl.viewpoint(userMatrix = rotation_matrix)


c6_mesh1_heat =  reconstruct_vertices_on_shape(dir = pset6$dirs,complex = c6_mesh_1,rate_vals = comp6$Rate2[,2],
                                               len = pset6$len,cuts = 100,cone_size = pset6$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)
c6_mesh2_heat =  reconstruct_vertices_on_shape(dir = pset6$dirs,complex = c6_mesh_2,rate_vals = comp6$Rate2[,2],
                                               len = pset6$len,cuts = 100,cone_size = pset6$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)

c6_mesh1_colors_heat=colfunc(1 + max(c6_mesh1_heat[,1]) - min(c6_mesh1_heat[,1]))[1 + c6_mesh1_heat[,1] - min(c6_mesh1_heat[,1])]
c6_mesh2_colors_heat=colfunc(1 + max(c6_mesh2_heat[,1]) - min(c6_mesh2_heat[,1]))[1 + c6_mesh2_heat[,1] - min(c6_mesh2_heat[,1])]

mfrow3d(1,2)
plot3d(c6_mesh1,col = c6_mesh1_colors_heat, back="lines", specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(c6_mesh2,col = c6_mesh2_colors_heat, back="lines", specular="black", axes = FALSE)
rgl.viewpoint(userMatrix = rotation_matrix)

c6_mesh1_colors_list = get_heat_colors(dir_name = pset6$dir1,cuts = 300,pset = pset6,comp = comp6, colfunc = colfunc)
c6_mesh2_colors_list = get_heat_colors(dir_name = pset6$dir2,cuts = 300,pset = pset6,comp = comp6, colfunc = colfunc)


#### Comp 7 ####

ind = 2
c7_mesh1 = vcgImport(list.files(pset7$dir1,full.names = TRUE)[ind])
c7_mesh2 = vcgImport(list.files(pset7$dir2,full.names = TRUE)[ind])
c7_mesh_1 = process_off_file_v3(list.files(pset7$dir1,full.names = TRUE)[ind])
c7_mesh_2 = process_off_file_v3(list.files(pset7$dir2,full.names = TRUE)[ind])

#Also finding the landmarks

mfrow3d(1,2)
c7_mesh1_vert = compute_selected_vertices_cones(dir = pset7$dirs, complex =c7_mesh_1, rate_vals = comp7$Rate2[,2], len = pset7$len, threshold = ((pset7$directions_per_cone)/(dim(comp7$Rate2)[1])),
                                                cone_size = pset7$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

c7_mesh2_vert = compute_selected_vertices_cones(dir = pset7$dirs, complex =c7_mesh_2, rate_vals = comp7$Rate2[,2], len = pset7$len, threshold = ((pset7$directions_per_cone)/(dim(comp7$Rate2)[1])),
                                                cone_size = pset7$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
c7_mesh1_colors = rep('white', dim(c7_mesh1$vb)[2])
c7_mesh1_colors[c7_mesh1_vert] = 'red'

plot3d(c7_mesh1,col = c7_mesh1_colors)
rgl.viewpoint(userMatrix = rotation_matrix)

c7_mesh2_colors = rep('white', dim(c7_mesh2$vb)[2])
c7_mesh2_colors[c7_mesh2_vert] = 'red'

plot3d(c7_mesh2,col = c7_mesh2_colors)
rgl.viewpoint(userMatrix = rotation_matrix)


c7_mesh1_heat =  reconstruct_vertices_on_shape(dir = pset7$dirs,complex = c7_mesh_1,rate_vals = comp7$Rate2[,2],
                                               len = pset7$len,cuts = 100,cone_size = pset7$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)
c7_mesh2_heat =  reconstruct_vertices_on_shape(dir = pset7$dirs,complex = c7_mesh_2,rate_vals = comp7$Rate2[,2],
                                               len = pset7$len,cuts = 100,cone_size = pset7$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)

c7_mesh1_colors_heat=colfunc(1 + max(c7_mesh1_heat[,1]) - min(c7_mesh1_heat[,1]))[1 + c7_mesh1_heat[,1] - min(c7_mesh1_heat[,1])]
c7_mesh2_colors_heat=colfunc(1 + max(c7_mesh2_heat[,1]) - min(c7_mesh2_heat[,1]))[1 + c7_mesh2_heat[,1] - min(c7_mesh2_heat[,1])]

mfrow3d(1,2)
plot3d(c7_mesh1,col = c7_mesh1_colors_heat, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(fruit1,fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(c7_mesh2,col = c7_mesh2_colors_heat, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(veg1,veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

c7_mesh1_colors_list = get_heat_colors(dir_name = pset7$dir1,cuts = 300,pset = pset7,comp = comp7, colfunc = colfunc)
c7_mesh2_colors_list = get_heat_colors(dir_name = pset7$dir2,cuts = 300,pset = pset7,comp = comp7, colfunc = colfunc)

####

save.image('chris_henry_analysis.Rdata')
