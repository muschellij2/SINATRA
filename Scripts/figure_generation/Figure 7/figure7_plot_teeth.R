set.seed(55)
library(sinatra)
library(FNN)
library(rgl)
library(Rvcg)
library(plyr)

#Parameters for the Analysis
pset1 = list(dir1 = '~/Documents/new_teeth_ect_aligned_by_species/Microcebus/',dir2 = '~/Documents/new_teeth_ect_aligned_by_species/Mirza/',
             num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
                dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset2 =  list(dir1 = '~/Documents/new_teeth_ect_aligned_by_species/Microcebus/',dir2 = '~/Documents/new_teeth_ect_aligned_by_species/Tarsius/',
              num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
              dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset3 = list(dir1 = '~/Documents/new_teeth_ect_aligned_by_species/Microcebus/',dir2 = '~/Documents/new_teeth_ect_aligned_by_species/Saimiri/',
             num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset4 = list(dir1 = '~/Documents/new_teeth_ect_aligned_by_species/Tarsius/',dir2 = '~/Documents/new_teeth_ect_aligned_by_species/Saimiri/',
             num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset5 = list(dir1 = '~/Documents/new_teeth_ect_aligned_by_species/Tarsius/',dir2 = '~/Documents/new_teeth_ect_aligned_by_species/Mirza/',
             num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset6 = list(dir1 = '~/Documents/new_teeth_ect_aligned_by_species/Saimiri/',dir2 = '~/Documents/new_teeth_ect_aligned_by_species/Mirza/',
             num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))


cap_radius = 0.25
num_cones = 35
directions_per_cone = 5
len = 75
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  cap_radius, directions_per_cone = directions_per_cone)
rotation_matrix=matrix(c(0.8065218,0.5911149,0.01028626,0,0.5583186,-0.7558216,0.34207344,0,0.1944301,-0.2816325,-0.93961692,0,0,0,0,1),ncol=4,byrow=TRUE)
ball = TRUE
ball_radius = 0.5
ec_type = 'ECT'
#### Start Comparison ####



load('~/Documents/new_teeth_ect_aligned_by_species/rate_analysis.Rdata')

#### Comp 2 ####
indices = c(261,295,833,224,939,278,252,593,177,357,4845,4215,4619,162,244,4618,4225,31,4705,4524,219,4301,4898,262,4860,69,64,4478,4452,4224,102,179,4176)

color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)

col_pal2=c(color1,color1,color2,color2,color3)
colfunc2 <- colorRampPalette(col_pal2)

ind = 11

vertex_index = indices[ind]
mesh2 = vcgImport(list.files(pset2$dir2,full.names = TRUE)[ind])
mesh_2 = process_off_file_v3(list.files(pset2$dir2,full.names = TRUE)[ind])

#Also finding the landmarks


mesh2_heat =  reconstruct_vertices_on_shape(dir = pset4$dirs,complex = mesh_2,rate_vals = comp4$Rate2,
                                            len = pset4$len,cuts = 1000,cone_size = pset2$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)
mesh4_2_heat =  reconstruct_vertices_on_shape(dir = pset5$dirs,complex = mesh_2,rate_vals = comp5$Rate2,
                                              len = pset5$len,cuts = 1000,cone_size = pset5$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)
mesh5_2_heat =  reconstruct_vertices_on_shape(dir = pset2$dirs,complex = mesh_2,rate_vals = comp2$Rate2,
                                              len = pset2$len,cuts = 1000,cone_size = pset5$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)

mesh2_colors_heat=colfunc(1 + max(mesh2_heat[,1]) - min(mesh2_heat[,1]))[1 + mesh2_heat[,1] - min(mesh2_heat[,1])]
mesh4_2_colors_heat=colfunc(max(mesh4_2_heat[,1]+1) - min(mesh4_2_heat[,1]))[mesh4_2_heat[,1] - min(mesh4_2_heat[,1])+1]
mesh5_2_colors_heat=colfunc(1 + max(mesh5_2_heat[,1]) - min(mesh5_2_heat[,1]))[1 + mesh5_2_heat[,1] - min(mesh5_2_heat[,1])]

closest_vertices = knnx.index(query = matrix(t(mesh2$vb[-4,vertex_index]),ncol = 3), data = t(mesh2$vb[-4,]), k = 100)
open3d()
rotation_matrix=matrix(c(0.919328629970551,0.389629483222961,0.0549933537840843,0,0.390242397785187,-0.920712232589722,-0.000437354668974876,0,0.0504627674818039,0.0218627471476793,-0.998486816883087,0,0,0,0,1), ncol = 4, byrow = TRUE)

mfrow3d(1,3)
colors = rep('white', dim(mesh2$vb)[2])
colors[closest_vertices] = 'red'
plot3d(mesh2,col = colors,  specular="black", axes = FALSE,xlab = '',ylab = '',zlab='')
#rgl.points(matrix(t(mesh2$vb[-4,vertex_index]),ncol = 3), col = 'blue', size = 15)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(mesh2,col = mesh2_colors_heat,  specular="black", axes = FALSE, xlab = '',ylab = '',zlab='')
rgl.viewpoint(userMatrix = rotation_matrix)


plot3d(mesh2,col = mesh4_2_colors_heat,  specular="black", axes = FALSE, xlab = '',ylab = '',zlab='')
rgl.viewpoint(userMatrix = rotation_matrix)


plot3d(mesh2,col = mesh5_2_colors_heat,specular="black", axes = FALSE,  xlab = '',ylab = '',zlab='')
rgl.viewpoint(userMatrix = rotation_matrix)

