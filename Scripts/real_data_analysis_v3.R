setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
setwd('~/Documents/')
library(R.utils)
sourceDirectory('~/Documents/SINATRA/SINATRA_Pipeline_Branch/')
sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')

#library(svd)
set.seed(15)

#Specifying Directories and the Associated Files

#old_angry_dir = 'Data/HDM/by_diet_scaled_aligned/Old_Frugivore/'
#old_smiling_dir='Data/HDM/by_diet_scaled_aligned/Old_Follivore/'
#old_bug_dir='Data/HDM//by_diet_scaled_aligned/Old_Insectivore/'
#new_angry_dir = 'Data/HDM/by_diet_scaled_aligned/New_Frugivore/'
#new_smiling_dir='Data/HDM/by_diet_scaled_aligned/New_Follivore/'
#new_bug_dir='Data/HDM//by_diet_scaled_aligned/New_Insectivore/'

smiling_dir = '~/Documents/face_meshes/new_aligned_meshes/smiling/'
angry_dir = '~/Documents/face_meshes/new_aligned_meshes/angry/'

angry_files =list.files(path=angry_dir,full.names = TRUE)
smiling_files =list.files(path=smiling_dir,full.names = TRUE)
#Parameters for the Analysis
cap_radius = 0.15
num_cones = 25
directions_per_cone = 5
len = 100
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  cap_radius, directions_per_cone = directions_per_cone)
ball = FALSE
ball_radius = 0.5
ec_type = 'ECT'
#### Start Comparison #### 

smile_angry=real_data_summary(dir1=smiling_dir,dir2 = angry_dir,direction=dirs,class1='Smiling', class2='Angry',
                                            radius=0,accuracy=FALSE,len = len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)

color1='blue'
color2 = 'blue'
color3 = 'lightblue'
color4='palegreen'
color5='orangered'
col_pal=c(color1,color1,color2,color2,color2, color3,color3,color3, color4,color4,color4, color5)
#col_pal =c(color1,color2,color3)
colfunc <- colorRampPalette(col_pal)

#### New angry smiling ####
ind = 25
rotation_matrix=matrix(c(0.99972576,0.02127766,0.00978078,0,0.01589701,-0.92330807,0.38373107,0,0.01719557,-0.38347024,-0.92339307,0,0,0,0,1),ncol=4,byrow=TRUE)
angry1 = vcgImport(angry_files[ind])
smiling1 = vcgImport(smiling_files[ind])
angry_1 = process_off_file_v3(angry_files[ind])
smiling_1 = process_off_file_v3(smiling_files[ind])

#Also finding the landmarks

mfrow3d(1,2)
angry1_vert = compute_selected_vertices_cones(dir = dirs, complex =angry_1, rate_vals = smile_angry$Rate2[,2], len = len, threshold = ((directions_per_cone)/(dim(smile_angry$Rate2)[1])),
                                              cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

smiling1_vert = compute_selected_vertices_cones(dir = dirs, complex = smiling_1, rate_vals = smile_angry$Rate2[,2], len = len, threshold = ((directions_per_cone)/(dim(smile_angry$Rate2)[1])),
                                            cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
angry_colors = rep('white', dim(angry1$vb)[2])
angry_colors[angry1_vert] = 'red'

plot3d(angry1,col = angry_colors)
#plot_selected_landmarks(angry1,angry1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

smiling_colors = rep('white', dim(smiling1$vb)[2])
smiling_colors[smiling1_vert] = 'blue'
plot3d(smiling1,col = smiling_colors)
#plot_selected_landmarks(smiling1,smiling1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

angry_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = angry_1,rate_vals = smile_angry$Rate2[,2],
                                            len = len,cuts = 200,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)
smiling_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = smiling_1,rate_vals = smile_angry$Rate2[,2],
                                          len = len,cuts = 200,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

angry_colors=colfunc(max(angry_heat[,1]))[angry_heat[,1]]
smiling_colors=colfunc(max(smiling_heat[,1]))[smiling_heat[,1]]

mfrow3d(1,2)
plot3d(angry1,col = angry_colors, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(angry1,angry1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(smiling1,col = smiling_colors, back="lines", specular="black", axes = FALSE)
#plot_selected_landmarks(smiling1,smiling1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)


#### New angry Bug ####
