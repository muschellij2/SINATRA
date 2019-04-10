setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
library(R.utils)
sourceDirectory('~/Documents/SINATRA/SINATRA_Pipeline_Branch/')
sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')
set.seed(55)
load('Figure1.Rdata')

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
cap_radius = 0.15
num_cones = 35
directions_per_cone = 5
len = 100 
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  cap_radius, directions_per_cone = directions_per_cone)
ball = TRUE
ball_radius = 0.5
ec_type = 'ECT'
#### Start Comparison #### 
load('Teeth_v2.Rdata')
color1='blue'
color2='lightgreen'
color3='orangered'
col_pal=c(color1,color2,color2,color2,color2,color3)
#col_pal =c(color1,color2,color3)
colfunc <- colorRampPalette(col_pal)

#### New Fruit Veg ####
ind = 13
rotation_matrix=matrix(c(0.8065218,0.5911149,0.01028626,0,0.5583186,-0.7558216,0.34207344,0,0.1944301,-0.2816325,-0.93961692,0,0,0,0,1),ncol=4,byrow=TRUE)
veg1 = vcgImport(new_veg_files[ind])
fruit1 = vcgImport(new_fruit_files[ind])
fruit_1 = process_off_file_v3(new_fruit_files[ind])
veg_1 = process_off_file_v3(new_veg_files[ind])

#Also finding the landmarks

#1a
mfrow3d(1,2)
fruit_colors = rep('white', dim(fruit1$vb)[2])

plot3d(fruit1,col = fruit_colors, back="lines", specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
#plot_selected_landmarks(fruit1,fruit1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

veg_colors = rep('white', dim(veg1$vb)[2])
plot3d(veg1,col = veg_colors, back="lines", specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
#plot_selected_landmarks(veg1,veg1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

#1b
arrow_type = 'rotation'
s = 1/10
width = 1/9
thickness = 1/20
plot3d(fruit1,col = fruit_colors,  specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
for (k in 101:105){
  arrow3d(-0.5*dirs[101,],0.5*dirs[k,],s = s, type = arrow_type,width= width,thickness = thickness,col = 'red')
}
for (k in 91:95){
  arrow3d(-0.5*dirs[90,],0.5*dirs[k,],s = s, type = arrow_type,width= width,thickness = thickness, col = 'blue')
}
for (k in 56:60){
  arrow3d(-0.5*dirs[56,],0.65*dirs[k,],s = s, type = arrow_type,width= width,thickness = thickness, col = 'green')
}
#arrow3d(0.5*-dirs[1,],0.5*dirs[1,],s = 1/9, type = 'rotation',width= 1/9,thickness = 1/20)
#plot_selected_landmarks(fruit1,fruit1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

#1d
fruit_colors_heat=colfunc(max(fruit_heat[,1]) - min(fruit_heat[,1]))[fruit_heat[,1] - min(fruit_heat[,1])]
veg_colors_heat=colfunc(max(veg_heat[,1]) - min(veg_heat[,1]))[veg_heat[,1] - min(veg_heat[,1])]

mfrow3d(1,2)
plot3d(fruit1,col = fruit_colors_heat, back="lines", specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
#plot_selected_landmarks(fruit1,fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(veg1,col = veg_colors_heat, back="lines", specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
#plot_selected_landmarks(veg1,veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

save.image('Figure1.Rdata')
