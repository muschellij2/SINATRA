setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
library(R.utils)
sourceDirectory('~/Documents/SINATRA/SINATRA_Pipeline_Branch/')
sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')

#library(svd)
#set.seed(15)
#set.seed(35)
set.seed(45)

#Specifying Directories and the Associated Files
old_bug_dir='Data/new_aligned_shapesv3/V13_v2/gp1'
new_bug_dir='Data/new_aligned_shapesv3/V13_v2/gp2'
old_bug_files=list.files(path=old_bug_dir,full.names = TRUE)

new_bug_files=list.files(path=new_bug_dir,full.names = TRUE)


files = c(old_fruit_files, old_veg_files, old_bug_files, new_fruit_files, new_veg_files, new_bug_files)
#Parameters for the Analysis
cap_radius = 0.25
num_cones = 35
directions_per_cone = 5
len = 75
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  cap_radius, directions_per_cone = directions_per_cone)
ball = TRUE
ball_radius = 0.5
ec_type = 'ECT'
#### Start Comparison #### 
#35, 75 length 0.25 cap radius
bug_bug_old=real_data_summary(dir1=new_bug_dir,dir2 = old_bug_dir,direction=dirs,class1='Bug', class2='Bug',
                                          radius=0,accuracy=FALSE,len = len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)

rotation_matrix=matrix(c(0.99972576,0.02127766,0.00978078,0,0.01589701,-0.92330807,0.38373107,0,0.01719557,-0.38347024,-0.92339307,0,0,0,0,1),ncol=4,byrow=TRUE)

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
new_bug1_vert = compute_selected_vertices_cones(dir = dirs, complex =new_bug_1, rate_vals = bug_bug_old$Rate2[,2], len = len, threshold = (directions_per_cone/(1*dim(bug_bug_old$Rate2)[1])),
                                                cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball)

old_bug1_vert = compute_selected_vertices_cones(dir = dirs, complex = old_bug_1, rate_vals = bug_bug_old$Rate2[,2], len = len, threshold = (directions_per_cone/(1*dim(bug_bug_old$Rate2)[1])),
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





top_dirs = c(1,3)
find_dirs = matrix(dirs[1,],ncol = 3)

d =  reconstruct_vertices_on_shape(dir = find_dirs,complex = new_bug_1,rate_vals = bug_bug_old$Rate2[,2],
                                   len = len,cuts = 100,cone_size = 1,ball_radius = ball_radius, ball = ball)

test_colors=colfunc(max(d[,1])- min(d[,1]))[(d[,1] -min(d[,1]))]
plot3d(new_bug1,col = test_colors, specular="black")

