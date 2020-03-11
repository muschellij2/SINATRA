setwd('/Users/brucewang/')
library(R.utils)
sourceDirectory('~/Documents/SINATRA/SINATRA_Pipeline_Branch/')
sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')

#library(svd)
#set.seed(15)
#set.seed(35)
set.seed(15)
load('Figure2.Rdata')
#Specifying Directories and the Associated Files

#old_fruit_dir = 'Data/HDM/by_diet_scaled_aligned/Old_Frugivore/'
#old_veg_dir='Data/HDM/by_diet_scaled_aligned/Old_Follivore/'
#old_bug_dir='Data/HDM//by_diet_scaled_aligned/Old_Insectivore/'
#new_fruit_dir = 'Data/HDM/by_diet_scaled_aligned/New_Frugivore/'
#new_veg_dir='Data/HDM/by_diet_scaled_aligned/New_Follivore/'
#new_bug_dir='Data/HDM//by_diet_scaled_aligned/New_Insectivore/'


old_veg_dir='Data/new_aligned_shapesv3/V13_v2/gp1'
new_veg_dir='Data/new_aligned_shapesv3/V13_v2/gp2'

old_veg_files=list.files(path=old_veg_dir,full.names = TRUE)

new_veg_files=list.files(path=new_veg_dir,full.names = TRUE)

veg_pset = list(num_cones = 35, cap_radius = 0.15, len = 100, directions_per_cone = 5,
                dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.15, directions_per_cone = 5))
#Parameters for the Analysis
ball = TRUE
ball_radius = 0.5
ec_type = 'ECT'
#### Start Comparison #### 

# 35, 0.25 cap radius

veg_veg_old=real_data_summary(dir1=new_veg_dir,dir2 = old_veg_dir,direction=veg_pset$dirs,class1='Vegetable', class2='Vegetable',
                                          radius=0,accuracy=FALSE,len = veg_pset$len, ball = ball, ball_radius = ball_radius, ec_type = ec_type)
#35, 75 length 0.25 cap radius
#0.15 , 75 length cap radius

color1='blue'
color2='lightgreen'
color3='orangered'
color3 = 'red'
col_pal=c(color1,color2,color2,color2,color3)
#col_pal =c(color1,color2,color3)
colfunc <- colorRampPalette(col_pal)

rotation_matrix=matrix(c(0.8065218,0.5911149,0.01028626,0,0.5583186,-0.7558216,0.34207344,0,0.1944301,-0.2816325,-0.93961692,0,0,0,0,1),ncol=4,byrow=TRUE)

#### Old Veg New Veg ####
#Figure 2c
class_1_probs = read.csv('Data/new_aligned_shapesv3/V13_v2/gp1_spt.csv', header = FALSE)
class_2_probs = read.csv('Data/new_aligned_shapesv3/V13_v2/gp2_spt.csv', header = FALSE)
ind = 3
old_veg1 = vcgImport(old_veg_files[ind])
new_veg1 = vcgImport(new_veg_files[ind])
new_veg_1 = process_off_file_v3(new_veg_files[ind])
old_veg_1 = process_off_file_v3(old_veg_files[ind])

mfrow3d(1,2)
new_veg_colors = rep('white', dim(new_veg1$vb)[2])
new_veg_colors[(which(class_2_probs > 0.05))] = 'red2'
plot3d(new_veg1,col = new_veg_colors, specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
rgl.viewpoint(userMatrix = rotation_matrix)

old_veg_colors = rep('white', dim(old_veg1$vb)[2])
old_veg_colors[(which(class_1_probs > 0.05))] = 'red2'

plot3d(old_veg1,col = old_veg_colors, specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
rgl.viewpoint(userMatrix = rotation_matrix)
rgl.postscript('~/Dropbox (DataPlusMath)/Sub-Image Analysis/Manuscript/Figures/figures/figure2/Figure_2c_v1.pdf',fmt = 'pdf')



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
                                              len = veg_pset$len,cuts = 300,cone_size = veg_pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)
veg_old_heat =  reconstruct_vertices_on_shape(dir = veg_pset$dirs,complex = old_veg_1,rate_vals = veg_veg_old$Rate2[,2],
                                              len = veg_pset$len,cuts = 300,cone_size = veg_pset$directions_per_cone,ball_radius = ball_radius, ball = ball, radius = 1)

#veg_new_colors=colfunc(max(veg_new_heat[,1]))[veg_new_heat[,1]]
#veg_old_colors=colfunc(max(veg_old_heat[,1]))[veg_old_heat[,1]]

#col_pal=c(color1,'lightblue','lightgreen',color5)
#colfunc <- colorRampPalette(col_pal)
veg_new_colors=colfunc(1 + max(veg_new_heat[,1]) - min(veg_new_heat[,1]))[1 + veg_new_heat[,1] - min(veg_new_heat[,1])]
veg_old_colors=colfunc(1 + max(veg_old_heat[,1]) - min(veg_old_heat[,1]))[1 + veg_old_heat[,1] - min(veg_old_heat[,1])]
mfrow3d(1,2)
plot3d(new_veg1,col = veg_new_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(new_veg1,new_veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)

plot3d(old_veg1,col = veg_old_colors, specular="black", axes = FALSE)
#plot_selected_landmarks(old_veg1,old_veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



class_1_probs = read.csv('Data/new_aligned_shapesv3/V13_v2/gp1_spt.csv', header = FALSE)
class_2_probs = read.csv('Data/new_aligned_shapesv3/V13_v2/gp2_spt.csv', header = FALSE)


roc_curve_caricature = compute_roc_curve_teeth(data_dir1 = old_veg_dir, data_dir2 = new_veg_dir, gamma = 0.5,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                             rate_values = veg_veg_old$Rate2[,2],directions_per_cone = veg_pset$directions_per_cone,curve_length = veg_pset$len,
                             directions = veg_pset$dirs,truncated = 300,ball_radius = ball_radius, ball = ball, radius = 1)

roc_curve_frame = data.frame(roc_curve_caricature)
ggplot(roc_curve_frame, aes(x = X1,y = X2,group = X3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X3) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Class') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "ROC Curve for Caricatured Teeth")

save.image('Figure2.Rdata')

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

mfrow3d(nr = 2, nc = 2)
plot3d(og_mesh,col = 'white', specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
rgl.viewpoint(userMatrix = rotation_matrix)
new_veg_colors[(which(class_2_probs > 0.05))] = 'red2'
plot3d(new_veg1,col = new_veg_colors, specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
rgl.viewpoint(userMatrix = rotation_matrix)
plot3d(og_mesh,col = 'white', specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
rgl.viewpoint(userMatrix = rotation_matrix)
old_veg_colors = rep('white', dim(old_veg1$vb)[2])
old_veg_colors[(which(class_1_probs > 0.05))] = 'red2'
plot3d(old_veg1,col = old_veg_colors, specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
rgl.viewpoint(userMatrix = rotation_matrix)
plot3d(old_veg1,col = veg_old_colors, specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
rgl.viewpoint(userMatrix = rotation_matrix)
