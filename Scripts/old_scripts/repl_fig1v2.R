setwd('/Users/brucewang/Dropbox (Princeton)/Data + Experiments Tim Sudijono/')
library(R.utils)
#sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')
set.seed(55)
library(sinatra)
library(FNN)
library(rgl)
library(Rvcg)
library(plyr)

#Parameters for the Analysis
pset1 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Microcebus/',dir2 = '~/Documents/doug_new_teeth_by_species/Mirza/',
             num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
                dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset2 =  list(dir1 = '~/Documents/doug_new_teeth_by_species/Microcebus/',dir2 = '~/Documents/doug_new_teeth_by_species/Tarsius/',
              num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
              dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset3 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Microcebus/',dir2 = '~/Documents/doug_new_teeth_by_species/Saimiri/',
             num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset4 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Tarsius/',dir2 = '~/Documents/doug_new_teeth_by_species/Saimiri/',
             num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset5 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Tarsius/',dir2 = '~/Documents/doug_new_teeth_by_species/Mirza/',
             num_cones = 35, cap_radius = 0.25, len = 75, directions_per_cone = 5,
             dirs = generate_equidistributed_cones(num_directions = 35, cap_radius =  0.25, directions_per_cone = 5))
pset6 = list(dir1 = '~/Documents/doug_new_teeth_by_species/Saimiri/',dir2 = '~/Documents/doug_new_teeth_by_species/Mirza/',
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

load('chris_henry_analysis_2019_07_02.Rdata')
#### Start Comparison ####

#### New Fruit Veg ####

k = 10
mesh = vcgImport(list.files(pset2$dir2,full.names = TRUE)[k])
plot3d(mesh, col = c2_mesh2_colors_list[[k]],back="lines", specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
rgl.viewpoint(userMatrix = rotation_matrix)
mesh$vb[1:3,] = mesh$vb[1:3,]*10
mesh$material = list(color = rep('white',dim(mesh$vb)[2]))
mesh$material = list(color = c2_mesh2_colors_list[[k]])
plot3d(mesh)
vcgWrlWrite(mesh,filename = 'tooth_prototype_wrlv3')
g = col2rgb(c2_mesh2_colors_list[[k]])
g = g/256
write.csv(t(g), col.names = FALSE, row.names = FALSE, file = 'colorsv2.csv')
fruit = vcgImport('tooth_prototype_wrl.wrl',readcolor = TRUE)
plot3d(fruit)
plot3d(fruit1,col = fruit_colors_heat, back="lines", specular="white", axes = FALSE,xlab = '', ylab = '',zlab='')
#plot_selected_landmarks(fruit1,fruit1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)



plot3d(veg1,col = veg_colors_heat, back="lines", specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
#plot_selected_landmarks(veg1,veg1_vert,num_landmarks = 20)
rgl.viewpoint(userMatrix = rotation_matrix)
rgl.snapshot(filename = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Figures/Figure Drafts//figure1/figure1e_v6.png',fmt = 'png')

col.table=colfunc(1 + max(fruit_heat[,1]) - min(fruit_heat[,1]))
color.bar(col.table,min=10000 * round(min(fruit_heat[,2]) - min(fruit_heat[,2]),4),max=10000 * round(max(fruit_heat[,2])-min(fruit_heat[,2]),4))
dev.copy2pdf(file='~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Figures/Figure Drafts/figure1/figure_1d_tooth1_colorbar_v6.pdf',out.type = 'pdf')
col.table=colfunc(1 + max(veg_heat[,1]) - min(veg_heat[,1]))
color.bar(col.table,min=10000*round(min(veg_heat[,2]) - min(veg_heat[,2]),4),max=10000*round(max(veg_heat[,2])-min(veg_heat[,2]),4))
dev.copy2pdf(file='~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Figures/Figure Drafts/figure1/figure_1d_tooth2_colorbar_v6.pdf',out.type = 'pdf')


save.image('Figure1.Rdata')
