setwd('/Users/timothysudijono/Dropbox/Data + Experiments/')
library(R.utils)
library(Rcpp)
library(Rvcg)
library(rgl)
library(pracma)
library(expm)
sourceDirectory('~/projects/Research/SINATRA/SINATRA_Pipeline_Branch/')
sourceCpp('~/projects/Research/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')
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
for (k in 101:103){
  normal_dir <- apply(dirs[101:103,],2, sum)/3
  arrow3d(-0.75*normal_dir,0.75*dirs[k,],s = s, type = arrow_type,width= width,thickness = thickness,col = 'blue')
}
for (k in 91:93){
  normal_dir <- apply(dirs[91:93,],2, sum)/3
  arrow3d(-0.75*normal_dir,0.75*dirs[k,],s = s, type = arrow_type,width= width,thickness = thickness, col = 'red')
}
for (k in 56:58){
  normal_dir <- apply(dirs[56:58,],2, sum)/3
  arrow3d(-0.75*normal_dir,0.75*dirs[k,],s = s, type = arrow_type,width= width,thickness = thickness, col = 'green')
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





### reconstruction
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  2*cap_radius, directions_per_cone = directions_per_cone)
arrow_type = 'rotation'
s = 1/10
width = 1/9
thickness = 1/20
sublevel_size <- 0.05
idxs_of_intersection <- list()
normal_dir <- apply(dirs[101:103,],2, sum)/3
for (k in 101:103){
  vec <- 0.75*dirs[k,] + 0.75*normal_dir
  idxs_of_intersection[[k - 100]] <- computeIntersectionsWithPlane(vec[1],vec[2],vec[3], c(0,0,0), sublevel_size, fruit1)
}

# Color the intersection
common_idxs <- Reduce(intersect,idxs_of_intersection)
cols = rep('white', dim(fruit1$vb)[2])
cols[common_idxs] = 'red'
plot3d(fruit1,col = cols,  specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')

for (k in 101:103){
  vec <- 0.75*dirs[k,] + 0.75*normal_dir
  offset <- vec - normal_dir*(vec*normal_dir)/(sqrt(sum(normal_dir^2)))
  arrow3d(-0.75*normal_dir + offset,0.75*dirs[k,] + offset,s = s, type = arrow_type,width= width,thickness = thickness,col = 'dark red')
  a <- 0.75*dirs[k,1] + 0.75*normal_dir[1]
  b <- 0.75*dirs[k,2] + 0.75*normal_dir[2]
  c <- 0.75*dirs[k,3] + 0.75*normal_dir[3]
  d <- c(0)
  
  addPlaneWithThickness(vec[1],vec[2],vec[3],c(0,0,0), sublevel_size, color = 'pink', 0.15)
}
rgl.viewpoint(userMatrix = rotation_matrix)
aspect3d("iso")


# another direction
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  2*cap_radius, directions_per_cone = directions_per_cone)
arrow_type = 'rotation'
s = 1/10
width = 1/9
thickness = 1/20
sublevel_size <- 0.075
idxs_of_intersection <- list()
for (k in 56:58){
  vec <- 0.75*dirs[k,] + 0.75*normal_dir
  idxs_of_intersection[[k - 55]] <- computeIntersectionsWithPlane(vec[1],vec[2],vec[3], c(0,0,0), sublevel_size, fruit1)
}

# Color the intersection
common_idxs <- Reduce(intersect,idxs_of_intersection)
cols = rep('white', dim(fruit1$vb)[2])
cols[common_idxs] = 'green'
plot3d(fruit1,col = cols,  specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')

normal_dir <- apply(dirs[56:58,],2, sum)/3
for (k in 56:58){
  vec <- 0.75*dirs[k,] + 0.75*normal_dir
  offset <- vec - normal_dir*(vec*normal_dir)/(sqrt(sum(normal_dir^2)))
  arrow3d(-0.75*normal_dir,0.75*dirs[k,],s = s, type = arrow_type,width= width,thickness =sublevel_size,color = 'dark green')
  a <- 0.75*dirs[k,1] + 0.75*normal_dir[1]
  b <- 0.75*dirs[k,2] + 0.75*normal_dir[2]
  c <- 0.75*dirs[k,3] + 0.75*normal_dir[3]
  d <- c(0)
  
  addPlaneWithThickness(vec[1],vec[2],vec[3], c(0,0,0), sublevel_size, 'light green',0.2)
  idxs_of_intersection[[k - 55]] <- computeIntersectionsWithPlane(vec[1],vec[2],vec[3], c(0,0,0), sublevel_size, fruit1)
}

rgl.viewpoint(userMatrix = rotation_matrix)
aspect3d("iso")

# another direction
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  2*cap_radius, directions_per_cone = directions_per_cone)
arrow_type = 'rotation'
s = 1/10
width = 1/9
thickness = 1/20
sublevel_size <- 0.05
idxs_of_intersection <- list()
normal_dir <- apply(dirs[91:93,],2, sum)/3
for (k in 91:93){
  vec <- 0.75*dirs[k,] + 0.75*normal_dir
  idxs_of_intersection[[k - 90]] <- computeIntersectionsWithPlane(vec[1],vec[2],vec[3], c(0,0,0), sublevel_size, fruit1)
}

# Color the intersection
common_idxs <- Reduce(intersect,idxs_of_intersection)
cols = rep('white', dim(fruit1$vb)[2])
cols[common_idxs] = 'red'
plot3d(fruit1,col = cols,  specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')

for (k in 91:93){
  vec <- 0.75*dirs[k,] + 0.75*normal_dir
  offset <- vec - normal_dir*(vec*normal_dir)/(sqrt(sum(normal_dir^2)))
  arrow3d(-0.75*normal_dir,0.75*dirs[k,],s = s, type = arrow_type,width= width,thickness =sublevel_size,color = 'dark red')
  a <- 0.75*dirs[k,1] + 0.75*normal_dir[1]
  b <- 0.75*dirs[k,2] + 0.75*normal_dir[2]
  c <- 0.75*dirs[k,3] + 0.75*normal_dir[3]
  d <- c(0)
  
  addPlaneWithThickness(vec[1],vec[2],vec[3], c(0,0,0), sublevel_size, 'pink',0.2)
}

rgl.viewpoint(userMatrix = rotation_matrix)
aspect3d("iso")

### Sequence of figures ###
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  2*cap_radius, directions_per_cone = directions_per_cone)
arrow_type = 'rotation'
s = 1/15
width = 1/12
thickness = 1/40
sublevel_size <- 0.05
idxs_of_intersection <- list()
normal_dir <- apply(dirs[56:58,],2, sum)/3

for (k in 56:58){
  vec <- 0.75*dirs[k,] + 0.75*normal_dir
  idxs_of_intersection[[k - 55]] <- computeIntersectionsWithPlane(vec[1],vec[2],vec[3], c(0,0,0), sublevel_size, fruit1)
}

rgl.viewpoint(userMatrix = rotation_matrix,zoom = 1)
aspect3d("iso")

# Color the intersection
common_idxs <- Reduce(intersect,idxs_of_intersection)
cols1 = rep('white', dim(fruit1$vb)[2])
cols1[common_idxs] = 'green'
plot3d(fruit1,col = 'white',  specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')

for (k in 56:58){
  arrow3d(-0.75*normal_dir,0.75*dirs[k,],s = s, type = arrow_type,width= width,thickness =sublevel_size,color = 'dark green')
}

for (k in 56:58){
  vec <- 0.75*dirs[k,] + 0.75*normal_dir
  addPlaneWithThickness(vec[1],vec[2],vec[3], c(0,0,0), sublevel_size, 'light green',0.2)
}

#### Overlap colors
idxs_of_intersection <- list()
normal_dir <- apply(dirs[91:93,],2, sum)/3

for (k in 91:93){
  vec <- 0.75*dirs[k,] + 0.75*normal_dir
  idxs_of_intersection[[k - 90]] <- computeIntersectionsWithPlane(vec[1],vec[2],vec[3], c(0,0,0), sublevel_size, fruit1)
}
common_idxs <- Reduce(intersect,idxs_of_intersection)

cols2 = cols
cols2[common_idxs] = 'red'
plot3d(fruit1,col = cols2,  specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')

rgl.viewpoint(userMatrix = rotation_matrix)
aspect3d("iso")

##### Helper Functions #####
##### plot intersections on the teeth #####
computeIntersectionsWithPlane <- function(a,b,c, center, thickness, complex){
  vertices <- complex$vb[1:3,]
  v <- c(a,b,c)
  projs <- colSums(v*vertices)
  d <- dot(v,center)
  idxs <- which(projs > d & projs < d + thickness)
  idxs
}

##### Add planes #####
addPlaneWithThickness <- function(a,b,c,center,thickness,color,alpha){
  bbox <- par3d("bbox")
  slab <- extrude3d(bbox[c(1,2,2,1)], bbox[c(3,3,4,4)], thickness)
  
  dir <- c(a,b,c)
  dir <- dir/sqrt(sum(dir^2))
  
  if(dir[3] < 0){
    dir <- -dir
  }
  
  r <- sqrt(dir[1]^2 + dir[2]^2)
  x <- dir[1]
  y <- dir[2]
  z <- dir[3]
  
  phi <- acos(z)
  theta <- acos(x/r)
  
  # accomodate arccos and quadrant
  if(y < 0){
    theta <- 2*pi - theta
  }
    
  slab <- rotate3d(slab, -phi, 0,1,0)
  slab <- rotate3d(slab, -theta, 0,0,1)
  slab <- translate3d(slab,center[1],center[2],center[3])
  
  shade3d(slab, col = color,alpha = alpha)
}

##### Testing Funcs ######
#Test plane plotting
thickness = 0.1
center = c(0,0,0)
example(plot3d)
a <- -1
b <- 3
c <- -1
dir <- c(a,b,c)
dir <- dir/sqrt(sum(dir^2))
addPlaneWithThickness(a,b,c,center,thickness)
arrow3d(-2*dir,2*dir, s = s, type = 'rotation',width= width,thickness = thickness,col = 'black')
aspect3d("iso")







