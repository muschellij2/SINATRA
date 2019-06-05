#### Aligning and Scaling the meshes ####
setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
library(Morpho)
library(rgl)
library(Rvcg)
library(FNN)
library(R.utils)
sourceDirectory('~/Documents/SINATRA/SINATRA_Pipeline_Branch/')
dir='~/Documents/doug_new_teeth_by_species/Tarsius/'
#Vertex 290,357 for Mesh 10 in Tarsius
rotation_matrix=matrix(c(0.8065218,0.5911149,0.01028626,0,0.5583186,-0.7558216,0.34207344,0,0.1944301,-0.2816325,-0.93961692,0,0,0,0,1),ncol=4,byrow=TRUE)
files = list.files(dir,full.names = TRUE)
n = 10
mesh = vcgImport(files[n])
plot3d(mesh, col = 'white')

rgl.viewpoint(userMatrix =rotation_matrix)


plot3d(mesh, col = 'white')
coord = matrix(c(-0.300355,0.180715,-0.147777),ncol = 3)

rgl.points(coord,size = 8, col = 'red')

vertex = t(mesh$vb[-4,])

indices = knnx.index(vertex,coord,k= 1)
rgl.points(matrix(vertex[indices,],ncol = 3),size = 10, col = 'blue')
b = knnx.index(vertex,matrix(vertex[indices,],ncol = 3),200)
cols = rep('white',dim(vertex)[1])
cols[b] = 'red'
plot3d(mesh, col = cols)

lm = find_landmarks(mesh, num_landmark = 30)

plot3d(mesh, col = 'white')
rgl.points(matrix(vertex[lm[22],],ncol = 3), size = 8, col = 'blue')





#ids1 =  list.files(paste(Output_dir2),full.names = TRUE)
ids1 =  list.files(paste(Data_dir),full.names = TRUE)

#Ids = list.files(paste(Data_dir,'lowres',sep = ''))
Ids = list.files(Data_dir)
Names =  list.files(Data_dir)
#Ids = c('AMNH-M-211491_M818')
#Ids = Ids[-98]
#Ids = Ids[-74]
for (i in 1:length(Ids)){
  Ids[i] = gsub("\\..*","",Ids[i])
}

Ids = Ids[-54]
Names = Ids

Levels=c(96,192)
#FULL is a list of 3 returned elements.  The user gets to specify what is returned.
FULL = align_shapes(Data_dir, Output_dir, Levels, Ids, Names)
mfrow3d(nr=5,nc = 10)
for (i in 1:50){
  file = vcgImport(ids1[i])
  plot3d(file, color = 'white', axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
  rgl.viewpoint(userMatrix =rotation_matrix)
}
ids2 = list.files(paste(Output_dir,'/Aligned_Shapes', sep = ''),full.names = TRUE)
mfrow3d(nr=5,nc = 10)
for (i in 1:50){
  file = vcgImport(ids2[i])
  plot3d(file, color = 'white', axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
  rgl.viewpoint(userMatrix =rotation_matrix)
}



