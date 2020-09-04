#### Aligning and Scaling the meshes ####
library(Morpho)
library(rgl)
library(Rvcg)
library(FNN)
library(R.utils)
library(sinatra)
sourceDirectory('~/gitrepos/SINATRA/SINATRA_Pipeline_Branch/')
dir='~/Documents/new_teeth_ect_aligned_by_species/Tarsius/'
# Mesh 1 has vertex 261
# Mesh 2 has vertex 295
# Mesh 3 has vertex 833
# Mesh 4 has vertex 224
# Mesh 5 has vertex 939
# Mesh 6 has vertex 278
# Mesh 7 has vertex 252
# Mesh 8 has vertex 593
# Mesh 9 has vertex 177
# Mesh 10 has vertex 357
# Mesh 11 has vertex 4845
# Mesh 12 has vertex 4215
# Mesh 13 has vertex 4619
# Mesh 14 has vertex 162
# Mesh 15 has vertex 244
# Mesh 16 has vertex 4618
# Mesh 17 has vertex 4225
# Mesh 18 has vertex 31
# Mesh 19 has vertex 4705
# Mesh 20 has vertex 4524
# Mesh 21 has vertex 219
# Mesh 22 has vertex 4301
# Mesh 23 has vertex 4898
# Mesh 24 has vertex 262
# Mesh 25 has vertex 4860
# Mesh 26 has vertex 69
# Mesh 27 has vertex 64
# Mesh 28 has vertex 4478
# Mesh 29 has vertex 4452
# Mesh 30 has vertex 4224
# Mesh 31 has vertex 102
# Mesh 32 has vertex 179
# Mesh 33 has vertex 4176
#Vertex 290,357 for Mesh 10 in Tarsius
rotation_matrix=matrix(c(0.8065218,0.5911149,0.01028626,0,0.5583186,-0.7558216,0.34207344,0,0.1944301,-0.2816325,-0.93961692,0,0,0,0,1),ncol=4,byrow=TRUE)
files = list.files(dir,full.names = TRUE)
n = 20
mesh = vcgImport(files[n])


#lm = find_landmarks(mesh, 10)

#rgl.points(coord,size = 8, col = 'red')

vertex = t(mesh$vb[-4,])

plot3d(mesh, col = 'white',alpha = 0.95)
#rgl.points(matrix(vertex[lm[2],],ncol = 3), size = 8, col = 'blue')
rgl.points(matrix(vertex[4524,],ncol = 3), size = 8, col = 'blue')
rgl.viewpoint(userMatrix = rotation_matrix)


g = knnx.index(vertex, query = matrix(vertex[lm[7],],ncol = 3),k=30)
plot3d(mesh, col = 'white')
rgl.points(matrix(vertex[g[4524],],ncol = 3), size = 8, col = 'blue')
rgl.viewpoint(userMatrix = rotation_matrix)




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



