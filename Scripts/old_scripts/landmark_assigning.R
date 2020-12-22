#### Aligning and Scaling the meshes ####
setwd('/Users/brucewang/Dropbox (Princeton)/Data + Experiments Tim Sudijono/')
library(Morpho)
library(rgl)
library(Rvcg)
dir='~/Documents/doug_new_teeth_raw/'


landmarks = read.csv(file = '~/Documents/doug_new_teeth_raw/landmarks_prelim.csv')
landmarks$museum_number = lapply(landmarks$museum_number, as.character)

n = 22
mesh = vcgImport(paste(dir,landmarks$museum_number[n][[1]], sep = ''))

landmark_coords = matrix(as.vector(landmarks[n,4:dim(landmarks)[2]]), ncol = 3, byrow = TRUE)


plot3d(mesh, col = 'white')

rgl.points(landmark_coords,size = 8, col = 'red')


rotation_matrix=matrix(c(0.99972576,0.02127766,0.00978078,0,0.01589701,-0.92330807,0.38373107,0,0.01719557,-0.38347024,-0.92339307,0,0,0,0,1),ncol=4,byrow=TRUE)


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



