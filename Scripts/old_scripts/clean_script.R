#### Aligning and Scaling the meshes ####
setwd('~/Dropbox/Data + Experiments Tim Sudijono/')
library(Morpho)
library(rgl)
library(Rvcg)
dir='~/Documents/face_meshes/off_converted_meshes/front_angry/'
out_path = '~/Documents/face_meshes/off_cleaned_meshes/front_angry/'

clean_face_files=function(input_dir,output_dir){
  files=list.files(path = dir,full.names = TRUE)
  filenames=list.files(path=dir,full.names = FALSE)
  num_files=length(files)
  for (i in 1:num_files){
    print(files[i])
    file=vcgImport(files[i])
    file = vcgIsolated(file)
    file = vcgClean(file, sel = c(0,1,2,3,4,5,6,7), iterate = TRUE)
    filename=paste(out_path,filenames[i],sep='')
    area=vcgArea(file)
    file2=scalemesh(file,1/sqrt(area),center='mean')
    centroid <- colMeans(vert2points(file2))
    file3=translate3d(file2,-centroid[1],-centroid[2],-centroid[3])
    print(filename)
    vcgOffWrite(file3,filename = filename)
  }
}
clean_face_files(dir, out_path)


new_files = list.files(out_path, full.names = TRUE)
d = vcgImport(new_files[5])
#m = vcgBallPivoting(radius = 0.02, d)

mfrow3d(3,10)
for (i in 121:150){
  d = vcgImport(new_files[i])
  plot3d(d, col = skin1, axes = FALSE)
  rgl.viewpoint(userMatrix = m)
}

library(auto3dgm)
library(rgl)
library(Rvcg)

Data_dir = "/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/Data/HDM/all_scaled"
Output_dir = "/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/Data/all_files_alignedv4/"
#Data_dir = '~/Documents/face_meshes/to_be_aligned_meshes'
#Output_dir = "~/Documents/face_meshes/newly_aligned_meshes"
#Output_dir2  =  "~/Documents/face_meshes/newly_aligned_meshes/Aligned_Shapes/"
rotation_matrix=matrix(c(0.99972576,0.02127766,0.00978078,0,0.01589701,-0.92330807,0.38373107,0,0.01719557,-0.38347024,-0.92339307,0,0,0,0,1),ncol=4,byrow=TRUE)


#ids1 =  list.files(paste(Output_dir2),full.names = TRUE)
ids1 =  list.files(paste(Data_dir),full.names = TRUE)

Ids = list.files(paste(Data_dir,'lowres',sep = ''))
Ids = list.files(Data_dir)
Names =  list.files(Data_dir)
#Ids = c('AMNH-M-211491_M818')
#Ids = Ids[-98]
#Ids = Ids[-74]
for (i in 1:length(Ids)){
  Ids[i] = gsub("\\..*","",Ids[i])
}
Names = Ids

Levels=c(96,192)
#FULL is a list of 3 returned elements.  The user gets to specify what is returned.
FULL = align_shapes(Data_dir, Output_dir, Levels, Ids, Names)
mfrow3d(nr=5,nc = 10)
for (i in 1:50){
  file = vcgImport(ids1[i])
  plot3d(file, color = 'white')
  rgl.viewpoint(userMatrix =rotation_matrix)
}



