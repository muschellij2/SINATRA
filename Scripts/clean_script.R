#### Aligning and Scaling the meshes ####
setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
library(Morpho)
library(rgl)
library(Rvcg)
dir='~/Documents/face_meshes/off_converted_meshes/front_angry/'
out_path = '~/Documents/face_meshes/off_cleaned_meshes/front_angry/'

b = list.files(dir, full.names = TRUE)

mesh1 = vcgImport(b[8])
plot3d(mesh1, col = 'white')
cleanmesh = vcgIsolated(mesh1)
plot3d(cleanmesh, col = 'white')
cleanmesh = vcgClean(cleanmesh, sel = c(0,1,2,3,4,5,7,0,1,2,3,4,5,6,7, 0,1,2,3,4,5,6,7))
plot3d(cleanmesh,col = 'white')

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

want_files = new_files[c(4,6,7,8,9,11,13,14,16,18,19,20,21,23,25,26,27,28,29,30,31,32,33,36,39,41,42,47,48,49,50,51,52,53,55,59,65,67,68,70,75,76,77,78,79,82,83,85,86,88,91,93,94,96,97,98,101,104,106,108,109,111,113,114,118,119,120,122,123,129,130,136,138,141,144)]
new_want_files = want_files[c(1,2,3,4,5,6,12,14,15,16,17,18,19,20,22,23,24,27,28,31,32,33,34,35,36,37,38,39,40,41,43,45,46,48,49,52,54,55,59,60,62,64,65,66,70,71,75)]
length(new_want_files)

mfrow3d(3,5)
for (i in 1:15){
  d = vcgImport(new_want_files[i])
  plot3d(d, col = skin1, axes = FALSE)
  rgl.viewpoint(userMatrix = m)
}

plot3d(d, col = 'white')
plot3d(m, col = 'white')

dir='~/Documents/face_meshes/off_converted_meshes/front_smiling/'
out_path = '~/Documents/face_meshes/off_cleaned_meshes/front_smiling/'
clean_face_files(dir, out_path)
new_files = list.files(out_path, full.names = TRUE)
want_files = new_files[c(1,4,6,7,9,11,12,17,20,21,23,24,26,27,28,31,36,41,42,51,57,58,64,65,74,75,79,80,82,83,84,85,86,88,90,91,93,96,97,99,100,103,106,108,110,111,112,114,118,119)]
length(want_files)
new_want_files = want_files[-c(10,37,41,45)]
write(new_want_files,'smiling_faces.txt',sep = ',')
mfrow3d(3,10)
for (i in 1:30){
  d = vcgImport(new_want_files[i])
  plot3d(d, col = skin1, axes = FALSE)
  rgl.viewpoint(userMatrix = m)
}


library(auto3dgm)
library(rgl)
library(Rvcg)

#Data_dir = "/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/Data/all_files_tingran_scaled/"
#Output_dir = "/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/Data/all_files_alignedv2"
Data_dir = '~/Documents/face_meshes/to_be_aligned_meshes'
Output_dir = "~/Documents/face_meshes/newly_aligned_meshes"
Output_dir2  =  "~/Documents/face_meshes/newly_aligned_meshes/Aligned_Shapes/"
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



