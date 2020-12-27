#Post ETC alignments script

# Setup: Grid of 2918 directions 
# And discretization of 0(3) of size 16,000
# Prerequisities: The rotations corresponding to the permutations have been computed
library(Rvcg)

# Step 0: Load in the permutation-rotation setup:


rotations=readRDS("~/Dropbox (Princeton)/SINATRA_Data/ECT_alignment_demo/res/rots")
permutations=readRDS("~/Dropbox (Princeton)/SINATRA_Data/ECT_alignment_demo/res/perms")

dir=readRDS("~/Dropbox (Princeton)/SINATRA_Data/ECT_alignment_demo/postECTmiscdata/directions")
GRID=readRDS("~/Dropbox (Princeton)/SINATRA_Data/ECT_alignment_demo/postECTmiscdata/grid")

# Step 2: compute the euler curves (note directions already specified)

frobnorm=function(x){
  return(sqrt(sum(x^2)))
}

# A script for reading .OFF files
process_off_file_v3=function(input_dir){
  off=vcgImport(input_dir,silent = TRUE)
  vertices=as.matrix(t(off$vb)[,1:3])
  vertices=t(t(vertices)-colMeans(vertices))
  vertices=vertices/frobnorm(vertices)
  faces=as.matrix(t(off$it))
  edges=vcgGetEdge(off)
  edges=as.matrix(edges[,1:2])
  complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
  return(complex)
}


library(sinatra)
# Step 1: Compute the Euler curves
meshes=c()
#data_dir="~/Documents/research/SINATRAReview/TOBEALIGNED/doug_new_teeth_scaled/doug_new_teeth_scaled/"
data_dir="subsample/"
mesh_files <- list.files(file.path(data_dir), full.names = TRUE)
for (i in 1:length(mesh_files)) {
  mesh_orig=mesh_files[i]
  meshes=c(meshes,process_off_file_v3(mesh_orig))
  
}

### scale and center ###



bradius=0
for (i in 1:length(mesh_files)) {
  index=1+3*(i-1)
  coords=meshes[[index]]
  #coords=coords-rowMeans(coords)
  #coords=coords/sqrt(sum(coords^2))
  bradius=max(bradius,max(abs(coords)))
  print(bradius)
  meshes[[index]]=coords
}

curve_length=200
### Scale and center ends ###
data=list()
#data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]) )
for (i in 1:length(mesh_files)) {
  index1=1+3*(i-1)
  index2=3+3*(i-1)
  ec_curve1=matrix(ncol=curve_length,nrow=dim(dir)[1])
  complex1=meshes[index1:index2]
  for (j in 1:dim(dir)[1]){
    vertex_function1 <- complex1$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])
    curve1 <- compute_standardized_ec_curve(complex1, vertex_function1, curve_length-1, FALSE,ball_radius = bradius)[,2]
    ec_curve1[j,]<-curve1
    
  }
  data[[i]]<-ec_curve1
}

# Step 3: Obtain the post ECT alignment:



tmpN=length(data)
combo=combn(tmpN,2)
alignments=list()
indices=list()
distancematrix=matrix(rep(0,tmpN^2),ncol=tmpN,nrow=tmpN,byrow=T)
# To Do: Script that computes pairwise alignments, and the penalties, and the MST

for (i in 1:dim(combo)[2]) {
  record=dim(data[[1]])[1]*dim(data[[1]])[2]*1000
  index=0
  index1=combo[1,i]
  index2=combo[2,i]
  for (j in 1:length(permutations)) {
    d1=data[[index1]][permutations[[j]],]
    d2=data[[index2]]
    diff=sum((d1-d2)^2)
    if(diff<record){
      record=diff
      index=j
    }
    alignments[[i]]=rotations[[index]]
    distancematrix[index1,index2]<-distancematrix[index2,index1]<-record
  }
  
}


# These are the alignments obtained from Quasi BnB
# Load alignments from Quasi BnB:

R12=matrix(c(0.9028,-0.0085,0.43,0.1872,0.9079,-0.3751,-0.3872,0.4191,0.8212),ncol=3,nrow=3,byrow=T)
R13=matrix(c(0.9027,-0.0695,-0.4246,0.4248,0.3009,0.8538,-0.0684,0.9511,-0.3011),ncol=3,nrow=3,byrow=T)
R14=matrix(c(-0.2418,-0.907,0.3448,0.9692,-0.2087,0.1306,0.0464,-0.3657,-0.9296),ncol=3,nrow=3,byrow=T)
R15=matrix(c(-0.5984,0.6485,0.4705,0.6543,0.0567,0.7541,-0.4624,-0.7591,0.4582),ncol=3,nrow=3,byrow=T)
R23=matrix(c(0.6025,0.356,-0.7144,0.7235,0.1344,0.6771,-0.337,0.9248,0.1766),ncol=3,nrow=3,byrow=T)
R24=matrix(c(-0.0266,0.9502,0.3104,0.9358,0.1329,-0.3265,-0.3515,0.2817,-0.8928),ncol=3,nrow=3,byrow=T)
R25=matrix(c(-0.4096,0.2166,0.8862,0.91,0.0293,0.4135,-0.0636,-0.9758,0.2091),ncol=3,nrow=3,byrow=T)
R34=matrix(c(-0.4337,-0.257,-0.8636,0.7629,0.4052,-0.5037,0.4794,-0.8774,0.0203),ncol=3,nrow=3,byrow=T)
R35=matrix(c(-0.2648,0.9242,0.2751,-0.8639,-0.3541,0.3582,-0.4285,0.1428,-0.8922),ncol=3,nrow=3,byrow=T)
R45=matrix(c(0.4628,-0.7034,-0.5395,0.318,0.6999,-0.6396,-0.8274,-0.1244,-0.5476),ncol=3,nrow=3,byrow=T)


# Plot the aligned meshes. Top: Unaligned. Middle: post ECT aligned. Bottom: Quasi-BnB

#A2=alignments[[1]]
#A3=alignments[[2]]
#A4=alignments[[3]]
#A5=alignments[[4]]


mfrow3d(3,5)
plot3d(meshes[[1]])
plot3d(meshes[[4]])
plot3d(meshes[[7]])
plot3d(meshes[[10]])
plot3d(meshes[[13]])


plot3d(meshes[[1]])
plot3d(t(solve(alignments[[1]])%*%t(meshes[[4]])))
plot3d(t(solve(alignments[[2]])%*%t(meshes[[7]])))
plot3d(t(solve(alignments[[3]])%*%t(meshes[[10]])))
plot3d(t(solve(alignments[[4]])%*%t(meshes[[13]])))

plot3d(meshes[[1]])
plot3d(t(solve(R12)%*%t(meshes[[4]])))
plot3d(t(solve(R13)%*%t(meshes[[7]])))
plot3d(t(solve(R14)%*%t(meshes[[10]])))
plot3d(t(solve(R15)%*%t(meshes[[13]])))
