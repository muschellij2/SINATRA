# Given a grid of elements in O(3) ("GRID"), and a set of n directions "dir")
# This script computes the optimal permutation p: [n] -> [n] for each provided element of O(3)
# These results are needed for post-ECT aligning the shape in postECTaligmnent.R
# NB. The provided grid massive and running the script may crash R. These may be computed in batches of
# 1000 elements. To do so, see the out commented lines after l. 108

# directions and grid are 2 pre-generated R dataframes that are described as follows:
# directions: Row for each direction, columns correspond to x,y, and z coordinates

# grid: A discretization of O(3). Each row of the dataframe is an element of O(3).
# The first 3 columns are angles about x,y and z coordinates. The last column indicates a reflection (1: no reflection, -1: reflection).

library(clue)
dir=readRDS("postECTmiscdata/directions")
GRID=readRDS("postECTmiscdata/grid")

# Some helper functions

cross<-function(x,y){
  a<-x
  a[1]<-x[2]*y[3]-x[3]*y[2]
  a[2]<-x[3]*y[1]-x[1]*y[3]
  a[3]<-x[1]*y[2]-x[2]*y[1]
  return(a)
}


computerotation=function(a,b){
  v=cross(a,b)
  c=sum(a*b)
  vx=matrix(c(0,-v[3],v[2],v[3],0,v[1],-v[2],v[1],0),ncol=3,nrow=3,byrow=T)
  rotation=diag(1,3)+vx+(1/(1+c))*(vx%*%vx)
  return(rotation)
}


anglenorm=function(R,Q){
  if(det(R)*det(Q)<0){
    return(100)
  }
  Qt=t(Q)
  P=R%*%Qt
  value = (sum(diag(P))-1)/2
  return(acos(value))
}

rotation = function(x,y){
  u=x/sqrt(sum(x^2))
  
  v=y-sum(u*y)*u
  if(sum(abs(v))<10^-12){
    return(diag(1,3))
  }
  v=v/sqrt(sum(v^2))
  
  cost=sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2))
  
  sint=sqrt(1-cost^2);
  
  diag(length(x)) - u %*% t(u) - v %*% t(v) + 
    cbind(u,v) %*% matrix(c(cost,-sint,sint,cost), 2) %*% t(cbind(u,v))
}


# Helper functions end

# apply functions



foorotation=function(foo){
  theta=foo[1]
  phi=foo[2]
  psi=foo[3]
  refl=foo[4]
  c=cos(theta)
  s=sin(theta)
  A=matrix(c(c,-s,0,s,c,0,0,0,1),nrow=3,ncol=3,byrow=T)
  c=cos(phi)
  s=sin(phi)
  B=matrix(c(c,0,s,0,1,0,-s,0,c),nrow=3,ncol=3,byrow=T)
  c=cos(psi)
  s=sin(psi)
  C=matrix(c(1,0,0,0,c,-s,0,s,c),nrow=3,ncol=3,byrow=T)
  D=diag(refl,3)
  R=D%*%C%*%B%*%A
  return(R)
}

foopermutationB=function(foo){
  tmpdir<<-t(foo%*%t(dir))
  NN<<-dim(dir)[1]
  AA<<-expand.grid(1:NN,1:NN)
  vex<<-apply(AA,FUN=fun2,MARGIN=1)
  XX<<-matrix(vex,ncol=NN,nrow=NN)
  perm1<<-solve_LSAP(XX)
  return(perm1[1:NN])
  
}


fun2=function(kk){
  return(sum((dir[kk[1],]-tmpdir[kk[2],])^2))
}

# apply functions end

rotations=list()
permutations=list()

#NB. Computing this is heavy and will likely crash if done in one go 
# Consider Subsetting GRID to 1000 entries or so at a time
GRID=GRID
# Subset to first 1000 entries
#GRID=GRID[1:1000,]
NN=dim(dir)[1]
AA=expand.grid(1:NN,1:NN)
library(parallel)
TMPGRID=split(GRID, rep(1:nrow(GRID)))
tmpGRID=lapply(TMPGRID,unlist)
list_of_rotations=lapply(tmpGRID,foorotation)

start_time <- Sys.time()
perms<<-mclapply(list_of_rotations,foopermutationB, mc.cores = detectCores())
end_time <- Sys.time()
print(end_time-start_time)

saveRDS(perms,"perms")
saveRDS(list_of_rotations,"rots")

#saveRDS(perms,"perms01000")
#saveRDS(list_of_rotations,"rots01000")


