setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
library(R.utils)
#sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')
set.seed(55)
library(sinatra)
library(FNN)
library(rgl)
library(Rvcg)
library(plyr)
library(imager)

#Parameters for the Analysis
file = vcgImport(file = '~/Downloads/morphosourceMedia_11_04_19_010903/AM:M:10434_M38516-69993_Xeromys_myoides_Cranium.stl')
g = vcgQEdecim(file,percent = 0.05)
mfrow3d(1,2)
plot3d(file,color = 'white')
plot3d(g,color = 'white')
g = grayscale(file) > 0.5

plot(g)

mesh_image = function(image,threshold){
  
  matrix = image[,,1,1]
  vertices <- matrix(NA,nrow = 0, ncol = 4)
  edges <- matrix(NA, nrow = 0, ncol = 2)
  faces <- matrix(NA,nrow = 0, ncol = 3)
  length_x <- dim(matrix)[1]
  length_y <- dim(matrix)[2]
  
  for(i in 1:(length_x - 1)){
    for(j in 1:(length_y - 1)){
      
      if (matrix[i,j] == TRUE){
        vertex = c(i,j)
      }
      
      if ((i < (length_x-1))&(j<(length(y-1)))) 
      #The last three columns are the coordinates of the vertices
      #indices of the vertices of this pixel:
      vertex1 <- (i-1)*(dim(matrix)[2]) + j
      vertex2 <- (i - 1)*(dim(matrix)[2]) + j+1
      vertex3 <- (i)*(dim(matrix)[2]) + j
      vertex4 <- (i)*(dim(matrix)[2]) + j+1
      
      # add the vertices to complex
      vertices <- rbind(vertices, c( vertex1, i-length_x/2, j-length_y/2, matrix[i,j]) )
      vertices <- rbind(vertices, c( vertex2, i-length_x/2, j+1-length_y/2, matrix[i,j+1]))
      vertices <- rbind(vertices, c( vertex3, i+1-length_x/2, j-length_y/2, matrix[i+1,j] ))
      vertices <- rbind(vertices, c( vertex4, i+1-length_x/2, j+1-length_y/2, matrix[i+1,j+1] ))
      
      # add the edges to complex
      edges <- rbind(edges, c(vertex1,vertex2))
      edges <- rbind(edges, c(vertex2,vertex4))
      edges <- rbind(edges, c(vertex1,vertex3))
      edges <- rbind(edges, c(vertex3,vertex4))
      #edges<-rbind(edges,c(vertex2,vertex3))
      # add the faces to complex
      #faces <- rbind(faces, c(vertex1,vertex2,vertex3,vertex4))
      faces <- rbind(faces, c(vertex1,vertex2,vertex3))
      faces <- rbind(faces, c(vertex2,vertex3,vertex4))
    }
  }
  
  vertices <- unique(vertices)
  vertices = vertices[order(vertices[,1]),]
  edges <- unique(edges)
  faces <- unique(faces)
  complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
  complex$Vertices[,2]=complex$Vertices[,2]/(grid_length/2)
  complex$Vertices[,3]=complex$Vertices[,3]/(grid_length/2)
  complex$Vertices=complex$Vertices[,2:4]
  vertex=complex$Vertices
  faces=complex$Faces
  ind=rep(1,dim(vertex)[1])
  vertex=cbind(vertex,ind)
  mesh=tmesh3d(t(vertex),indices=t(faces))
  vertices=as.matrix(t(mesh$vb)[,1:3])
  faces=as.matrix(t(mesh$it))
  edges=vcgGetEdge(mesh)
  edges=as.matrix(edges[,1:2])
  complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
  return(complex)
  
}