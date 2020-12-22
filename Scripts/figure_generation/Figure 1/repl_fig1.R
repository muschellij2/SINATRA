setwd('/Users/brucewang/Dropbox (Princeton)/Data + Experiments Tim Sudijono/')
library(R.utils)
library(Rcpp)
sourceDirectory('~/Documents/SINATRA/SINATRA_Pipeline_Branch/')
sourceCpp('~/Documents/SINATRA/SINATRA_Pipeline_Branch/BAKRGibbs.cpp')
set.seed(55)
load('Figure1.Rdata')

#Specifying Directories and the Associated Files


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
color1='blue'
color2='lightgreen'
color3='orangered'
color3 = 'red'
col_pal=c(color1,color2,color2,color3)
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
mfrow3d(2,1)
fruit_colors = rep('white', dim(fruit1$vb)[2])

plot3d(fruit1,col = fruit_colors, back="lines", specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
#plot_selected_landmarks(fruit1,fruit1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)

veg_colors = rep('white', dim(veg1$vb)[2])
plot3d(veg1,col = veg_colors, back="lines", specular="black", axes = FALSE,xlab = '', ylab = '',zlab='')
#plot_selected_landmarks(veg1,veg1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)
rgl.snapshot(filename = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Figures/Figure Drafts/figure1/figure1a_v3.png',fmt = 'png')
rgl.postscript(filename = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Figures/figures/figure1/figure1a_v2.pdf',fmt = 'pdf')
#1b
arrow_type = 'rotation'
s = 1/10
width = 1/9
thickness = 1/20
plot3d(fruit1,col = fruit_colors,  specular="black", axes = FALSE,xlab = '', ylab = '',zlab='',alpha = 0.95)
for (k in 101:105){
  arrow3d(-0.5*dirs[101,],0.5*dirs[k,],s = s, type = arrow_type,width= width,thickness = thickness,col = 'blue')
}
for (k in 91:95){
  arrow3d(-0.5*dirs[90,],0.5*dirs[k,],s = s, type = arrow_type,width= width,thickness = thickness, col = 'red')
}
for (k in 56:60){
  arrow3d(-0.5*dirs[56,],0.65*dirs[k,],s = s, type = arrow_type,width= width,thickness = thickness, col = 'green')
}
#arrow3d(0.5*-dirs[1,],0.5*dirs[1,],s = 1/9, type = 'rotation',width= 1/9,thickness = 1/20)
#plot_selected_landmarks(fruit1,fruit1_vert)
rgl.viewpoint(userMatrix = rotation_matrix)
rgl.snapshot(filename = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Figures/Figure Drafts/figure1/figure1b_v4.png',fmt = 'png')
rgl.postscript(filename = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Figures/figures/figure1/figure1b_v1.pdf',fmt = 'pdf')

#1c

#### Create the Filtartion ####

dir_index = 101
projections <- fruit_1$Vertices[,1:3]%*%dirs[dir_index,]
buckets <- seq(-ball_radius,ball_radius,length.out = len+1)

#bucket these projections into curve_length number of groups; could have also solved this with the cut function
step_length <- (max(buckets) - min(buckets))/len
#Replace projections by buckets
projection_buckets <- apply((projections - min(buckets))/step_length,1, function(float) as.integer(float)) + (len+1)*(dir_index-1)
#print(step_l)
projection_buckets=projection_buckets+1
# print(paste(min(projection_buckets), max(projection_buckets)))
selected_vertices[[i]]=which(projection_buckets %in% indices)
# print(paste(min(projection_buckets), max(projection_buckets)))
min(projection_buckets)
max(projection_buckets)


total_projections = cbind(projections,projection_buckets)
mfrow3d(1,5)
s = 1/6
width = 1/6
thickness = 1/20

for (k in seq(min(projection_buckets),max(projection_buckets),length.out = 5)){
  alphas = rep(0, dim(fruit1$vb)[2])
  filtration_vertex = which(projection_buckets < k)
  alphas[filtration_vertex] = 0.97
  if (k == max(projection_buckets)){
    alphas[filtration_vertex] = 1
  }
  
  filtered_projections = total_projections[which(total_projections[,2] < k+1),]
  projection = 1.25*max(filtered_projections[,1])
  colors = rep('white', dim(fruit1$vb)[2])
  plot3d(fruit1,col = colors,  specular="black", axes = FALSE,xlab = '', ylab = '',zlab='',alpha = alphas)
  arrow3d(-0.5*dirs[dir_index,],0.75*dirs[dir_index,],s = s, type = arrow_type,width= width,thickness = thickness,col = 'red')
  rgl.viewpoint(userMatrix = rotation_matrix)
}
rgl.snapshot(filename = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Figures/Figure Drafts/figure1/figure1c_v3.png',fmt = 'png')
#rgl.postscript(filename = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Figures/figures/figure1/figure1c_v1.pdf',fmt = 'pdf')
rgl.snapshot(filename = '~/Dropbox (Princeton)/Sub-Image Analysis/Manuscript/Figures/figures/figure1/figure1c_v1.png',fmt = 'png')

#EC curve
dir_index = 6
projections <- fruit_1$Vertices[,1:3]%*%dirs[dir_index,]
projections2 <- veg_1$Vertices[,1:3]%*%dirs[dir_index,]

ec_curve = compute_standardized_ec_curve(complex = fruit_1, vertex_function = projections,curve_length = len,first_column_index = FALSE ,ball_radius = ball_radius)
ec_curve2 = compute_standardized_ec_curve(complex = veg_1, vertex_function = projections2,curve_length = len,first_column_index = FALSE ,ball_radius = ball_radius)

#ec_curve = integrate_ec_curve(ec_curve)
#ec_curve2 = integrate_ec_curve(ec_curve2)
par(oma = c(4, 1, 1, 1))
par(mfrow  = c(1,5))

t=seq(0,10,0.1)
y=sin(t)*3
x=sin(t)*3*cos(t)
y[1:30] = max(y[1:30]) 
y[16:30] = max(y[16:30]) 
y[31:60] = min(y[31:60]) 
y[46:60] = min(y[46:60]) 
y[61:90] = max(y[61:90]) 
y[76:101] = max(y[76:101]) 
y[91:101] = min(y[91:101]) 


x[1:15] = max(x[1:15]) 
x[16:30] = min(x[16:30]) 
x[31:45] = max(x[31:45]) 
x[46:60] = min(x[46:60]) 
x[61:75] = max(x[61:75]) 
x[76:90] = min(x[76:90]) 
x[91:101] = max(x[91:101]) 

ec_curve = cbind(t,y)
ec_curve2 = cbind(t,x)

sig = c(t[15:25],t[80:88])
#sig = c(ec_curve[,1][10:20],ec_curve[,1][44:54],ec_curve[,1][69:83])
par(oma = c(4, 1, 1, 1))
par(mfrow  = c(1,5))

for (k in seq(min(ec_curve[,1]),max(ec_curve[,1]),length.out = 5)){
  curve = ec_curve[which(ec_curve[,1]<=k),]
  curve2 = ec_curve2[which(ec_curve2[,1]<=k),]
  if (k == 0){
    plot(c(0,0),type = 'l',xlab="Filtration Steps", ylab="Topological Summary Statistics",bty="n",lwd=3,ylim=c(-4,4))
  }
  else{
  plot(curve,type="l", xlab="Filtration Steps", ylab="Topological Summary Statistics",bty="n",lwd=3,ylim=c(-4,4))
  lines(curve2, col = 'black',lwd=3,lty=2) 
  }
}
plot(curve,type="l", xlab="Filtration Steps", ylab="Topological Summary Statistics",bty="n",lwd=3,ylim=c(-4,4))
lines(curve2, col = 'black',lwd=3,lty=2) 
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top",legend = c("Species #1","Species #2"),lty=c(1,2),lwd=2, bty = "n",horiz = TRUE)

for (k in seq(min(ec_curve[,1]),max(ec_curve[,1]),length.out = 5)){
  curve = ec_curve[which(ec_curve[,1]<=k),]
  curve2 = ec_curve2[which(ec_curve2[,1]<=k),]
  if (k == 0){
    plot(c(0,0),type = 'l',xlab="Filtration Steps", ylab="Topological Summary Statistics",bty="n",lwd=3,ylim=c(-4,4))
  }
  else{
  plot(curve,type="l", xlab="Filtration Steps", ylab="Topological Summary Statistics",bty="n",lwd=3,ylim=c(-4,4))
  if (length(curve)> 2){
    for(i in 1:(dim(curve)[1]-1)){
      segments(curve[,1][i],curve[,2][i],curve[,1][i+1],curve[,2][i+1],col=ifelse(curve[i,1]%in%sig,"red","black"),lwd = 3)
    }
  }
  lines(curve2, col = 'black',lwd=3,lty=2) 
  if (length(curve)> 2){
    for(i in 1:(dim(curve2)[1]-1)){
      segments(curve2[,1][i],curve2[,2][i],curve2[,1][i+1],curve2[,2][i+1],col=ifelse(curve2[i,1]%in%sig,"red","black"),lwd = 3,lty=3)
    }
  }
  }
}

plot(curve,type="l", xlab="Filtration Steps", ylab="Topological Summary Statistics",bty="n",lwd=3,ylim=c(-4,4))
if (length(curve)> 2){
  for(i in 1:(dim(curve)[1]-1)){
    segments(curve[,1][i],curve[,2][i],curve[,1][i+1],curve[,2][i+1],col=ifelse(curve[i,1]%in%sig,"red","black"),lwd = 3)
  }
}
lines(curve2, col = 'black',lwd=3,lty=2) 
if (length(curve)> 2){
  for(i in 1:(dim(curve2)[1]-1)){
    segments(curve2[,1][i],curve2[,2][i],curve2[,1][i+1],curve2[,2][i+1],col=ifelse(curve2[i,1]%in%sig,"red","black"),lwd = 3,lty=3)
  }
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top",legend = c("Species #1","Species #2","Significant"),lty=c(1,2),lwd=2, bty = "n",horiz = TRUE, col = c("black","black","red"))




open3d()
show2d({
  par(mar=c(-1,1,-1,1))
  plot(curve, axes=TRUE, type="l", xlab="Filtration Steps", ylab="Topological Summary Statistics",bty="n",lwd=3,xlim = c(-0.5,0.5),ylim=c(-0.5,1.25))
})


#1d
fruit_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = fruit_1,rate_vals = fruit_veg_new$Rate2[,2],
                                            len = len,cuts = 300,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)
veg_heat =  reconstruct_vertices_on_shape(dir = dirs,complex = veg_1,rate_vals = fruit_veg_new$Rate2[,2],
                                          len = len,cuts = 300,cone_size = directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)

fruit_colors_heat=colfunc(1 + max(fruit_heat[,1]) - min(fruit_heat[,1]))[1 + fruit_heat[,1] - min(fruit_heat[,1])]
veg_colors_heat=colfunc(1 + max(veg_heat[,1]) - min(veg_heat[,1]))[1 + veg_heat[,1] - min(veg_heat[,1])]

mfrow3d(2,1)
fruit2 = fruit1
fruit2$vb[1:3,] = fruit2$vb[1:3,]*10
fruit2$material = list(color = rep('white',dim(fruit2$vb)[2]))
plot3d(fruit2)
vcgStlWrite(fruit2,filename = 'tooth_prototype_stl')
g = col2rgb(fruit_colors_heat)
g = g/256
write.csv(t(g), col.names = FALSE, row.names = FALSE, file = 'colors.csv')
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
