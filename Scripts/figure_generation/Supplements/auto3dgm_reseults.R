#### Aligning and Scaling the meshes ####
library(Morpho)
library(rgl)
library(Rvcg)
dir='Data/doug_new_teeth_scaled'
out_path = 'Data/doug_new_teeth_scaled_aligned'

b = list.files(dir, full.names = TRUE)

Data_dir = dir
Output_dir = out_path
rotation_matrix=matrix(c(0.8065218,0.5911149,0.01028626,0,0.5583186,-0.7558216,0.34207344,0,0.1944301,-0.2816325,-0.93961692,0,0,0,0,1),ncol=4,byrow=TRUE)


open3d()
rgl.bg(color = 'white')
#rgl.bg(color = 'pink')
mfrow3d(nr=6,nc = 10)
rgl.bg(color = 'white')
tarsius_files = list.files('~/Documents/doug_new_teeth_raw_by_species/Tarsius/',full.names = TRUE)
microcebus_files = list.files('~/Documents/doug_new_teeth_raw_by_species/Microcebus/',full.names = TRUE)
mirza_files = list.files('~/Documents/doug_new_teeth_raw_by_species/Mirza/',full.names = TRUE)
saimiri_files = list.files('~/Documents/doug_new_teeth_raw_by_species/Saimiri/',full.names = TRUE)

for (file_name in tarsius_files){
  file = vcgImport(file_name)
  plot3d(file, color = 'white', axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
  rgl.bg(color = 'lightblue')
  rgl.viewpoint(userMatrix =rotation_matrix)
}
for (file_name in microcebus_files){
  file = vcgImport(file_name)
  plot3d(file, color = 'white', axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
  rgl.bg(color = 'rosybrown1')
  rgl.viewpoint(userMatrix =rotation_matrix)
}

for (file_name in mirza_files){
  file = vcgImport(file_name)
  plot3d(file, color = 'white', axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
  rgl.bg(color = 'lightgreen')
  rgl.viewpoint(userMatrix =rotation_matrix)
}
for (file_name in saimiri_files){
  file = vcgImport(file_name)
  plot3d(file, color = 'white', axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
  rgl.bg(color = 'gold')
  rgl.viewpoint(userMatrix =rotation_matrix)
}

open3d()
rgl.bg(color = 'white')
#rgl.bg(color = 'pink')
mfrow3d(nr=6,nc = 10)
rgl.bg(color = 'white')


tarsius_files = list.files('~/Dropbox (Princeton)/Data + Experiments Tim Sudijono/Data/doug_new_teeth_by_species/Tarsius/',full.names = TRUE)
microcebus_files = list.files('~/Dropbox (Princeton)/Data + Experiments Tim Sudijono/Data/doug_new_teeth_by_species/Microcebus/',full.names = TRUE)
mirza_files = list.files('~/Dropbox (Princeton)/Data + Experiments Tim Sudijono/Data/doug_new_teeth_by_species/Mirza/',full.names = TRUE)
saimiri_files = list.files('~/Dropbox (Princeton)/Data + Experiments Tim Sudijono/Data/doug_new_teeth_by_species/Saimiri/',full.names = TRUE)

for (file_name in tarsius_files){
  file = vcgImport(file_name)
  plot3d(file, color = 'white', axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
  rgl.bg(color = 'blue')
  rgl.viewpoint(userMatrix =rotation_matrix)
}
for (file_name in microcebus_files){
  file = vcgImport(file_name)
  plot3d(file, color = 'white', axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
  rgl.bg(color = 'red')
  rgl.viewpoint(userMatrix =rotation_matrix)
}

for (file_name in mirza_files){
  file = vcgImport(file_name)
  plot3d(file, color = 'white', axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
  rgl.bg(color = 'green')
  rgl.viewpoint(userMatrix =rotation_matrix)
}
for (file_name in saimiri_files){
  file = vcgImport(file_name)
  plot3d(file, color = 'white', axes = FALSE, xlab = '',ylab = '',zlab='',specular = 'black')
  rgl.bg(color = 'orange')
  rgl.viewpoint(userMatrix =rotation_matrix)
}
