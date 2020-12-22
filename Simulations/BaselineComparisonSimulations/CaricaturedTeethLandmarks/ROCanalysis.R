

# use Euclidean farthest point sampling to generate the landmarks on the caricatured teeth

# grab the coordinates of the landmarks, run lasso / elastic net on these

# also grab the points that are causal for the shape.

# visualize the ROC

# build a training / test set framework?

### Update mesh functions to generate data. -- in the function real data summary, create_comparison_matrix

library(sinatra)
library(Rvcg)
set.seed(15)

#Parameters for the Analysis
num_landmarks <- 500

### Helper functions
get_euclidean_fps_landmarks = function(mesh, num_landmarks){
  ### Gets landmarks using farthest point sampling
  landmarks = rdist::farthest_point_sampling(t(mesh$vb)[,1:3], metric = 'euclidean', num_landmarks)
  landmarks
}

original_tooth_path = "/Users/timothysudijono/Dropbox/Data + experiments/Data/all_files_realignedv3/Aligned_Shapes/clean_V13_sas.off"
mesh = Rvcg::vcgImport(original_tooth_path)

landmark_indices = get_euclidean_fps_landmarks(mesh, num_landmarks)

### inspect / plot the teeth in the directories..

base_dir = "~/Dropbox/CaricaturedTeeth/new_aligned_shapesv4/"
data_dirs = list.dirs(base_dir,recursive = FALSE)

roc_curves <- list()

for (i in 1:20) {
#for (i in 1:5) {
  #i=6 example of NAs - need to debug this.
  cat("on directory", i)
  dir = data_dirs[i]
  old_data_dir = paste(dir,'/mesh/gp1',sep='')
  new_data_dir = paste(dir,'/mesh/gp2',sep='')
  
  old_data_files = list.files(old_data_dir, full.names = TRUE)
  new_data_files = list.files(new_data_dir, full.names = TRUE)
  
  class_1_probs = read.csv(paste(dir,'/gp1_spt.csv',sep=''), header = FALSE) # lets us get the teue vertices, > 0.25 = gamma
  class_2_probs = read.csv(paste(dir,'/gp2_spt.csv',sep=''), header = FALSE)

  # varying the number of directions, but want to do something totally different
  #total_dirs = c(1)
  
  landmark_pset = list(base_shape_dir = "~/Dropbox/Data + experiments/Data/all_files_realignedv3/Aligned_Shapes/clean_V13_sas.off",
                       num_landmarks = num_landmarks)
  
  data_summary = real_data_summary(shape_transformation = 'landmark', dir1=new_data_dir,dir2 = old_data_dir,
                                   base_shape_dir = landmark_pset[[1]], num_landmarks = landmark_pset[[2]], 
                                   mode = "elastic_net", alpha = 0.25) 
  
  ### These aren't necessarily rate_values they're just importance values -- we can replace  them with the lasso coefficients? How do these even serve as importance values? 
  ### Seems like just a jeuristic
  
  roc_curve = compute_roc_curve_teeth(data_dir1 = old_data_dir, data_dir2 = new_data_dir, 
                                      gamma = 0.25,class_1_probs = class_1_probs,class_2_probs = class_2_probs,
                                      rate_values = data_summary$Rate2, 
                                      directions_per_cone = NULL, curve_length = NULL, directions = NULL,
                                      mode = "landmark", base_shape_dir = landmark_pset[[1]])
  #print(roc_curve)

  roc_frame = data.frame(roc_curve)
  roc_curves[[i]] = roc_frame
  
  ### visualize / debug
  
  # which(abs(data_summary$inference_result)[,2] > 0)
  # 
  # base_shape = Rvcg::vcgImport(landmark_pset[[1]])
  # landmark_indices = get_euclidean_fps_landmarks(base_shape, landmark_pset[[2]])
  # 
  # gamma = 0.25
  # class1 = intersect(which(class_1_probs > gamma), landmark_indices)
  # class2 = intersect(which(class_2_probs > gamma), landmark_indices)
  # true_vertices = union(class1, class2)
  # 
  # chosen_vertices = which(data_summary$Rate2 > 0)
  # 
  # mesh = vcgImport(new_data_files[2])
  # rgl::plot3d(mesh,color= 'white')
  # rgl::points3d(t(mesh$vb)[true_vertices,1:3], color = "#2a6bbe", size = 5)
  # rgl::points3d(t(mesh$vb)[chosen_vertices,1:3], color = '#ff4d4d', size = 10)
  # 
  # rgl::plot3d(mesh,color= 'white')
  # for(k in 5:6){
  #   rgl::points3d(t(vcgImport(new_data_files[k])$vb)[1:5,1:3], color = '#ff4d4d', size = 10)
  # }
  
}

# handle NaNs -- just omit them for now
total_roc = matrix(0, nrow = num_landmarks, ncol = 2)
not_null = 0

for (j in 1:length(roc_curves)){
  if (!is.nan(roc_curves[[j]][1,2])){
    total_roc = total_roc + roc_curves[[j]][,1:2]
    not_null  = not_null+1
  }
}
total_roc = total_roc/not_null

### Save and plot ###

#write.csv(total_roc,
          file = '~/projects/Research/SINATRA/Results/LandmarkSimulations/5peaks_500landmarks_caricature_roc.csv',
          row.names = FALSE)

library(ggplot2)
ggplot() + 
  geom_line(data = total_roc, aes(x = X1,y = X2),alpha = 0.75,  size = 1.5, color = 'red') +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)") +
  ggtitle(sprintf("5 Caricatured Peaks: Group Lasso + Landmark")) +
  coord_cartesian(xlim= c(0,1.0))+
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold")) +
  scale_colour_hue(l=40)
ggsave('~/projects/Research/SINATRA/Results/LandmarkSimulations/5peaks_500landmarks_caricature_roc.pdf')


## Ideas about why Lasso doesn't work  -- there;s much more noise on the vertex level; this is removed using  the ECT transform
## Lasso doesn't work as a feature selection method as well as RATE in the presence of correlated noise.

 ## Throw out the curves with NA? what does this mean?

### Debugging tools:
#visualize reconstructed points
#plot the vectorizations, they should be different

# another thing: when debugging the model, have to debug the group lasso cross validated coefficients.
# what is margin based loss versus imsclassification error?


# double check if the teeth are aligned. 
# Seemingly, NaNs randomly appear. Does CV introduce some sort of noise?