# library(grplasso)
library(gglasso) # might be easier to use
library(glmnet)
library(sinatra)
library(pdist)
library(scatterplot3d)
library(ggplot2)

library(rgl)
library(Rvcg)
library(FNN)
#source('~/projects/Research/SINATRA/Simulations/BaselineComparisonSimulations/GroupLassoSimulation/GenerateSimulatedShapes.R')



#Parameters for the Analysis
num_causal_region = 5
num_shared_region = 10
causal_region_size = 20
shared_region_size = 20

subdivision = 4
nsim = 25

num_shape_scaffold = 100
landmarks =1000
# what is cusps? number of equidistributed points around the sphere.
#plot3d(causal_dirs[,1],causal_dirs[,2],causal_dirs[,3])
causal_regions_1 = sample(1:num_shape_scaffold,num_causal_region)
causal_regions_2 = sample((1:num_shape_scaffold)[-causal_regions_1],num_causal_region)
shared_regions = sample((1:num_shape_scaffold)[-c(causal_regions_1,causal_regions_2)],num_shared_region)


#### Load in simulated shapes

# There's a function to generate these shapes!

shape_data <- generate_landmark_data_sphere_simulation(nsim = nsim, noise_region_size = shared_region_size,
                                            causal_region_size = causal_region_size, subdivision = subdivision,
                                            num_shape_scaffold = num_shape_scaffold, num_landmarks = landmarks,
                                            causal_regions_1 = causal_regions_1, causal_regions_2 = causal_regions_2,
                                            shared_regions = shared_regions)

data = shape_data[[1]]
class1_positive_landmarks = shape_data[[2]]
class2_positive_landmarks = shape_data[[3]]
# which regions are causal though?

#### Don't encode EC -- instead vectorize via the landmarks on the shape


#### Apply Group Lasso using the grplasso package to the data, defining groups as the coordinates of the landmark points.
# Do both group lasso and elastic net with group parameters.
groups <- rep(1:(dim(data)[2]/3),each = 3)
#lasso <- cv.gglasso(x = data[,-1], y = data[,1], group = groups, loss = "logit")

lasso = cv.glmnet(x = data[,-1], y = data[,1], alpha = 1.0,  family = "binomial",
                  intercept = TRUE, loss = "logit", nlambda = 500, pred.loss = "loss",nfolds = 10) # how do these work?
coefs <- coef(lasso, s = 'lambda.min')
# pick a cross validation parameter.
#coefs <- coef(lasso, s = c(0,0.0001,0.001,0.01,0.10,0.5,1.0,10))
coefs

#Binary Option
# normalize the coefficients. Before sweeping.
nonzero = which(abs(coefs)>0) # why is this 0.8?
results = unique(groups[nonzero])
results

causal_regions_1
causal_regions_2


#landmarks = matrix(data[1,-1],ncol=3,byrow = TRUE)
#rgl::plot3d(landmarks[,1], landmarks[,2], landmarks[,3])
#need to scatter plot with RGL.

## Can also run GPs

#### Vary the regularization parameter, see which groups the algorithm picks out.
# sumamrize by picking the mean / median


#### Can compare these ROC curves directly to the sinatra simulation ...
causal_regions = 1
shared_regions = 2
generate_ROC_baseline(num_shape_scaffold = 50,
                      num_landmarks = 500,
                      num_causal_region = causal_regions,
                      num_shared_region = shared_regions)


# is the issuewith LASSO and highly correlated variables? Because landmarks are in close proximity, the variables are
# hence quite correlated. How does Lasso select  between two correlated variables?

# vary the number of cusps...
# how do you make the TPR / FPR fair?
causal_regions = 1
shared_regions = 2
averaged_roc_curve <- generate_averaged_ROC_baseline(runs = 10,
                                                     num_shape_scaffold = 25,
                                                     num_landmarks = 500,
                                                     num_causal_region = causal_regions,
                                                     num_shared_region = shared_regions)
data <- data.frame(averaged_roc_curve)

### Save the data

#write.csv(data,
          file = '~/projects/Research/SINATRA/Results/LandmarkSimulations/group_lasso_1causal_2shared_2000landmarks_sphere_simulated_roc.csv',
          row.names = FALSE)

#### Summarize plots....
ROC_curve_plt <- ggplot(data <- data, aes(x = X1, y = X2)) +
  geom_line(stat = "identity") +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)",
       title = sprintf("Landmark Representation ROC curve, %s Causal, %s Shared", causal_regions, shared_regions )) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt)

ggsave('~/projects/Research/SINATRA/Results/LandmarkSimulations/group_lasso_1causal_2shared_2000landmarks_sphere_simulated_roc.pdf')



