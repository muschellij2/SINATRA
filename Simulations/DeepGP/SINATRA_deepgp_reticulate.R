library(sinatra)
library(Rvcg)
library(rgl)
library(ggplot2)

set.seed(17)

### set parameters
nsim = 25
curve_length = 25
distance_to_causal_point = 0.1
causal_points = 3
shared_points = 6
num_cones = 25
truncated = 300
two_curves = TRUE
ball = TRUE
ball_radius = 2.5
type = 'vertex'
directions_per_cone = 10
cap_radius = 0.15
radius = 0
ec_type = 'ECT'
mode = 'sphere'

subdivision = 3
num_causal_region = 5
num_shared_region = 5
grid_size = 25



### generate sphere data
cusps = 50
causal_dirs = generate_equidistributed_points(cusps,cusps)
causal_regions_1 = sample(1:cusps,num_causal_region)
causal_regions_2 = sample((1:cusps)[-causal_regions_1],num_causal_region)
shared_regions = sample((1:cusps)[-c(causal_regions_1,causal_regions_2)],num_shared_region)
directions <- generate_equidistributed_cones(num_cones,cap_radius,directions_per_cone)
data <- generate_data_sphere_simulation(nsim = nsim,dir = directions, curve_length = curve_length,noise_points = shared_points,
                                        causal_points = causal_points, ball_radius = ball_radius, subdivision = subdivision,
                                        cusps = cusps, causal_regions_1 = causal_regions_1, causal_regions_2 = causal_regions_2,
                                        shared_regions = shared_regions, ec_type = ec_type)
ec_curve_data <- data$data
num_cones <- dim(directions)[1]/directions_per_cone


### Do the inference


# let's use reticulate
library(reticulate)
# change EC curve data to proper form
write.matrix(ec_curve_data, file = "Results/sinatra_ECT_dataROC.txt", sep = ",")

# run python script
conda_install("pandas")
py_run_file("sghmc_dgp/sinatra_rate.py")

# read in posterior samples
# read in data matrix
posterior_samples <- read.table('./posteriorsamples_sinatraROC.txt', header = FALSE)
posterior_samples <- as.matrix(posterior_samples)

X <-  read.table('./designmatrix_sinatraROC.txt', header = FALSE)
X <- as.matrix(X)

### get the rate_values
RATE <- RATE(X, posterior_samples)
rate_values <- RATE$RATE 

# or get rate values using GPC model:
# standardize ec_curve_data
rate_values_GPC <- find_rate_variables_with_other_sampling_methods(ec_curve_data, bandwidth = 0.01, type = 'ESS')[,2]

#standardized data
cleaned_data <- ec_curve_data
cleaned_data <- ec_curve_data - matrix(rep(apply(ec_curve_data,2,mean),dim(ec_curve_data)[1]), nrow = dim(ec_curve_data)[1], byrow = TRUE)
dev <- apply(cleaned_data,2,sd)
dev[dev == 0] = 1
cleaned_data <- t(t(cleaned_data) / dev)

rate_values_cleaned_GPC <- find_rate_variables_with_other_sampling_methods(cleaned_data, bandwidth = 0.01, type = 'EP')[,2]


#Indices for Two Classes
index1 = seq(1,nsim,2)
complex_points1 = data$complex_points[index1]

index2 = seq(2,nsim,2)
complex_points2 = data$complex_points[index2]




roc_curve1 =  compute_roc_curve_vertex(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                       curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                       eta = eta, directions_per_cone = directions_per_cone, directions = directions, class = 1,truncated = truncated,
                                       ball_radius = ball_radius, radius = radius, mode = mode,subdivision = subdivision)
roc_curve1 = cbind(roc_curve1, rep(1,dim(roc_curve1)[1]))
roc_curve1 = cbind(roc_curve1,(1:dim(roc_curve1)[1]))

roc_curve2 =  compute_roc_curve_vertex(data = data, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                       curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                       eta = eta, directions_per_cone = directions_per_cone, directions = directions,class = 2,truncated = truncated,
                                       ball_radius = ball_radius, radius = radius, mode = mode, subdivision = subdivision)
roc_curve2 = cbind(roc_curve2, rep(2,dim(roc_curve2)[1]))
roc_curve2 = cbind(roc_curve2,(1:dim(roc_curve2)[1]))

roc_curve = rbind(roc_curve1,roc_curve2)

roc_curve <- data.frame(roc_curve1)
plot <- ggplot(d <- roc_curve, aes(x = X1,y = X2)) +
  ggtitle(sprintf("DGPC: Medium Simulation Size 5")) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)") +
  geom_line(stat = 'Identity') +
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
print(plot)












