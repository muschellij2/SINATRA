set.seed(55)
library(sinatra)
library(FNN)
library(rgl)
library(Rvcg)
library(plyr)

#Parameters for the Analysis

cap_radius = 0.25
num_cones = 5
directions_per_cone = 5
len = 75
dirs = generate_equidistributed_cones(num_directions = num_cones, cap_radius =  cap_radius, directions_per_cone = directions_per_cone)
rotation_matrix=matrix(c(0.8065218,0.5911149,0.01028626,0,0.5583186,-0.7558216,0.34207344,0,0.1944301,-0.2816325,-0.93961692,0,0,0,0,1),ncol=4,byrow=TRUE)
ball = TRUE
ball_radius = 0.5
ec_type = 'ECT'

num_causal_region = 1
num_shared_region = 9
causal_points = 10
shared_points = 10

subdivision = 3

nsim = 25

cusps = 2*num_causal_region + num_shared_region + 1
causal_dirs = generate_equidistributed_points(cusps, cusps +1)
causal_regions_1 = sample(1:cusps,num_causal_region)
causal_regions_2 = sample((1:cusps)[-causal_regions_1],num_causal_region)
shared_regions = sample((1:cusps)[-c(causal_regions_1,causal_regions_2)],num_shared_region)


data = generate_data_sphere_simulation(nsim = nsim,dir = dirs, curve_length = len,noise_points = shared_points,
                                       causal_points = causal_points,ball_radius = ball_radius, subdivision = subdivision,
                                       cusps = cusps, causal_regions_1 = causal_regions_1, causal_regions_2 = causal_regions_2,
                                       shared_regions = shared_regions, ec_type = ec_type)


#### Sanity Checks ####

j = cor(t(data$data[,-1]))
heatmap(j)



#### Comp 2 ####
indices = c(261,295,833,224,939,278,252,593,177,357,4845,4215,4619,162,244,4618,4225,31,4705,4524,219,4301,4898,262,4860,69,64,4478,4452,4224,102,179,4176)

color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)

col_pal2=c(color1,color1,color2,color2,color3)
colfunc2 <- colorRampPalette(col_pal2)

