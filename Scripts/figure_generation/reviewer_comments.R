set.seed(4913, kind = "L'Ecuyer-CMRG")
library(sinatra)
library(FNN)
library(rgl)
library(Rvcg)
library(plyr)
library(pdist)
library(gglasso)
library(numbers)
library(data.table)

#Parameters for the Analysis

cap_radius = 0.15
num_cones = 5
directions_per_cone = 5
len = 75

num_causal_region = 2
num_shared_region = 1
causal_points = 10
shared_points = 10

subdivision = 3

nsim = 25

cusps = 50
causal_dirs = generate_equidistributed_points(cusps, cusps +1)
causal_regions_1 = sample(1:cusps,num_causal_region)
causal_regions_2 = sample((1:cusps)[-causal_regions_1],num_causal_region)
shared_regions = sample((1:cusps)[-c(causal_regions_1,causal_regions_2)],num_shared_region)

#### Fig 2.a. ####
g1 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                   mode = 'sphere_baseline',num_cusps = cusps,
                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region)


write.csv(g1,'~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared.csv')

roc_curve_frame = data.frame(g1)
library(ggplot2)
roc_curve_frame = roc_curve_frame[roc_curve_frame[,3] == 1,]
roc_curve_frame[,4] = rep('Baseline',dim(roc_curve_frame)[2])
library(ggplot2)
ggplot(roc_curve_frame, aes(x = X1,y = X2,group = X4)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X4) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Scenario') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "ROC Curve  2 Causal Region, 1 Shared Region")
ggsave('~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared.pdf')

#### Fig 2.b. ####
num_causal_region = 6
num_shared_region = 3
causal_points = 10
shared_points = 10
g2 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                   mode = 'sphere_baseline',num_cusps = cusps,
                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region)


write.csv(g2,'~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared.csv')

roc_curve_frame = data.frame(g2)
library(ggplot2)
roc_curve_frame = roc_curve_frame[roc_curve_frame[,3] == 1,]
roc_curve_frame[,4] = rep('Baseline',dim(roc_curve_frame)[2])
library(ggplot2)
ggplot(roc_curve_frame, aes(x = X1,y = X2,group = X4)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X4) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Scenario') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "ROC Curve  6 Causal Region, 3 Shared Region")
ggsave('~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared.pdf')
#### Fig 2.c. ####
num_causal_region = 10
num_shared_region = 5
causal_points = 10
shared_points = 10
g3 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                   mode = 'sphere_baseline',num_cusps = cusps,
                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region)


write.csv(g3,'~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared.csv')

roc_curve_frame = data.frame(g3)
roc_curve_frame = roc_curve_frame[roc_curve_frame[,3] == 1,]
roc_curve_frame[,4] = rep('Baseline',dim(roc_curve_frame)[2])
library(ggplot2)
ggplot(roc_curve_frame, aes(x = X1,y = X2,group = X4)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(X4) )) +
  geom_abline(slope = 1,intercept = 0) +
  scale_x_continuous(name='False Positive Rate',limits=c(0,1)) +
  scale_y_continuous(name='True Positive Rate', limits=c(0,1)) +
  labs(color='Scenario') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size=10)) + ggtitle(label = "ROC Curve  10 Causal Region, 5 Shared Region")
ggsave('~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared.pdf')

#### Test New Function ####
library(glmnet)
data2 = generate_data_sphere_simulation_new(nsim = nsim,dir = dirs, curve_length = len,noise_points = shared_points,
                                       causal_points = causal_points,ball_radius = ball_radius, subdivision = subdivision,
                                       cusps = cusps, causal_regions_1 = causal_regions_1, causal_regions_2 = causal_regions_2,
                                       shared_regions = shared_regions, ec_type = ec_type)
groups = rep(1:(dim(data2$data)[2]/3),each = 3)
#lasso = cv.gglasso(x = data2$data[,-1], y = data2$data[,1], group = groups, loss = "logit")
lasso = cv.glmnet(x = data2$data[,-1], y = data2$data[,1], alpha = 0.0,  family = "binomial")
coefs  = coef(lasso, s = 'lambda.min')

#Binary Option
nonzero = which(abs(coefs[-1])>0.8)
results = unique(groups[nonzero])


# Heatmap Option


g = cbind(groups, abs(coefs[-1]))

library(dplyr)
df = as.data.table(g)

new_df = aggregate(df[,2],list(df$groups),mean)

abs_coefs = new_df$V2
results = groups[which(abs_coefs > 0.5)]
cuts = cut(abs_coefs,50,labels = FALSE)


color1='blue'
color2='lightgreen'
color3='red'
#col_pal = c("#00007F", "blue", "#007FFF", "cyan",  "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
col_pal=c(color1,color1,color2,color2,color2,color3)
col_pal=c(color1,color2,color3)

colfunc <- colorRampPalette(col_pal)
heat_colors=colfunc(max(cuts) - min(cuts))[cuts - min(cuts)]

mesh = vcgSphere(subdivision = 3)
mesh$vb[1:3,] = t(data2$complex_points[[1]])
cols = rep('white', dim(mesh$vb)[2])
cols[data$causal_points1] = 'red'
cols[data$causal_points2] = 'red'
cols[data$noise] = 'blue'
cols[results] = 'green'

mfrow3d(1,2)
plot3d(mesh, col = cols)
plot3d(mesh, col = heat_colors)
#### Sanity Checks ####

j = cor(t(data_summary$data[,-1]))
heatmap(j)


generate_ROC_with_coned_directions(nsim = 10, curve_length = 50, grid_size = 25, distance_to_causal_point = 0.1, 
                                   causal_points = causal_points,shared_points = shared_points,  eta = 0.1, 
                                   truncated = -1, two_curves = TRUE, ball = TRUE, ball_radius = 2.5, type = 'vertex',
                                   min_points = 3,directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                   mode = 'sphere_baseline', 
                                   subdivision = 3,num_causal_region = 2, num_shared_region = 1)


#### Comp 2 ####
indices = c(261,295,833,224,939,278,252,593,177,357,4845,4215,4619,162,244,4618,4225,31,4705,4524,219,4301,4898,262,4860,69,64,4478,4452,4224,102,179,4176)

color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)

col_pal2=c(color1,color1,color2,color2,color3)
colfunc2 <- colorRampPalette(col_pal2)

