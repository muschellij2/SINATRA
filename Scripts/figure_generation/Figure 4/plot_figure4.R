##############################
### Template for new baseline plots ###
##############################
library(ggplot2)
library(data.table)

### load in previous results
causal_points = 5
shared_points = 10

roc_curve_frame = read.csv(sprintf('~/projects/Research/SINATRA/Scripts/Data/sphere_roc_%scausal_%sshared2.csv',
                                   shared_points,
                                   causal_points))
head(roc_curve_frame)

# Add in Group Lasso to these results
group_lasso_results_500 = read.csv(sprintf('%sgroup_lasso_%scausal_%sshared_500landmarks_sphere_simulated_roc.csv',
                                "~/projects/Research/SINATRA/Results/LandmarkSimulations/",
                                causal_points, shared_points))
colnames(group_lasso_results_500) <- c("V1","V2")
group_lasso_results_500$V3 <- "Group Lasso (500 Landmarks)"
roc_curve_frame <- rbind(roc_curve_frame, group_lasso_results_500)


group_lasso_results_2000 = read.csv(sprintf('%sgroup_lasso_%scausal_%sshared_2000landmarks_sphere_simulated_roc.csv',
                                           "~/projects/Research/SINATRA/Results/LandmarkSimulations/",
                                           causal_points, shared_points))
colnames(group_lasso_results_2000) <- c("V1","V2")
group_lasso_results_2000$V3 <- "Group Lasso (2000 Landmarks)"
roc_curve_frame <- rbind(roc_curve_frame, group_lasso_results_2000)


### ggplot code
ROC_curve_plot <- ggplot(roc_curve_frame, aes(x = V1,y = V2,group = V3)) + 
  geom_line(data = subset(roc_curve_frame, V3 == "SINATRA"),
            aes(x = V1,y = V2,group = V3, color = factor(V3)),alpha = 0.75,  size = 1.5) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "Limit"),
            aes(x = V1,y = V2,group = V3, color = factor(V3)),alpha = 0.75,  size = 1.5, linetype = 4) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "EN"),
            aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5, linetype = 4) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "Group Lasso"),
            aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5, linetype = 4) +
  scale_color_manual(values = c("SINATRA" = "Black",
                                "Group Lasso (500 Landmarks)" = "#C49A00",
                                "Group Lasso (2000 Landmarks)" = "#53B400",
                                "Limit Shapes" = "#F8766D",
                                "Limit Shapes (Misspecified)" = "#00B6EB",
                                "Baseline (EN Max)" = "#A58AFF",
                                "Baseline (EN Mean)" = "#FB61D7")) + 
  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", color = "Method") +
  ggtitle(sprintf("%s Causal Regions, %s Shared Regions, Size 10",causal_points,shared_points)) +
  #coord_cartesian(xlim=c(0, 0.2)) + 
  geom_abline(intercept = 0, slope = 1, alpha =0.5) + 
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=16,face="bold")) +
  guides(color = guide_legend(override.aes = list(size = 1.5))) 
  #scale_colour_hue(l=40)
print(ROC_curve_plot)


ggsave(sprintf("~/Dropbox/Sub-Image Analysis/Manuscript/bioRxiv/New_Figures/SphereSim_%scausal_%sshared_baseline_group_lasso.pdf",
               causal_points, shared_points))


##############################
### Template for Old plots ###
##############################

causal_points = 5
shared_points = 10
reconstruction_type = 'EC'
load(sprintf("~/projects/Research/SINATRA/Simulations/Sphere_Simulation/new/df_ROC_causal%s_shared%s_10.RData",causal_points,shared_points))

### Plot results ###
# plot the first class for simplicity
class_one_ROC <- rdfmeans[rdfmeans$Class == 1,]
ROC_curve_plt <- ggplot(data <- class_one_ROC, aes(x = FPR, y = TPR, color = Num_Cones)) +
  geom_line(stat = "identity",size = 0.6) +
  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", color = "# Cones") +
  ggtitle(sprintf("%s Causal Regions, %s Shared Regions, Size 10", causal_points, shared_points)) +
  coord_cartesian(xlim=c(0, 1.0)) + 
  geom_abline(intercept = 0, slope = 1, alpha =0.5) + 
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=16,face="bold")) +
  guides(color = guide_legend(override.aes = list(size = 1.5)))
  #scale_colour_hue(l=40)
print(ROC_curve_plt)

ggsave(
  sprintf("~/Dropbox/Sub-Image Analysis/Manuscript/bioRxiv/New_Figures/Sphere_Sims/causal%s_shared%s_10.pdf",
          causal_points,shared_points))

##########################
### Varying parameters ###
##########################
causal_points = 5
shared_points = 10
#can't recall difference between underscore 5 and 10 ... 

load(sprintf("~/projects/Research/SINATRA/Simulations/Sphere_Simulation/df_ROC_coneangle_causal%s_shared%s_10.RData",
             causal_points, shared_points))

### Plot results ###
class_one_ROC <- rdfmeans[rdfmeans$Class == 1,]
ROC_curve_plt <- ggplot(data <- class_one_ROC, aes(x = FPR, y = TPR, color = Radii)) +
  geom_line(stat = "identity", size = 0.6) +
  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", color = "Cone Angle") +
  ggtitle(sprintf("Varying Cone Angle")) +
  geom_abline(intercept = 0, slope = 1) + 
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=16,face="bold")) +
  guides(color = guide_legend(override.aes = list(size = 1.5))) 
  #scale_colour_hue(l=40)
print(ROC_curve_plt)
ggsave(
  sprintf("~/Dropbox/Sub-Image Analysis/Manuscript/bioRxiv/New_Figures/Sphere_Sims/Varying_Params/ConeAngle/causal%s_shared%s_10.pdf",
          causal_points,shared_points))


load(sprintf("~/projects/Research/SINATRA/Simulations/Sphere_Simulation/df_ROC_curvelength_causal%s_shared%s_10.RData",
             causal_points, shared_points))

### Plot results ###
class_one_ROC <- rdfmeans[rdfmeans$Class == 1,]
ROC_curve_plt <- ggplot(data <- class_one_ROC, aes(x = FPR, y = TPR, color = Lengths)) +
  geom_line(stat = "identity",size = 0.6) +
  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", color = "Lengths") +
  ggtitle(sprintf("Varying Curve Length")) +
  geom_abline(intercept = 0, slope = 1) + 
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=16,face="bold")) +
  guides(color = guide_legend(override.aes = list(size = 1.5)))  
  #scale_colour_hue(l=40)
print(ROC_curve_plt)
ggsave(
  sprintf("~/Dropbox/Sub-Image Analysis/Manuscript/bioRxiv/New_Figures/Sphere_Sims/Varying_Params/CurveLength/causal%s_shared%s_10.pdf",
          causal_points,shared_points))


load(sprintf("~/projects/Research/SINATRA/Simulations/Sphere_Simulation/df_ROC_dirpercone_causal%s_shared%s_10.RData",
             causal_points, shared_points))


### Plot results ###
class_one_ROC <- rdfmeans[rdfmeans$Class == 1,]
ROC_curve_plt <- ggplot(data <- class_one_ROC, aes(x = FPR, y = TPR, color = Dir_Per_Cone)) +
  geom_line(stat = "identity", size = 0.6) +
  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", color = "# Directions") +
  ggtitle(sprintf("Varying Directions Per Cone")) +
  geom_abline(intercept = 0, slope = 1)+
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=16,face="bold")) +
  guides(color = guide_legend(override.aes = list(size = 1.5))) 
  #scale_colour_hue(l=40)
print(ROC_curve_plt)
ggsave(
  sprintf("~/Dropbox/Sub-Image Analysis/Manuscript/bioRxiv/New_Figures/Sphere_Sims/Varying_Params/DirectionsPerCone/causal%s_shared%s_10.pdf",
          causal_points,shared_points))

#######################################
### Reformat raw simulation results ###

load("~/projects/Research/SINATRA/Simulations/Sphere_Simulation/new/data_ROC_curve_length_causal5_shared10_region_size10.RData")

rdfmeans <- aggregate(simulation_results[c("FPR","TPR")],
                      by = list("Lengths" = simulation_results$Lengths,
                                "Index" = simulation_results$Index,
                                "Class" = simulation_results$Class), mean)

rdfmeans$Lengths <- as.factor(rdfmeans$Lengths)

df_results_file = sprintf("~/projects/Research/SINATRA/Simulations/Sphere_Simulation/df_ROC_curvelength_causal%d_shared%d_%s.RData",
                          5,10,10)
save(rdfmeans, file = df_results_file)



