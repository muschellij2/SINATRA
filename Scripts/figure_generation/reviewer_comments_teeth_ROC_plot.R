##############################
### Template for new baseline plots ###
##############################

peaks = 5

roc_curve_frame = read.csv(sprintf('~/projects/Research/SINATRA/Scripts/Data/%speaks_caricature_roc.csv',peaks))
head(roc_curve_frame)


colnames(roc_curve_frame) <- c('V1','V2','V3')

### Load in the group lasso results

# Add in Group Lasso to these results
group_lasso_results_500 = read.csv(
  sprintf("~/projects/Research/SINATRA/Results/LandmarkSimulations/%speaks_1500landmarks_caricature_roc.csv",
          peaks))

colnames(group_lasso_results_500) <- c("V1","V2")
group_lasso_results_500$V3 <- "Group Lasso (1500 Landmarks)"
roc_curve_frame <- rbind(roc_curve_frame, group_lasso_results_500)


group_lasso_results_2000 = read.csv(
  sprintf("~/projects/Research/SINATRA/Results/LandmarkSimulations/%speaks_2000landmarks_caricature_roc.csv",
          peaks))
colnames(group_lasso_results_2000) <- c("V1","V2")
group_lasso_results_2000$V3 <- "Group Lasso (2000 Landmarks)"
roc_curve_frame <- rbind(roc_curve_frame, group_lasso_results_2000)


### GGplot2 code

plt <- ggplot(roc_curve_frame, aes(x = V1,y = V2,group = V3)) + 
  geom_line(data = subset(roc_curve_frame, V3 == "SINATRA"),
            aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "Limit"),
            aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5, linetype = 4) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "EN"),
            aes(x = V1,y = V2,group = V3,color = factor(V3)),
            alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "Group Lasso"),
            aes(x = V1,y = V2,group = V3,color = factor(V3)),
            alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Method") +
  ggtitle(sprintf("%s Peaks Caricatured Teeth Results", peaks)) +
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
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  scale_colour_hue(l=40)
print(plt)

ggsave(sprintf("~/Dropbox/Sub-Image Analysis/Manuscript/bioRxiv/New_Figures/caricature_%speaks_group_lasso.pdf",peaks))

##############################
### Template for old plots ###
##############################

### In Bruce's file



############################################
### Template for varying parameter plots ###
############################################


### In Bruce's file


