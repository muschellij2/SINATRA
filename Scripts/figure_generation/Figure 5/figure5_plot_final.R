##############################
### Template for new baseline plots ###
##############################

peaks = 3

roc_curve_frame = read.csv(sprintf('~/projects/Research/SINATRA/Scripts/Data/%speaks_caricature_roc.csv',peaks))
head(roc_curve_frame)


colnames(roc_curve_frame) <- c('V1','V2','V3')

### Load in the group lasso results

# Add in Group Lasso to these results
group_lasso_results_500 = read.csv(
  sprintf("~/projects/Research/SINATRA/Results/LandmarkSimulations/%speaks_500landmarks_caricature_roc.csv",
          peaks))

colnames(group_lasso_results_500) <- c("V1","V2")
group_lasso_results_500$V3 <- "Group Lasso (500 Landmarks)"
roc_curve_frame <- rbind(roc_curve_frame, group_lasso_results_500)


group_lasso_results_2000 = read.csv(
  sprintf("~/projects/Research/SINATRA/Results/LandmarkSimulations/%speaks_2000landmarks_caricature_roc.csv",
          peaks))
colnames(group_lasso_results_2000) <- c("V1","V2")
group_lasso_results_2000$V3 <- "Group Lasso (2000 Landmarks)"
roc_curve_frame <- rbind(roc_curve_frame, group_lasso_results_2000)


### GGplot2 code
library(ggplot2)
library(scales)
library(data.table)
show_col(hue_pal()(7))
hue_pal()(7)

plt <- ggplot(roc_curve_frame, aes(x = V1,y = V2,group = V3)) + 
  geom_line(data = subset(roc_curve_frame, V3 == "SINATRA"),
            aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5) +
  geom_line(data = subset(roc_curve_frame, V3 == "Group Lasso (500 Landmarks)"),
            aes(x = V1,y = V2,group = V3, color = factor(V3)),
            alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
  geom_line(data = subset(roc_curve_frame, V3 == "Group Lasso (2000 Landmarks)"), 
            aes(x = V1,y = V2,group = V3, color = factor(V3)),
            alpha = 0.75,  size = 1.5, linetype = 4, position=position_jitter(w=0.01, h=0)) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "Limit"),
            aes(x = V1,y = V2,group = V3,color = factor(V3)),alpha = 0.75,  size = 1.5, linetype = 4) +
  geom_line(data = subset(roc_curve_frame, V3 %like% "EN"),
            aes(x = V1,y = V2,group = V3,color = factor(V3)),
            alpha = 0.75,  size = 1.5, linetype = 4,position=position_jitter(w=0.01, h=0)) +
  scale_color_manual(values = c("SINATRA" = "Black",
                                "Group Lasso (500 Landmarks)" = "#C49A00",
                                "Group Lasso (2000 Landmarks)" = "#53B400",
                                "Limit Shapes" = "#F8766D",
                                "Limit Shapes (Misspecified)" = "#00B6EB",
                                "Baseline (EN Max)" = "#A58AFF",
                                "Baseline (EN Mean)" = "#FB61D7")) + 
  
  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", color = "Method") +
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
  guides(color = guide_legend(override.aes = list(size = 1.5))) 
  #scale_color_viridis(discrete = TRUE) 
  #scale_color_brewer(palette="Accent") 
  #scale_colour_hue(l=50, c = 100)
  
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


