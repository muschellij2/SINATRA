### Clear Console ###
cat("\014")

### Clear Environment ###
rm(list = ls(all = TRUE))

### Load SINATRA functions ###
source("Final Code for Functions (Really)/RATEv2.R")
source("Final Code for Functions (Really)/sinatra/R/GPC_Approximate_Inference.R")
source("Final Code for Functions (Really)/Simulation_Functions.R")
source("Final Code for Functions (Really)/SINATRA_Simulation_Balled_EC copy.R")
sourceCpp("Final Code for Functions (Really)/sinatra/R/BAKRGibbs.cpp")

### Load in the R Libraries ###
library(truncnorm)
library(doParallel)

### Set the parameters for the analysis ###
set.seed(4913, kind = "L'Ecuyer-CMRG")
n.simulations <- 100
reconstruction_type <- "vertex" #set how we assess each reconstruction
causal_points <- 15
shared_points <- 5
num_cones <- 15

######################################################################################
######################################################################################
######################################################################################

### Setup DoParallel ###
no_cores <- detectCores() - 1 
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl)  

### Run the analysis in Parallel ###
simulation_results <- foreach(i=1:n.simulations, .combine = 'rbind', .noexport = c('GaussKernel')) %dopar% {
  set.seed(i)
  results <- generate_metric_curve(nsim = 100, curve_length = 20, grid_size = 25,
                        distance_to_causal_point = 0.1, causal_points = causal_points,
                        shared_points = shared_points, num_cones = num_cones,
                        directions_per_cone = 4, eta = 0.1, ball_radius = 2.5,type = reconstruction_type)
   
  rdf <- data.frame(cbind(results,rep(i,dim(results)[1])) )
  rdf <- cbind(rdf,1:num_cones)
  rdf <- plyr::rename(rdf,c("X1" = "TPR", "X2" = "Intersection" , "X3" = "Trial", "1:num_cones" = "Cone"))
  
}

stopCluster(cl)

######################################################################################
######################################################################################
######################################################################################


# Recast results into dataframe
mdf <- melt(simulation_results, id.vars = c('Trial','Cone'),variable.name = 'metrics')

######################################################################################
######################################################################################
######################################################################################

### Average the results together over trials ###

mdfmeans <- aggregate(mdf["value"], by = list(Cone = mdf$Cone, Metric = mdf$metric), mean)
mdfstd <- aggregate(mdf["value"], by = list(Cone = mdf$Cone, Metric = mdf$metric), sd)
mdfmeans <- cbind(mdfmeans,mdfstd$value)
mdfmeans <- plyr::rename(mdfmeans, c("value" = "Mean","mdfstd$value" = "Std"))
limits <- aes(ymax = mdfmeans$Mean + mdfmeans$Std, 
              ymin = mdfmeans$Mean - mdfmeans$Std)

### Plot the averaged results ###
avg_curve_plt <- ggplot(data <- mdfmeans,aes(x = Cone, y = Mean, fill = Metric)) +
  geom_bar(stat = "identity",position = position_dodge(0.9)) +
  geom_errorbar(limits, position = position_dodge(0.9),
                width = 0.25) + 
  labs(x = "Number of Cones", y = "Metric Value") +
  ggtitle(sprintf("causal:%d,shared:%d,type:%s",causal_points,shared_points,reconstruction_type)) +
  scale_fill_discrete(name = "Metric Type") +
  scale_x_continuous(breaks = 1:max(mdfmeans$Cone))

######################################################################################
######################################################################################
######################################################################################

### Save the results ###

# save the raw results
sim_results_file = "~/Desktop/Simulations/Results/data_metric_curve_causal15_shared5_vertex.RData"
save(simulation_results, file = sim_results_file)

# save the final plot
avg_curve_file = "~/Desktop/Simulations/Results/plot_metric_curve_causal15_shared5_vertex.pdf"
pdf(file = avg_curve_file)
plot(avg_curve_plt)
dev.off()


