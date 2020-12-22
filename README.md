# SINATRA 
Sub-Image Analysis using Topological Summary Statistics.

## Purpose 

The purpose of this branch is to provide the exact scripts/code used to produce to the simulation and real data results in the SINATRA manuscript.

## Scripts
The relevant locations for the scripts used to generate each figure, and process the data for each figure are located under Scripts/figure_generation. 

## Simulations

## Further additions

This branch contains further inference functions for implementing the SINATRA pipeline, including the use of Deep Gaussian Process Classification as an inference tool. The model and inference algorithm are obtained from *Havasi et al. 2018, Inference in Deep Gaussian Processes using Stochastic Gradient Hamiltonian Monte Carlo.* Example cases using this method can be found in the script `SINATRA_deepgp_ROC.R`, which relies on the python file `sinatra_rate.py` to generate posterior samples from the Deep GP. These are then used to generate feature importance values via RATE. Code for inference in the Deep GP model is contained in `Simulations\sghmc_dgp`.

## Relevant Citations

B. Wang*, T. Sudijono*, H. Kirveslahti*, T. Gao, D.M. Boyer, S. Mukherjee, and L. Crawford. SINATRA: a sub-image analysis pipeline for selecting features that differentiate classes of 3D shapes. _bioRxiv_. 701391.
