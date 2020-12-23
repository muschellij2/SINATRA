# SINATRA 
Sub-Image Analysis using Topological Summary Statistics.

## Purpose 

The purpose of this branch is to provide the exact scripts/code used to produce to the simulation and real data results in the SINATRA manuscript.

## Scripts
The relevant locations for the scripts used to generate each figure, and process the data for each figure are located under Scripts/figure_generation. To run the Limit Shapes algorithm, we refer to (Huang et. al 2019). http://www.lix.polytechnique.fr/~maks/papers/limit_shapes_SGP19.pdf . We will provide a link to the gitrepo/driver script when available. 

## Raw Data
For the locations of the data used in the studies, please use this link: https://www.dropbox.com/sh/rs8pjmhrwcdcuxk/AAC3Fj2_RNZLTVR_XhN4jiGxa?dl=0

## Data Generation
To generate the caricatured data, run the MATLAB scripts in Scripts/Data_Generation/GHdist/, originally sourced from : https://github.com/shaharkov/GPLmkBDMatch. The scripts to generate the data in the simulations are in Scripts/Figxx . 

## Simulations

## Further additions

This branch contains further inference functions for implementing the SINATRA pipeline, including the use of Deep Gaussian Process Classification as an inference tool. The model and inference algorithm are obtained from *Havasi et al. 2018, Inference in Deep Gaussian Processes using Stochastic Gradient Hamiltonian Monte Carlo.* Example cases using this method can be found in the script `SINATRA_deepgp_ROC.R`, which relies on the python file `sinatra_rate.py` to generate posterior samples from the Deep GP. These are then used to generate feature importance values via RATE. Code for inference in the Deep GP model is contained in `Simulations\sghmc_dgp`.

## Relevant Citations

B. Wang*, T. Sudijono*, H. Kirveslahti*, T. Gao, D.M. Boyer, S. Mukherjee, and L. Crawford. SINATRA: a sub-image analysis pipeline for selecting features that differentiate classes of 3D shapes. _bioRxiv_. 701391.
