###########################################################################################
### File:    null results n30.R
### Authors: Miranda Kroehl, miranda.kroehl@ucdenver.edu
###          
### Date:    April 2018
###
### Purpose: To summarize simulation results, including new supremum test
###
### Change log:
###########################################################################################

setwd("C:/Miranda/Dissertation, BW/Part1/Permutation, Single X/Programs with Bootstrap/For Github")
setwd("")

source("results helper_functions bootstrap.R")

null.results <- NULL

# X2 confounder b/w X1 and Y
load("sim_n_30_corX1X2_0.6_corX1M_0_corX2M_0_corX1Y_0.6_corX2Y_0.6_corMY_0.6.RData")
sim <- SIMULATION.RESULTS
sim.results <- summarize.allsims(sim, sims.per.seed)
null.results <- summarize.sims2(sim.results, sim.seeds, null.results)
