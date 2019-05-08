################################################################################
### file:    simulator file.R
### authors: Miranda Kroehl, miranda.kroehl@ucdenver.edu
### date:    April 2018
###
### purpose: run simulations with bootstraps
###
### change log:
### 04/12/18 - file created. adapted from simulatorPerms files for dissertation work
################################################################################

rm(list = ls())

# Load required packages
library(parallel)
library(mediation)
library(boot)
library(dplyr)
library(snow)


# Set your working directory
setwd("")


# As long as all the files are in the above directory, this should call up all the functions
source("helper_functions.R")
source("mediationTestPermMethods.R")


###################################################################################
# Simulator
# size:          N (sample size)
# n.perms:       number of permutations (10k recommended)
# n.boots:       number of bootstrap redraws (5k recommended)
# sims.per.seed: number of simulations to run 
# cor.x1.x2:	 specified correlation between x1 and x2
# cor.x1.m:	 specified correlation between x1 and m
# cor.x2.m:	 specified correlation between x2 and m
# cor.x1.y:	 specified correlation between x1 and y
# cor.x2.y:	 covariance between x2 and y
# cor.m.y:	 specified correlation between m and y
###################################################################################

size <- 30
n.perms <- 10000
n.boots <- 5000
sims.per.seed <- 1000

cor.x1.x2 <- 0.6
cor.x1.m  <- 0 
cor.x2.m  <- 0
cor.x1.y  <- 0.6
cor.x2.y  <- 0.6
cor.m.y   <- 0.6

# Set up the initial seed values for a simulation run:
sim.seeds <- expand.grid(size, n.perms, n.boots, cor.x1.x2, cor.x1.m, cor.x2.m,
                         cor.x1.y, cor.x2.y, cor.m.y)

sim.seeds <- data.frame(sim.seeds)

options(cl.cores = detectCores() - 1)

# set up a cluster, 
this.cluster <- makeCluster(getOption("cl.cores", 2))
clusterCall(cl = this.cluster, fun = function(){library(lasso2); library(boot)})
clusterExport(cl = this.cluster, list = objects())
SIMULATION.RESULTS <- 
  parLapply(cl = this.cluster,
            1:nrow(sim.seeds),
            function(idx) {
              args <- as.list(sim.seeds[idx, ])
              # args$conf.levels = conf.levels
              formals(mediationTestPermMethods) <- args
              rtn <- 
                replicate(sims.per.seed, mediationTestPermMethods(), simplify = FALSE)
              return(rtn)
            })
stopCluster(this.cluster)

### Save the simulation results - Peter set it up to save as an image, which can then
# be loaded (below) and summarized using the helper functions
save.image(file = paste("sim_n_", size, "_corX1X2_", cor.x1.x2, "_corX1M_", cor.x1.m,
                        "_corX2M_", cor.x2.m, "_corX1Y_", cor.x1.y, "_corX2Y_", cor.x2.y,
                        "_corMY_", cor.m.y, ".RData", sep = ""))



