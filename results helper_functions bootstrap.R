###########################################################################################
### File:    Summarizing Perm Results NEW.R
### Authors: Miranda Kroehl, miranda.kroehl@ucdenver.edu
###          
### Date:    April 2018
###
### Purpose: To summarize simulation results, including new supremum test
###
### The functions contained in this script are:
###
### summarize.Results - generates datasets for simulations
###
### Change log:
###########################################################################################


###########################################################################################
###########################################################################################

summarize.allsims <- function(sim, sims.per.seed){
  sim.results <- NULL

  for (i in 1:sims.per.seed){
    mean.y 	<- sim[[1]][[i]]$mean[1]
    mean.x1 	<- sim[[1]][[i]]$mean[2]
    mean.x2 	<- sim[[1]][[i]]$mean[3]
    mean.m 	<- sim[[1]][[i]]$mean[4]
    sd.y	 	<- sim[[1]][[i]]$sds[1]
    sd.x1 	<- sim[[1]][[i]]$sds[2]
    sd.x2 	<- sim[[1]][[i]]$sds[3]
    sd.m 		<- sim[[1]][[i]]$sds[4]
    
    cor.x1.x2	<- sim[[1]][[i]]$cormat[2,3]
    cor.x1.m	<- sim[[1]][[i]]$cormat[2,4]
    cor.x2.m	<- sim[[1]][[i]]$cormat[3,4]
    cor.x1.y	<- sim[[1]][[i]]$cormat[1,2]
    cor.x2.y	<- sim[[1]][[i]]$cormat[1,3]
    cor.m.y	<- sim[[1]][[i]]$cormat[1,4]
    
    
    t.alpha1 	<- sim[[1]][[i]]$tt[1,1]
    t.alpha2 	<- sim[[1]][[i]]$tt[1,2]
    t.gamma1 	<- sim[[1]][[i]]$tt[1,3]
    t.gamma2 	<- sim[[1]][[i]]$tt[1,4]
    t.gamma3 	<- sim[[1]][[i]]$tt[1,5]
    t.beta1 	<- sim[[1]][[i]]$tt[1,6]
    t.beta2 	<- sim[[1]][[i]]$tt[1,7]

    beta1.ols	  <- sim[[1]][[i]]$beta.ols[2]
    beta2.ols	  <- sim[[1]][[i]]$beta.ols[3]
    gamma1.ols	  <- sim[[1]][[i]]$gamma.ols[2]
    gamma2.ols	  <- sim[[1]][[i]]$gamma.ols[3]
    gamma3.ols	  <- sim[[1]][[i]]$gamma.ols[4]
    alpha1.ols 	  <- sim[[1]][[i]]$alpha.ols[2]
    alpha2.ols 	  <- sim[[1]][[i]]$alpha.ols[3]
    se.beta1.ols  <- sim[[1]][[i]]$se.beta.ols[2]
    se.beta2.ols  <- sim[[1]][[i]]$se.beta.ols[3]
    se.alpha1.ols <- sim[[1]][[i]]$se.alpha.ols[2]
    se.alpha2.ols <- sim[[1]][[i]]$se.alpha.ols[3]
    se.gamma1.ols <- sim[[1]][[i]]$se.gamma.ols[2]
    se.gamma2.ols <- sim[[1]][[i]]$se.gamma.ols[3]
    se.gamma3.ols <- sim[[1]][[i]]$se.gamma.ols[4]
    
    t.new.ols		      <- sim[[1]][[i]]$tx1.orig.ols
    p.tx1.null  	    <- sim[[1]][[i]]$p.tx1.null
    p.alpha.null  	  <- sim[[1]][[i]]$p.alpha.null
    p.gamma.null  	  <- sim[[1]][[i]]$p.gamma.null
    p.supremum.null   <- max(p.tx1.null, p.alpha.null, p.gamma.null)
    lci.tx1star.alt.ols <- sim[[1]][[1]]$ci.tx1star.alt.ols[1] 
    uci.tx1star.alt.ols <- sim[[1]][[1]]$ci.tx1star.alt.ols[2] 
    
    accept.ho.tx1 	<- sim[[1]][[i]]$accept.ho.tx1
    accept.ho.joint 	<- sim[[1]][[i]]$accept.ho.joint
    accept.ho.supremum <- as.numeric(p.supremum.null > 0.05) 
    accept.ho.ci 	<- sim[[1]][[i]]$accept.ho.ci
    
    lci.perc <- sim[[1]][[i]]$bootci.perc[1] 
    uci.perc <- sim[[1]][[i]]$bootci.perc[2] 
    lci.bc <- sim[[1]][[i]]$bootci.bc[1] 
    uci.bc <- sim[[1]][[i]]$bootci.bc[2] 
    lci.bca <- sim[[1]][[i]]$bootci.bca[1] 
    uci.bca <- sim[[1]][[i]]$bootci.bca[2] 
    
    accept.ho.boot.perc <- sim[[1]][[i]]$accept.ho.boot.perc
    accept.ho.boot.bc   <- sim[[1]][[i]]$accept.ho.boot.bc
    accept.ho.boot.bca  <- sim[[1]][[i]]$accept.ho.boot.bca
    
    
    sim.results <- rbind(sim.results, 
                         cbind(i, mean.y, mean.x1, mean.x2, mean.m,
                               sd.y, sd.x1, sd.x2, sd.m,
                               cor.x1.x2, cor.x1.m, cor.x2.m, 
                               cor.x1.y, cor.x2.y, cor.m.y,
                               
                               t.alpha1, t.alpha2, t.gamma1, t.gamma2, t.gamma3, t.beta1, t.beta2,
                               
                               beta1.ols, beta2.ols, alpha1.ols, alpha2.ols,
                               gamma1.ols, gamma2.ols, gamma3.ols,
                               se.beta1.ols, se.beta2.ols, se.alpha1.ols, se.alpha2.ols, 
                               se.gamma1.ols, se.gamma2.ols, se.gamma3.ols,  
                               
                               t.new.ols, p.tx1.null, p.alpha.null, p.gamma.null, p.supremum.null,
                               lci.tx1star.alt.ols, uci.tx1star.alt.ols,
                               lci.perc, uci.perc, lci.bc, uci.bc, lci.bca, uci.bca,
                               accept.ho.tx1, accept.ho.joint, accept.ho.supremum, accept.ho.ci, 
                               accept.ho.boot.perc, accept.ho.boot.bc, accept.ho.boot.bca 
                         ))
    
  } ####### end for loop #######
  
  return(sim.results)

} ####### end summarize.allsims #######
  

###########################################################################################
###########################################################################################

summarize.sims2 <- function(sim.results, sim.seeds, null.results){

  sim.results2 <- data.frame(sim.results)
  
  n     <- sim.seeds$Var1
  x1.x2	<- sim.seeds$Var10
  x1.m 	<- sim.seeds$Var11
  x2.m 	<- sim.seeds$Var12
  x1.y 	<- sim.seeds$Var13
  x2.y 	<- sim.seeds$Var14
  m.y 	<- sim.seeds$Var15
  
  # Get Means 
  avg.cor.x1.x2 	<- mean(sim.results2$cor.x1.x2)
  avg.cor.x1.m 	<- mean(sim.results2$cor.x1.m)
  avg.cor.x2.m	<- mean(sim.results2$cor.x2.m)
  avg.cor.x1.y 	<- mean(sim.results2$cor.x1.y)
  avg.cor.x2.y 	<- mean(sim.results2$cor.x2.y)
  avg.cor.m.y 	<- mean(sim.results2$cor.m.y)
  
  avg.alpha1.ols	<- mean(sim.results2$alpha1.ols)
  avg.alpha2.ols	<- mean(sim.results2$alpha2.ols)
  avg.gamma3.ols	<- mean(sim.results2$gamma3.ols)
  min.alpha1.ols	<- min(sim.results2$alpha1.ols)
  min.alpha2.ols	<- min(sim.results2$alpha2.ols)
  min.gamma3.ols	<- min(sim.results2$gamma3.ols)
  max.alpha1.ols	<- max(sim.results2$alpha1.ols)
  max.alpha2.ols	<- max(sim.results2$alpha2.ols)
  max.gamma3.ols	<- max(sim.results2$gamma3.ols)
  
  avg.t.new.ols 	  <- mean(sim.results2$t.new.ols)
  avg.p.tx1.null    <- mean(sim.results2$p.tx1.null)
  avg.p.alpha.null  <- mean(sim.results2$p.alpha.null)
  avg.p.gamma.null  <- mean(sim.results2$p.gamma.null)
  avg.p.supremum.null <- mean(sim.results2$p.supremum.null)
  
  avg.lci.tx1star.alt.ols <- mean(sim.results2$lci.tx1star.alt.ols)
  avg.uci.tx1star.alt.ols <- mean(sim.results2$uci.tx1star.alt.ols)
  
  avg.lci.perc <- mean(sim.results2$lci.perc)
  avg.uci.perc <- mean(sim.results2$uci.perc)
  avg.lci.bc <- mean(sim.results2$lci.bc)
  avg.uci.bc <- mean(sim.results2$uci.bc)
  avg.lci.bca <- mean(sim.results2$lci.bca)
  avg.uci.bca <- mean(sim.results2$uci.bca)
  
  type1.ho.tx1  	<- 1 - mean(sim.results2$accept.ho.tx1)
  type1.ho.joint 	<- 1 - mean(sim.results2$accept.ho.joint )
  type1.ho.supremum <- 1 - mean(sim.results2$accept.ho.supremum )
  type1.ho.ci  	<- 1 - mean(sim.results2$accept.ho.ci)
  
  type1.ho.perc  	<- 1 - mean(sim.results2$accept.ho.boot.perc)
  type1.ho.bc   	<- 1 - mean(sim.results2$accept.ho.boot.bc)
  type1.ho.bca  	<- 1 - mean(sim.results2$accept.ho.boot.bca)
  
  null.results <- rbind(null.results, cbind(n, x1.x2, x1.m, x2.m, x1.y, x2.y, m.y,
                                            avg.cor.x1.x2, avg.cor.x1.m, avg.cor.x2.m, 
                                            avg.cor.x1.y, avg.cor.x2.y, avg.cor.m.y,
                                            avg.alpha1.ols, avg.alpha2.ols, avg.gamma3.ols, 
                                            min.alpha1.ols, min.alpha2.ols, min.gamma3.ols, 
                                            max.alpha1.ols, max.alpha2.ols, max.gamma3.ols, 			 
                                            
                                            avg.t.new.ols, avg.p.tx1.null, avg.p.alpha.null, 
                                            avg.p.gamma.null, avg.p.supremum.null,
                                            avg.lci.tx1star.alt.ols, avg.uci.tx1star.alt.ols,
                                            avg.lci.perc, avg.uci.perc, avg.lci.bc, avg.uci.bc, avg.lci.bca, avg.uci.bca,
                                            type1.ho.tx1, type1.ho.joint, type1.ho.supremum, type1.ho.ci,
                                            type1.ho.perc, type1.ho.bc, type1.ho.bca
  ))
  
  return(null.results)
  
}  ####### end summarize.sims2 ######
  
  
###########################################################################################
###########################################################################################


# This function summarizes the alt sims
summarize.alts <- function(sim.results, sim.seeds, alt.results){
  
  sim.results2 <- data.frame(sim.results)
  
  n     <- sim.seeds$Var1
  x1.x2	<- sim.seeds$Var10
  x1.m 	<- sim.seeds$Var11
  x2.m 	<- sim.seeds$Var12
  x1.y 	<- sim.seeds$Var13
  x2.y 	<- sim.seeds$Var14
  m.y 	<- sim.seeds$Var15
  
  # Get Means 
  avg.cor.x1.x2 	<- mean(sim.results2$cor.x1.x2)
  avg.cor.x1.m 	<- mean(sim.results2$cor.x1.m)
  avg.cor.x2.m	<- mean(sim.results2$cor.x2.m)
  avg.cor.x1.y 	<- mean(sim.results2$cor.x1.y)
  avg.cor.x2.y 	<- mean(sim.results2$cor.x2.y)
  avg.cor.m.y 	<- mean(sim.results2$cor.m.y)
  
  avg.alpha1.ols	<- mean(sim.results2$alpha1.ols)
  avg.alpha2.ols	<- mean(sim.results2$alpha2.ols)
  avg.gamma3.ols	<- mean(sim.results2$gamma3.ols)
  min.alpha1.ols	<- min(sim.results2$alpha1.ols)
  min.alpha2.ols	<- min(sim.results2$alpha2.ols)
  min.gamma3.ols	<- min(sim.results2$gamma3.ols)
  max.alpha1.ols	<- max(sim.results2$alpha1.ols)
  max.alpha2.ols	<- max(sim.results2$alpha2.ols)
  max.gamma3.ols	<- max(sim.results2$gamma3.ols)
  
  avg.t.new.ols 	  <- mean(sim.results2$t.new.ols)
  avg.p.tx1.null    <- mean(sim.results2$p.tx1.null)
  avg.p.alpha.null  <- mean(sim.results2$p.alpha.null)
  avg.p.gamma.null  <- mean(sim.results2$p.gamma.null)
  avg.p.supremum.null <- mean(sim.results2$p.supremum.null)
  
  avg.lci.tx1star.alt.ols <- mean(sim.results2$lci.tx1star.alt.ols)
  avg.uci.tx1star.alt.ols <- mean(sim.results2$uci.tx1star.alt.ols)
  
  avg.lci.perc <- mean(sim.results2$lci.perc)
  avg.uci.perc <- mean(sim.results2$uci.perc)
  avg.lci.bc <- mean(sim.results2$lci.bc)
  avg.uci.bc <- mean(sim.results2$uci.bc)
  avg.lci.bca <- mean(sim.results2$lci.bca)
  avg.uci.bca <- mean(sim.results2$uci.bca)
  
  power.ho.tx1  	<- 1 - mean(sim.results2$accept.ho.tx1)
  power.ho.joint 	<- 1 - mean(sim.results2$accept.ho.joint )
  power.ho.supremum 	<- 1 - mean(sim.results2$accept.ho.supremum )
  power.ho.ci  	<- 1 - mean(sim.results2$accept.ho.ci)
  
  power.ho.perc  	<- 1 - mean(sim.results2$accept.ho.boot.perc)
  power.ho.bc   	<- 1 - mean(sim.results2$accept.ho.boot.bc)
  power.ho.bca  	<- 1 - mean(sim.results2$accept.ho.boot.bca)
  
  alt.results <- rbind(alt.results, cbind(n, x1.x2, x1.m, x2.m, x1.y, x2.y, m.y,
                                            avg.cor.x1.x2, avg.cor.x1.m, avg.cor.x2.m, 
                                            avg.cor.x1.y, avg.cor.x2.y, avg.cor.m.y,
                                            avg.alpha1.ols, avg.alpha2.ols, avg.gamma3.ols, 
                                            min.alpha1.ols, min.alpha2.ols, min.gamma3.ols, 
                                            max.alpha1.ols, max.alpha2.ols, max.gamma3.ols, 			 
                                            
                                            avg.t.new.ols, avg.p.tx1.null, avg.p.alpha.null, 
                                            avg.p.gamma.null, avg.p.supremum.null,
                                            avg.lci.tx1star.alt.ols, avg.uci.tx1star.alt.ols,
                                            avg.lci.perc, avg.uci.perc, avg.lci.bc, avg.uci.bc, avg.lci.bca, avg.uci.bca,
                                            power.ho.tx1, power.ho.joint, power.ho.supremum, power.ho.ci,
                                            power.ho.perc, power.ho.bc, power.ho.bca
  ))
  
  
  
  return(alt.results)
  
}  ####### end summarize.alts ######
