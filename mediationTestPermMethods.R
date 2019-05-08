################################################################################
### file:    mediation.testPermMethods.R
### authors: Miranda Kroehl, miranda.kroehl@ucdenver.edu
### date:    31 May 2013
###
### purpose: Evaluate tests for mediation using different permutation methods
###          
###
### change log:
### 31 May 2013 - file created.  original code from newMediationTest function
### March/Apr 2018 - edit file to incorporate bootstraps
################################################################################

#################################################################################################
#  mediation.testPermMethods: evaluated 3 permutation methods for testing mediation
#
#  Adjusted Model: Y = gamma0 + gamma1*x1 + gamma2*x2 + gamma3*m + ey
#  IE Model:	     M = alpha0 + alpha1*x1 + alpha2*x2
#
#  Arguments:
#  size:        sample size
#  n.perms      number of times permutation is performed 
#  n.boots      number of bootstrap samples created to construct CIs
#  mu.x1, mu.x2, mu.em, mu.ey: means for random variable generation
#  sd.x1, sd.x2, sd.em, sd.ey: sd for random variable generation
#  cor.x1.x2, cor.x1.m, cor.x2.m: Correlations for data generation
#  cor.x1.y, cor.x2.y, cor.m.y:   Correlations for data generaion
#
# Returns:
#   a list (rtn)
#################################################################################################

mediationTestPermMethods <- 
  function(size, n.perms, n.boots,
           cor.x1.x2, cor.x1.m, cor.x2.m,
           cor.x1.y, cor.x2.y, cor.m.y){

    rtn <-  list(tt	  	       =  NULL,
		            means          =  NULL,
                 sds           =  NULL,
                 cormat        =  NULL,
 
                 beta.ols      =  NULL,
                 gamma.ols     =  NULL,
                 alpha.ols     =  NULL,
                 se.beta.ols   =  NULL,
                 se.gamma.ols  =  NULL,
                 se.alpha.ols  =  NULL,

                 tx1.orig.ols  =  NULL,

        		     p.tx1.null 	 	  = NULL,
        		     p.tx1.null.Mstar	= NULL,
        		     p.alpha.null 	  = NULL,
        		     p.gamma.null 	  = NULL,
        		     ci.tx1star.alt.ols	= NULL,
        		     ci.tx1star.alt.ols.Mstar	= NULL,
        
        		     accept.ho.tx1   	   = NULL,
        		     accept.ho.tx1.Mstar = NULL,
        		     accept.ho.joint 	   = NULL,
        		     accept.ho.ci	    	 = NULL,
        		     accept.ho.ci.Mstar	 = NULL,
        		     ci.covered 	   	   = NULL,
        		     ci.covered.Mstar  	 = NULL,
		            
		             bootci.perc = NULL,
		             bootci.bc   = NULL,
		             bootci.bca  = NULL,
		             accept.ho.boot.perc = NULL,
		             accept.ho.boot.bc   = NULL,
		             accept.ho.boot.bca  = NULL
    )

  simdat <- rdata(size,
                  cor.x1.x2, cor.x1.m, cor.x2.m,
                  cor.x1.y, cor.x2.y, cor.m.y) 

  ############################################################
  # true regression coefficient values
  tt <- true.Coef(cor.x1.x2, cor.x1.m, cor.x2.m,
           cor.x1.y, cor.x2.y, cor.m.y)

  rtn$tt <- tt

  ############################################################

  ############################################################
  # summary statistics about the generated data set
  rtn$means  <- colMeans(simdat)
  rtn$sds    <- apply(simdat, 2, sd)
  rtn$cormat <- cor(simdat)

  ############################################################
  # model fits
  crude.ols <- lm(y ~ x1 + x2, data = simdat, model = FALSE)
  adj.ols   <- update(crude.ols, formula = . ~ . + m)
  ie.ols    <- update(crude.ols, formula = m ~ .)

  # extract coeffiencts
  rtn$beta.ols  <- coef(crude.ols)
  rtn$gamma.ols <- coef(adj.ols)
  rtn$alpha.ols <- coef(ie.ols)

  # standard errors
  rtn$se.beta.ols <- sqrt(diag(vcov(crude.ols)))
  rtn$se.gamma.ols   <- sqrt(diag(vcov(adj.ols)))
  rtn$se.alpha.ols    <- sqrt(diag(vcov(ie.ols)))

  ############################################################
  # Data needs to be in matrix form for the following functions
  my.matrix <- as.matrix(simdat)
  y 	     <- as.matrix(my.matrix[,1])
  adj.mod0 <- as.matrix(my.matrix[, 2:3])
  adj.mod  <- as.matrix(my.matrix[, 2:4])
  m 	     <- as.matrix(my.matrix[,4])
  ie.mod0  <- as.matrix(my.matrix[,3])
  ie.mod   <- as.matrix(my.matrix[,2:3])

  # Set design matrices
  adj.design 		<- cbind(1, adj.mod)
  adj.reduced.design 	<- cbind(1, adj.mod0)
  ie.design 		<- cbind(1, ie.mod)
  ie.reduced.design 	<- cbind(1, ie.mod0)

  # Obtain fitted outcomes and residuals
  adj.fit    <- lm.fit(adj.design, y)
  adj.fitted <- t(adj.fit$fitted.values)
  adj.resid  <- t(adj.fit$residuals)

  adj.fit0 		   <- lm.fit(adj.reduced.design, y)
  adj.reduced.fitted <- t(adj.fit0$fitted.values)
  adj.reduced.resid  <- t(adj.fit0$residuals)

  ie.fit    <- lm.fit(ie.design, m)
  ie.fitted <- t(ie.fit$fitted.values)
  ie.resid  <- t(ie.fit$residuals)

  ie.fit0 		  <- lm.fit(ie.reduced.design, m)
  ie.reduced.fitted <- t(ie.fit0$fitted.values)
  ie.reduced.resid  <- t(ie.fit0$residuals)

  ############################################################
  # Calculate test statistic
  adj.stat 		<- setdiff(colnames(adj.mod), colnames(adj.mod0))
  gamma.orig 	<- adj.fit$coef[adj.stat]

# ie.stat 		<- setdiff(colnames(ie.mod), colnames(ie.mod0))
  ie.stat		<- "x1"
  alpha.orig 	<- ie.fit$coef[ie.stat]

  tx1.orig.ols 	 <- alpha.orig * gamma.orig
  rtn$tx1.orig.ols <- alpha.orig * gamma.orig


  ############################################################
  # Permutations 
  n.greater.tx1		<- 0
  n.greater.tx1.Mstar 	<- 0
  n.greater.alpha 	<- 0
  n.greater.gamma 	<- 0
  n.greater.gamma.Mstar	<- 0
 
  li <- freedman_lane_sim(adj.design, adj.reduced.design, ie.design, ie.reduced.design,
			 	    adj.fitted, adj.reduced.fitted, adj.resid, adj.reduced.resid,
			 	    ie.fitted, ie.reduced.fitted, ie.resid, ie.reduced.resid,
				    n.greater.tx1, n.greater.tx1.Mstar, n.greater.alpha, n.greater.gamma, 
				    n.perms, adj.stat, ie.stat, 
				    gamma.orig, alpha.orig, tx1.orig.ols) 

  p.tx1.null  	     <- (1 + li$n.greater.tx1) / (1 + n.perms)
  p.tx1.null.Mstar   <- (1 + li$n.greater.tx1.Mstar) / (1 + n.perms)
  p.alpha.null 	     <- (1 + li$n.greater.alpha) / (1 + n.perms)
  p.gamma.null 	     <- (1 + li$n.greater.gamma) / (1 + n.perms)
  ci.tx1star.alt.ols <- li$ci.tx1star.alt.ols
  ci.tx1star.alt.ols.Mstar <- li$ci.tx1star.alt.ols.Mstar

  rtn$p.tx1.null		   <- p.tx1.null
  rtn$p.tx1.null.Mstar <- p.tx1.null.Mstar
  rtn$p.alpha.null	   <- p.alpha.null
  rtn$p.gamma.null 	   <- p.gamma.null
  rtn$ci.tx1star.alt.ols <- ci.tx1star.alt.ols
  rtn$ci.tx1star.alt.ols.Mstar <- ci.tx1star.alt.ols.Mstar

  rtn$accept.ho.tx1   		<- as.numeric(p.tx1.null > 0.05)
  rtn$accept.ho.tx1.Mstar	<- as.numeric(p.tx1.null.Mstar > 0.05)
  rtn$accept.ho.joint 		<- as.numeric(p.alpha.null > 0.05 | p.gamma.null > 0.05) 
  rtn$accept.ho.ci    		<- as.numeric(coverZero(rbind(li$ci.tx1star.alt.ols)))
  rtn$accept.ho.ci.Mstar 	<- as.numeric(coverZero(rbind(li$ci.tx1star.alt.ols.Mstar)))

  # Evaluate Coverage for ho.ci
  t.alpha1  <- tt[1]
  t.gamma3  <- tt[5]
  t.product <- t.alpha1*t.gamma3

  rtn$ci.covered <- as.numeric(t.product >= ci.tx1star.alt.ols[1] && t.product <= ci.tx1star.alt.ols[2])
  rtn$ci.covered.Mstar <- as.numeric(t.product >= ci.tx1star.alt.ols.Mstar[1] && t.product <= ci.tx1star.alt.ols.Mstar[2])

  ############################################################
  # Bootstrap
  ############################################################
  
  ie.fcn <- function(simdat)
  {
    # Data needs to be in matrix form    
    my.matrix <- as.matrix(simdat)
    y 	     <- as.matrix(my.matrix[,1])
    m 	     <- as.matrix(my.matrix[,4])
    ie.mod   <- as.matrix(my.matrix[,2:3])
    adj.mod  <- as.matrix(my.matrix[, 2:4])
    
    # Set design matrices
    adj.design 		<- cbind(1, adj.mod)
    ie.design 		<- cbind(1, ie.mod)
    
    # Calculate regression coefficients
    adj.fit    <- lm.fit(adj.design, y)
    ie.fit    <- lm.fit(ie.design, m)
    
    # Calculate test statistic
    adj.stat 		<- setdiff(colnames(adj.mod), colnames(adj.mod0))
    gamma.orig 	<- adj.fit$coef[adj.stat]
    
    # ie.stat 		<- setdiff(colnames(ie.mod), colnames(ie.mod0))
    ie.stat		<- "x1"
    alpha.orig 	<- ie.fit$coef[ie.stat]
    
    tx1.orig.ols 	 <- alpha.orig * gamma.orig
    return(tx1.orig.ols)
  }
  
  #create function that calculates indirect effect for a given set of data
  medi<-function(idf,i,...)
  {
    dt<-idf[i,]
    bsmdl<-ie.fcn(dt)
  }
  
  #Run bootstrap (Data, function, Num of iterations, ordinary type, using indicies)
  boot_qty <- n.boots
  bs<-boot(simdat,medi,R=boot_qty,sim="ordinary", stype="i") 
  
  # Pull off CI for percentile, bias corrected, and bias corrected accelerated
  perc <- boot.ci(bs, type="perc")
  bootci.perc <- perc$percent[4:5]
  bc <- boot.ci(bs, type="bca", L=1)
  bootci.bc <- bc$bca[4:5]
  bca <- boot.ci(bs, type="bca")
  bootci.bca <- bca$bca[4:5]
  
  rtn$bootci.perc <- bootci.perc
  rtn$bootci.bc <- bootci.bc
  rtn$bootci.bca <- bootci.bca
  
  rtn$accept.ho.boot.perc <- as.numeric(coverZero(rbind(bootci.perc)))
  rtn$accept.ho.boot.bc   <- as.numeric(coverZero(rbind(bootci.bc)))
  rtn$accept.ho.boot.bca  <- as.numeric(coverZero(rbind(bootci.bca)))
  
  ############################################################
  
  return(rtn) 

}

# Example run with default settings
# system.time({
# testrslt <-  mediation.testPermMethods(size=100, n.perms=1000, n.boot=1000,
#           cor.x1.x2=0, cor.x1.m=0, cor.x2.m=0,
#           cor.x1.y=0.3, cor.x2.y=0.3, cor.m.y=0.3)
# })


###################
### End of File ###
###################
