###########################################################################################
### File:    helper_functions.R
### Authors: Miranda Kroehl, miranda.kroehl@ucdenver.edu
###          Peter DeWitt,   peter.dewitt@ucdenver.edu
### Date:    21 Feb 2013
###
### Purpose: Helper functions to generate the statistics needed for
###          permutation testing in the mediation analysis simulations
###
### The functions contained in this script are:
###
### rdata - generates datasets for simulations
### true.Coef - calculates the true regression coefficients from the correlation structure
### lm_fit_fast - solves regresson models quickly, returning beta coefficients
### .freedman_lane_sim - permutes residuals
### coverZero - deterimes if a ci covers zero
###
### Change log:
### 21 Feb 2013 - File created.  Original code copied from Kroehl's work and
###               edited by DeWitt.
### 6 Mar 2013 - File edited.  
###		    (1) Change equations for simulating data with correlations
###		    (2) Add test statistics for lasso estimates to tstat function
###		    (3) Modify lassoEstimates function to allow for variable to be swept out
###
### 9 May 2013 - File edited.
###		    (1) Changed cover function from coverZero to coverT 
### 21 May 2013 - File edited.
###		    (1) Permutations will be used for determining type 1 and power
###			  Added functions .freedman_lane_sim, and lm_fit_fast
###		    (2) Went back to coverZero function
###
### 1 Oct 2013 - File edited.
### 		    (1) Added true.Coeff function
###########################################################################################


############################################################
  # rdata: generates random sets of correlated data
  #
  # Arguments:
  #   size:        sample size
  #   cor.x1.x2:	 specified correlation between x1 and x2
  #   cor.x1.m:	 specified correlation between x1 and m
  #   cor.x2.m:	 specified correlation between x2 and m
  #   cor.x1.y:	 specified correlation between x1 and y
  #   cor.x2.y:	 covariance between x2 and y
  #   cor.m.y:	 specified correlation between m and y
  #
  # Returns:
  #   a data.frame
############################################################

getMultivariateGaussianSample <- function(rows, correlation) {
  cols = ncol(correlation)
  cholMatrix = chol(correlation)
  randomNormalMatrix = matrix(nrow=rows,rnorm(rows*cols));
  return (randomNormalMatrix %*% cholMatrix)
}

rdata <- function(size, 
                  cor.x1.x2, cor.x1.m, cor.x2.m,
                  cor.x1.y, cor.x2.y, cor.m.y) {
 
  cor1 	  <- matrix(c(1, 		cor.x1.x2,  cor.x1.m,	cor.x1.y, 
			cor.x1.x2, 	1, 		cor.x2.m, 	cor.x2.y, 
			cor.x1.m, 	cor.x2.m, 	1,		cor.m.y,
			cor.x1.y,	cor.x2.y,	cor.m.y,	1), 		
			nrow=4, ncol=4)

  sim <- getMultivariateGaussianSample(size, cor1)

  y  <- sim[,4]
  x1 <- sim[,1]
  x2 <- sim[,2]
  m  <- sim[,3]

  return(data.frame(y, x1, x2, m))

}



############################################################
  # true.Coef: calculates the true regression coefficients
  #   to be used in the simulations, based on correlations
  #
  # Arguments:
  #   cor.x1.x2:	 specified correlation between x1 and x2
  #   cor.x1.m:	 specified correlation between x1 and m
  #   cor.x2.m:	 specified correlation between x2 and m
  #   cor.x1.y:	 specified correlation between x1 and y
  #   cor.x2.y:	 covariance between x2 and y
  #   cor.m.y:	 specified correlation between m and y
  #
  # Returns:
  #   t.alpha1: 	true alpha1
  # 	t.alpha2:	true alpha2
  #	t.gamma1:	true gamma1
  #	t.gamma2:	true gamma2
  #	t.gamma3:	true gamma3
  #	t.beta1:	true beta1
  #	t.beta2:	true beta2

############################################################

true.Coef <- function(cor.x1.x2, cor.x1.m, cor.x2.m,
                      cor.x1.y, cor.x2.y, cor.m.y) {
 
  x.ie <- matrix(c(1, cor.x1.x2, cor.x1.x2, 1), nrow=2, ncol=2) 
  y.ie <- matrix(c(cor.x1.m, cor.x2.m), nrow=2, ncol=1)
  model.ie <- solve(x.ie, y.ie)
  t.alpha1 <- model.ie[1]
  t.alpha2 <- model.ie[2]

  x.adj <- matrix(c(1, cor.x1.x2, cor.x1.m, cor.x1.x2, 1, cor.x2.m, cor.x1.m, cor.x2.m, 1), nrow=3, ncol=3)
  y.adj <- matrix(c(cor.x1.y, cor.x2.y, cor.m.y), nrow=3, ncol=1)
  model.adj <- solve(x.adj, y.adj)
  t.gamma1  <- model.adj[1]
  t.gamma2  <- model.adj[2]
  t.gamma3  <- model.adj[3]

  x.crude <- matrix(c(1, cor.x1.x2, cor.x1.x2, 1), nrow=2, ncol=2) 
  y.crude <- matrix(c(cor.x1.y, cor.x2.y), nrow=2, ncol=1)
  model.crude <- solve(x.crude , y.crude )
  t.beta1  <- model.crude[1]
  t.beta2  <- model.crude[2]

  return(cbind(t.alpha1, t.alpha2, t.gamma1, t.gamma2, t.gamma3,
		   t.beta1, t.beta2))

}



###############################################################################
  # lm_fit_fast: solves regresson models quickly, returning beta coefficients
  #
  # Arguments:
  #   dat:	outcome matrix
  #	X:	predictor matrix
  # 	out:  coefficient of interest
  #
  # Returns:
  #   beta coefficient 
###############################################################################

lm_fit_fast <- function(dat, X, out) {
  # originally from R charm package.
  beta = solve(crossprod(X), t(crossprod(dat, X)))
  return(beta[out,1])
}





###############################################################################
  # Freedman Lane Function
  #
  # Performs permutations for test for medition using the models:
  # Adjusted Model: 	Y = gamma0 + gamma1*x1 + gamma2*x2 + gamma3*m
  # Reduced Adj Model: 	Y = gamma0.r + gamma1.r*x1 + gamma2.r*x2
  # IE Model:		M = alpha0 + alpha1*x2 + alpha2*x2
  # Reduced IE Model:	M = alpha0.r + alpha2.r*x2
  #
  # Arguments:
  #   adj.design:			design matrix for adjusted model (x1, x2, m)
  #	adj.reduced.design:	reduced design matrix (x1, x2)
  #	ie.design:			design matrix for ie model (x1, x2)
  #	ie.reduced.design:	reduced design matrix (x2)
  #	adj.fitted:			fitted values, yhat, from adjusted model
  #	adj.reduced.fitted:	fitted values, yhat.r, from reduced adj model
  #	adj.resid:			residuals from adjusted model
  #	adj.reduced.resid:	residulas from reduced adj model
  #	ie.fitted:			fitted values, mhat, from IE model
  #	ie.reduced.fitted:	fitted values, mhat.r, from reduced IE model 
  #	ie.resid:			residuals from IE model
  #	ie.reduced.resid:		residuals from reduced IE model 
  #	n.greater.tx1:		number of times permuted test stat is as extreme or more 
  #					extreme than the original, tx1.orig.ols
  #	n.greater.alpha:		number of times permuted coefficient, alpha1, is as extreme or more 
  #					extreme than the original, alpha.orig
  #	n.greater.gamma:		number of times permuted coefficient, gamma3, is as extreme or more 
  #					extreme than the original, gamma.orig
  #	n.perms:			number of permutations
  #	adj.stat:			identifies the mediating variable
  # 	ie.stat:			identifies the predictor variable
  #	gamma.orig:			original coefficient from the adjusted model, gamma3
  #   alpha.orig:			original coefficient from the ie model, alpha1
  #	tx1.orig.ols:		original test stat, alpha1 * gamma3
  #
  # Returns:
  #   a list with the following:
  #	n.greater.tx1		
  #  	n.gerater.alpha
  #   n.greater.gamma
  #  	ci.tx1star.alt.ols 	95% confidence interval for alpha1*gamma3
###############################################################################

freedman_lane_sim <- function(adj.design, adj.reduced.design, ie.design, ie.reduced.design,
			 	    adj.fitted, adj.reduced.fitted, adj.resid, adj.reduced.resid,
			 	    ie.fitted, ie.reduced.fitted, ie.resid, ie.reduced.resid,
				    n.greater.tx1, n.greater.tx1.Mstar, n.greater.alpha, n.greater.gamma, 
				    n.perms, adj.stat, ie.stat, 
				    gamma.orig, alpha.orig, tx1.orig.ols) {

  # number of simulations with a stat greater than the observed.
  nc <- ncol(adj.reduced.resid)

  adj.resid  <- t(adj.resid)
  adj.fitted <- t(adj.fitted)

  adj.reduced.resid  <- t(adj.reduced.resid)
  adj.reduced.fitted <- t(adj.reduced.fitted)

  ie.resid  <- t(ie.resid)
  ie.fitted <- t(ie.fitted)

  ie.reduced.resid  <- t(ie.reduced.resid)
  ie.reduced.fitted <- t(ie.reduced.fitted)

  tx1star.alt.ols   <- NULL
  tx1star.alt.ols.i <- NULL

  tx1star.alt.ols.Mstar   <- NULL
  tx1star.alt.ols.i.Mstar <- NULL


  for(i in 1:n.perms){

	ystar.null <- NULL
	mstar.null <- NULL
	ystar.alt  <- NULL
	mstar.alt  <- NULL

	perm <- sample(1:nc)

	for(j in 1:nc){
	      ystar.null[j] <- adj.reduced.fitted[j] + adj.reduced.resid[perm[j]]
		mstar.null[j] <- ie.reduced.fitted[j] + ie.reduced.resid[perm[j]]

	      ystar.alt[j]  <- adj.fitted[j] + adj.resid[perm[j]]
		mstar.alt[j]  <- ie.fitted[j] + ie.resid[perm[j]]
		}

	mstar <- as.matrix(mstar.null)
	colnames(mstar) <- "m"
	adj.design.Mstar <- cbind(ie.design, as.matrix(mstar))

	gamma.star.null 		<- lm_fit_fast(as.matrix(ystar.null), adj.design, adj.stat)
	gamma.star.null.Mstar 	<- lm_fit_fast(as.matrix(ystar.null), adj.design.Mstar, adj.stat)
	gamma.star.alt  		<- lm_fit_fast(as.matrix(ystar.alt), adj.design, adj.stat)
	gamma.star.alt.Mstar	<- lm_fit_fast(as.matrix(ystar.alt), adj.design.Mstar, adj.stat)

	alpha.star.null <- lm_fit_fast(as.matrix(mstar.null), ie.design, ie.stat)
	alpha.star.alt  <- lm_fit_fast(as.matrix(mstar.alt), ie.design, ie.stat)

	tx1star.null.ols 		<- alpha.star.null * gamma.star.null
	tx1star.null.ols.Mstar 	<- alpha.star.null * gamma.star.null.Mstar 
	tx1star.alt.ols.i    	<- alpha.star.alt * gamma.star.alt
	tx1star.alt.ols.i.Mstar	<- alpha.star.alt * gamma.star.alt.Mstar
	
	n.greater.tx1   		<- n.greater.tx1 + (abs(tx1star.null.ols) >= abs(tx1.orig.ols))
	n.greater.tx1.Mstar	<- n.greater.tx1.Mstar + (abs(tx1star.null.ols.Mstar) >= abs(tx1.orig.ols))
	n.greater.alpha 		<- n.greater.alpha + (abs(alpha.star.null) >= abs(alpha.orig))
	n.greater.gamma 		<- n.greater.gamma + (abs(gamma.star.null) >= abs(gamma.orig))

	tx1star.alt.ols[i] <- tx1star.alt.ols.i
	tx1star.alt.ols.Mstar[i] <- tx1star.alt.ols.i.Mstar

   }

  sorted.tx1star.alt.ols <- sort(tx1star.alt.ols)
  ci.tx1star.alt.ols <- c(sorted.tx1star.alt.ols[0.025*n.perms], sorted.tx1star.alt.ols[0.975*n.perms])

  sorted.tx1star.alt.ols.M <- sort(tx1star.alt.ols.Mstar)
  ci.tx1star.alt.ols.Mstar <- c(sorted.tx1star.alt.ols.M[0.025*n.perms], sorted.tx1star.alt.ols.M[0.975*n.perms])

  return(list("n.greater.tx1"=n.greater.tx1, "n.greater.tx1.Mstar"=n.greater.tx1.Mstar, 
		"n.greater.alpha"=n.greater.alpha, "n.greater.gamma"=n.greater.gamma, 
		"ci.tx1star.alt.ols" = ci.tx1star.alt.ols, 
		"ci.tx1star.alt.ols.Mstar" = ci.tx1star.alt.ols.Mstar))

}




###############################################################################
  # coverZero: determines whether a confidence interval covers the value 0
  #
  #   ci:  a p by 2 matrix.  column 1 is the lower bound and column 2 is the
  #        upper bound.
  #
  # Return:
  #   a vector of length p of logical values, TRUE is the ci covers 0 and FALSE
  #   if the CI does not cover 0.  
###############################################################################

coverZero <- function(ci) {
  apply(ci, 1, function(x){ prod(sign(x)) < 0 })
}






#################################################################################
#################################################################################
#################################################################################
#################################################################################
tStat <- function(data, indices, tnew.ols, tx1.ols, tx2.ols, 
			tnew.lasso, tx1.lasso, tx2.lasso) {

  # Generate the needed statistics for the boot call in the mediation analysis
  # simulations
  #
  # Arguments:
  #   data:
  #   indices: 
  #
  # Returns:
  #   

  d <- data[indices, ]

  ### OLS Models ###
  # Adjusted Model
  adjOLS <- lm(y ~ x1 + x2 + m, data = d)
  gamma3.OLS <- adjOLS$coeff[4]

  # Mediation Model
  ieOLS <- lm(m ~ x1 + x2, data = d)
  alpha1.OLS <- ieOLS$coeff[2]
  alpha2.OLS <- ieOLS$coeff[3]

  origTnew.OLS <- ((abs(alpha1.OLS) + abs(alpha2.OLS)) * gamma3.OLS) - tnew.ols
  origTx1.OLS  <- (alpha1.OLS * gamma3.OLS) - tx1.ols
  origTx2.OLS  <- (alpha2.OLS * gamma3.OLS) - tx2.ols



  ### LASSO Models ###
  # Adjusted Model
  adjLASSOall <- l1ce(y ~ x1 + x2 + m, 
                      data = d,
                      bound = (1:50/50), 
                      absolute.t = FALSE, 
                      standardize = FALSE, 
				sweep.out = ~ m)
  adjTuning <- lasso2:::gcv.l1celist(adjLASSOall)
  idx <- which.min(adjTuning[, "gcv"])
  adjLASSO <- l1ce(y ~ x1 + x2 + m, 
                      data = d,
                      bound = adjTuning[idx, "rel.bound"],
                      absolute.t = FALSE, 
                      standardize = FALSE, 
				sweep.out = ~ m)
  gamma3.LASSO <- adjLASSO$coeff[4]

  # Mediation Model
  ieLASSOall <- l1ce(m ~ x1 + x2, 
                     data = d,
                     bound = (1:50/50), 
                     absolute.t = FALSE, 
                     standardize = FALSE, 
				sweep.out = ~ 1)
  ieTuning <- lasso2:::gcv.l1celist(ieLASSOall)
  idx <- which.min(ieTuning[, "gcv"])
  ieLASSO <- l1ce(m ~ x1 + x2 , 
                     data = d,
                     bound = ieTuning[idx, "rel.bound"],
                     absolute.t = FALSE, 
                     standardize = FALSE, 
				sweep.out = ~ 1)
  alpha1.LASSO <- ieLASSO$coeff[2]
  alpha2.LASSO <- ieLASSO$coeff[3]

  # Calculate Test Statistics
  origTnew.LASSO <- ((abs(alpha1.LASSO) + abs(alpha2.LASSO))*gamma3.LASSO) - tnew.lasso
  origTX1.LASSO  <- (alpha1.LASSO*gamma3.LASSO) - tx1.lasso
  origTX2.LASSO  <- (alpha2.LASSO*gamma3.LASSO) - tx2.lasso


  allT <- cbind(origTnew.OLS, origTx1.OLS, origTx2.OLS, origTnew.LASSO, origTX1.LASSO, origTX2.LASSO)
  return(allT)
}





coverT <- function(ci, testStat) {
  # deterime if a ci covers the value 0
  #
  # Args:
  #   ci:  a p by 2 matrix.  column 1 is the lower bound and column 2 is the
  #        upper bound.
  #
  # Return:
  #   a vector of length p of logical values, TRUE is the ci covers the original 
  #	test statistic and FALSE if the CI does not.  

  testStat >= ci[1] && testStat <= ci[2]

}


###################
### End of File ###
###################
