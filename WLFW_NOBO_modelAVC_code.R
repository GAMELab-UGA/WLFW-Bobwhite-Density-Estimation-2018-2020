#---------------------------------------------------------------------------------------------------------------
## SCRIPT TO FIT AC DOSER AVC MODEL TO WLFW DATA FOR NORTHERN BOBWHITES (Colinus virginianus).
## From: Nolan, V., J. Yeiser, B. Costanzo, M. Martin, J. McGuire, C. Delancey, W. Lewis, and J. Martin.
##        2024. Effects of management practices Northern Bobwhite (Colinus virginianus) density in
##        privately owned working forests across the southeastern United States. Ecological Solutions
##        and Evidence.
#--------------------------------------------------------------------------------------------------------------


require(rjags)

# Data
load("data.WLFW.NOBO.gzip")
 
# c: Formatted point count detection data. Rows correspond to sites while columns correspond to
#     repeat visits. Values represent total number of bobwhites detected per survey.
# R: Number of survey sites.
# n.count: The number of repeat point count surveys performed at each site.
# J: Number of days of recording for ARUs at each site.
# v: Formatted ARU detection data. Rows correspond to sites while columns correspond to repeat
#     ARU recordings across days. Values represent the number of bobwhite covey calls detected
#     on the recording buy the automated classifier.
# y: Formatted ARU presence/absence detection data. Rows correspond to sites while columns
#     correspond to repeat ARU recordings across days. Values represent the presence (1) or
#     absence (0) of bobwhite covey call detections on the recording by the automated classifier.
# J.r: Number of ARU recordings at each site which detected at least one bobwhite call.
# A.times: Matrix indexing the specific ARU recordings at each site that have at least one
#     bobwhite covey call detected. Rows correspond to sites and columns are a sequential list
#     of the columns in v/y which have at least one detection.
# nyrs: Number of years of data.
# nstates: Number of states in the US which data came from.
# J.A: Maximum number of ARU recordings at a single site. y-dimension of y/x.
# yr: Vector giving the year each study site was sampled (1=2018, 2=2019, 3-2020)
# state: Vector giving the state of each study sites (1=Alabama, 2=Georgia, 3=Florida, 4=North Carolina, 5=South Carolina)
# days: Matrix giving the unique day of each ARU recording. 1 corresponds to the first calendar
#       day with a recording across years, 2 corresponds to the second calendar day with a recording
#       across years, etc. 1-31 correspond roughly to October 31st to December 9th 2018,
#       32-98 correspond roughly to October 8th to December 13th 2019, and 99-157 correspond
#       roughly to October 13th - December 15th 2020.
# n.days: Total number of days with ARU recordings across years.
# acres_brush: Scaled acreage of brush management at each site.
# acres_fire: Scaled acreage of prescribed fire at each site.
# abundance_covariates: 3-D array giving scaled landscape covariates at each site. x-dimension
#       corresponds to sites, y-dimension corresponds to the percent cover of the four
#       landscape covariates (forest, herbaceous, agriculture, burned) around each point,
#       and the z-dimension corresponds to landscape scale (500, 2000, 4000, 8000, and 10000 m).
# ncovariates_abundance: The number of landscape scales assessed. Corresponds to length of
#       of z-dimension in abundance_covariates.
# ydb: Formatted distance-detection point count data. The x-dimension corresponds to sites,
#       y-dimension corresponds to distance bands of detections, and z-dimension corresponds
#       to repeat visits. Values are the number of bobwhite detections on point counts in
#       each distance bands.
# db: Boundaries of distance bands used in the point count distance sampling model.
# pix: Relative area within each distance band used in the point count distance sampling model.
# nB: Number of distance bands used in the point count distance sampling model.
# mean_omega: Mean number of false positives per ARU recording. Assessed via manual validation 
#       of a subset of ARU detections.
# point.areaARU: circular area (square meters) surveyed by ARUs.
# point.areaPC: circular area (square meters) surveyed by point counts.

# Initial Values
set.seed(12)
N.init <- rep(1, data.WLFW.NOBO$R)
c.max <- apply(data.WLFW.NOBO$c, 1, max, na.rm = TRUE)
N.init <- ifelse(c.max > 1, c.max + 1, 1)
inits <- function() {
  list(
    N_ARU = N.init,
    N_PC = N.init,
    mean_alpha0 = runif(1,-1,1),
    sd_alpha0 = runif(1,0.001,10),
    mean_beta0 = runif(1,0,10),
    sd_beta0 = runif(1,0.001,10),
    beta0 = array(runif(1,0,10), dim=c(data.WLFW.NOBO$nyrs,data.WLFW.NOBO$nstates)),
    sd_omega = runif(1,0.001,10),
    tau.day = runif(1, 0.1, 1),
    a.phi = runif(1, 0.5, 1),
    mu.p.aru = runif(1, 0, 1),
    p.aru.1 = runif(1,0,1)
  )
}

# Parameters monitored
params <- c('abundance_scale','mean_alpha0','prec_alpha0','alpha.0','alpha1','alpha2','alpha6','alpha6q','alpha7','alpha7q','alpha8','alpha8q','alpha9','alpha9q',
                'p.aru.0','p.aru.1','mu.p.aru','mean_omega','tau.day','a.phi','mean_beta0','prec_beta0','beta0','pd',"gamma.1","phi","delta","omega")

# MCMC settings
n.iter <- 120000
n.thin <- 25
n.burn <- 10000
n.chain <- 3
n.adapt <- 10000


sink("WLFW_NOBO_integratedmodel.jags")
  cat("
      model{
      
      #-------------------------------------------------------------------------
      ## Priors ##
      #-------------------------------------------------------------------------


      ###### ABUNDANCE ######

      ## Global Abundance Intercept ##
     	mean_alpha0 ~ dnorm(0, 0.1)    	        # Global mean abundance
     	sd_alpha0 ~ dunif(1, 10)                # Standard deviation
     	prec_alpha0 <- 1/(sd_alpha0^2)  	# Precision
     	for(t in 1:nyrs){  
       	for(s in 1:nstates){
          alpha.0[t,s] ~ dnorm(mean_alpha0, prec_alpha0)
       	}#s
      }#t

      ## Covariates on abundance - model includes both linear and quadratic effects for the landscape-scale covariates (3-6) ##

		  #1) effect of management - acres brush
			alpha1 ~ dnorm(0, 0.001)
		  #2) effect of management - acres fire
			alpha2 ~ dnorm(0, 0.001)
		  #6) habitat 1 - forest 
			alpha6 ~ dnorm(0, 0.001)
			alpha6q ~ dnorm(0, 0.001)
		  #7) habitat 2 - herb 
			alpha7 ~ dnorm(0, 0.001)
			alpha7q ~ dnorm(0, 0.001)
		  #8) habitat 3 - agri 
    	alpha8 ~ dnorm(0, 0.001)
			alpha8q ~ dnorm(0, 0.001)
		  #9) Fire - fire frequency
    	alpha9 ~ dnorm(0, 0.001)
			alpha9q ~ dnorm(0, 0.001)

      ## 15 scales of landscape-scale covariates ranging from 500m - 10km ##

 	 	  for (c in 1:(ncovariates_abundance-2)){
			  # 6 scales for covariates
			  abundance_scale[c] ~ dcat(c(1/6, 1/6, 1/6, 1/6, 1/6,1/6))
		  }#c   


      ###### DETECTION - POINT COUNTS ######

      ## PC Detection - Distance Sampling ##
		  mean_beta0 ~ dunif(0, 20)	                # Global mean detection from point counts
		  sd_beta0 ~ dunif(1, 10)		        	# Standard deviation
		  prec_beta0 <- 1/(sd_beta0^2)		        # Precision
	   	for(t in 1:nyrs){
    	  for(s in 1:nstates){
    		  beta0[t,s] ~ dnorm(mean_beta0,prec_beta0)
    	  }#s
    	}#t


      ###### DETECTION - ARUS ######

      ## Mean number of false positive acoustic detections ##
 	 	  sd_omega ~ dunif(0, 10)
      prec_omega <- 1/(sd_omega^2)
	   	for(t in 1:nyrs){
    	  for(s in 1:nstates){ 
  	 	    omega[t,s] ~ dnorm(mean_omega,prec_omega)
    	  }#s
  	 	 }#t   
  
  	 	mu.p.aru ~ dunif(0,1)    		 # mu.alpha for p.aru.0 prior
 	  	p.aru.0 <- logit(mu.p.aru)     		 # prob (logit scale) of detecting at least one vocalization at an unoccupied site
  	  p.aru.1 ~ dunif(0, 1000)      		 # Additional prob (logit scale) of detecting at least one vocalization for every 1 increase in covey abundance. Constrained to be positive
		  tau.day ~ dunif(0.01, 1)    	         # Precision for random day effect on true vocalization/detection rate. 
  	  a.phi ~ dunif(0,100)          	         # Overdispersion parameter for zero-truncated negative binomial. 

      ## Loop over each recording day - gamma.1 = random day effect on true vocalization detection rate
    	for (i in 1:n.days) {                
       	gamma.1[i] ~ dnorm(0, tau.day)
      }#i

      ## Loop over number of total sites (R) and max number of repeat visits at each acoustic data site (J.A) 
    	for (i in 1:R) {
    	  for (j in 1:J.A) {
		      phi[i, j] ~ dnorm(a.phi, 100)
		    }#j
       }#i



      #-------------------------------------------------------------------------
      ## Liklihood and process model ##
      #-------------------------------------------------------------------------
  
      for (i in 1:R) {

        N_PC[i] ~ dpois(lambda_PC[i])   # Number of birds from PC additional area
        log(lambda_PC[i]) <- alpha.0[yr[i],state[i]] + alpha1*acres_brush[i] + alpha2*acres_fire[i] + 	 
 				                     	alpha6 * abundance_covariates[i,1,abundance_scale[1]] +
 				 		                  alpha6q * pow(abundance_covariates[i,1,abundance_scale[1]],2) +
 				 		                  alpha7 * abundance_covariates[i,2,abundance_scale[2]] +
 				 		                  alpha7q * pow(abundance_covariates[i,2,abundance_scale[2]],2) +
 						                  alpha8 * abundance_covariates[i,3,abundance_scale[3]] +
 				 		                  alpha8q * pow(abundance_covariates[i,3,abundance_scale[3]],2) +
 				 		                  alpha9 * abundance_covariates[i,4,abundance_scale[4]] +
 				 		                  alpha9q * pow(abundance_covariates[i,4,abundance_scale[4]],2) +
 				 	 	                  log(point.areaPC/10000)

        N_ARU[i] ~ dpois(lambda_ARU[i])   # Number of birds
        log(lambda_ARU[i]) <- alpha.0[yr[i],state[i]] + alpha1*acres_brush[i] + alpha2*acres_fire[i] +
 	  					                 alpha6 * abundance_covariates[i,1,abundance_scale[1]] +
 				 		                   alpha6q * pow(abundance_covariates[i,1,abundance_scale[1]],2) +
 				 		                   alpha7 * abundance_covariates[i,2,abundance_scale[2]] +
 				 		                   alpha7q * pow(abundance_covariates[i,2,abundance_scale[2]],2) +
 						                   alpha8 * abundance_covariates[i,3,abundance_scale[3]] +
 				 		                   alpha8q * pow(abundance_covariates[i,3,abundance_scale[3]],2) +
 				 		                   alpha9 * abundance_covariates[i,4,abundance_scale[4]] +
 				 		                   alpha9q * pow(abundance_covariates[i,4,abundance_scale[4]],2) +
 						                   log(point.areaARU/10000)

						
        ## Prob of detecting at least one vocalization in an acoustic recording 
        logit(p.a[i]) <-  (p.aru.0 + p.aru.1 * N_ARU[i])



        ##### Acoustic Data #####

        ## Loop over J = Number of repeat visits for acoustic data at each site
    	  for (j in 1:J[i]) {
      	  log(delta[i, j]) <- gamma.1[days[i, j]]
      		y[i, j] ~ dbin(p.a[i], 1)
        } #j


        ## Loop over max number of visits per site for v (acoustic vocalization data from clustering algorithm)
        ## A.times = indexing variable used to determine specific indexes with v > 0. 
   	    for (j in 1:J.r[i]) {
     		  v[i, A.times[i, j]] ~ dpois((delta[i, A.times[i, j]] * N_ARU[i] + omega[yr[i],state[i]]) * phi[i, A.times[i, j]] * y[i, A.times[i, j]]) 					T(1, )
	      } #j



        ##### Count Data #####

        ## Loop over number of repeat visits for count data for each site
   	    for (j in 1:n.count[i]) {
   		    pMarg[i,j] <- pd[yr[i],state[i]]
    		  ydb[i,1:nB,j] ~ dmulti(pd_bins_adj[1:nB,yr[i],state[i]], c[i,j])
    		  c[i,j] ~ dbin(pMarg[i,j], N_PC[i])
   	    }#j 
   	  }#i

   	  for(t in 1:nyrs){
   	    for(s in 1:nstates){
    	    log(sigma[t,s]) <- beta0[t,s]
       	  for(b in 1:nB){
            pd_bins[b,t,s] <- (sigma[t,s]^2*(1-exp(-db[b+1]^2/(2*sigma[t,s]^2)))-sigma[t,s]^2*(1-exp(-db[b]^2/
  	   					                (2*sigma[t,s]^2))))*2*3.1416/(point.areaPC*pix[b])
           	pd_bins_adj[b,t,s] <- pd_bins[b,t,s]*pix[b]+0.0000001
           }#b  
   	   	  pd[t,s] <- sum(pd_bins_adj[1:nB,t,s])
  	  	}#s
      }#t

  }
  ", fill=TRUE)
sink()


## Takes a few minutes to complete
out.model.A.V.C <- jags.model(data=data.WLFW.NOBO, inits = inits, file='WLFW_NOBO_integratedmodel.jags',n.chain = n.chain, n.adapt = n.adapt)
cs1 <- coda.samples(out.model.A.V.C, params, thin=n.thin, n.iter=n.iter, n.burn=n.burn)
summary(cs1)
gelman.diag(cs1, multivariate = FALSE) # All should be below 1.1