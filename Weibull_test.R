
#' Weibull tests
#'
#' @param data 1st column contains timepoints, 2nd column contains status
#' @param censor Largest time or maximum of observation period
#'
#' @return A table summarising the results of the corresponding Weibull tests
#' @export
#'
#' @examples
weibull_test = function(data, censor = max(data[,1], na.rm = T)){
  if(missing(data)){stop("Argument data has to be provided.")}
  if(!((is.data.frame(data))||(is.matrix(data)))){stop("Argument data has to be a matrix or data frame with the time-to-event in the first column and the event in the second column.")}
  if(ncol(data)<2) stop("Argument data has to be a matrix or data frame with the time-to-event in the first column and the event in the second column.")
  if(ncol(data)>2){ warning('Data matrix contains more than 2 columns. Are the columns in correct order?') }

  # censor data for censored weibull test
  data.c = data
  data.c[data$time > censor/2,1] = ceiling(censor/2)
  data.c[data$time > censor/2,2] = 0

  ## model estimation
  res.pgW = try(stats::nlm(logl,rep(0,3),data,hessian=T)) # power generalised Weibull
  res.W = summary(survival::survreg(survival::Surv(time = data$time, event = data$status)~1, dist = "weibull")) # Weibull
  res.c.W = summary(survival::survreg(survival::Surv(time = data.c$time, event = data.c$status)~1, dist = "weibull")) # Weibull on censored data

  ## tests
  # pgW test
  rej.pgW = test.pgW(res.pgW)
  rej.W = test.W(res.W)   # Weibull-test
  rej.c.W = test.W(res.c.W)   # Weibull-test on the cencored data
  rej.dW = rej.W + rej.c.W >= 1   # double Weibull
  rej.dWpW = rej.W + rej.c.W + rej.pgW >= 1  # dWSP-pWSP

  ## formatting
  out.W = naming(rej.W)
  out.c.W = naming(rej.c.W)
  out.pgW = naming(rej.pgW)
  out.dW = naming(rej.dW)
  out.dWpW = naming(rej.dWpW)

  # Significance levels
  levels = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10)

  out.table = as.data.frame(rbind(levels, out.W, out.c.W, out.pgW, out.dW, out.dWpW))
  rownames(out.table) = c("Significance level", "WSP", "cWSP", "pWSP", "dWSP", "dWSP-pWSP")
  colnames(out.table) = NULL

  cat("Weibull Shape Parameter tests for signal detection:\n")
  print(out.table)

  if(sum(is.na(rej.pgW))>1) cat("\nTest failed for the power generalized Weibull parameters due to an error in the numerical optimization.")

  invisible(out.table)
}





# Helper functions --------------------------------------------------------


################################
#### log-Likelihood for pgW ####
################################

# Input
#  ## par = Parameters for the pgW
#         = log(theta) (scale parameter),
#           log(nu) (1st shape parameter),
#           log(gamma) (2nd shape parameter)
#     -> theta, nu, gamma >0, therefore transformation via exp()
#  ## dat = 1st column contains timepoints, 2nd column contains status

# further explanation
#  ## x = time and status of one individual

# output
# ## value of the minus-loglikelihood value

logl = function(par, dat){
  theta = exp(par[1])
  nu = exp(par[2])
  gamma = exp(par[3])
  dens.survi = function(x){
    ((nu/gamma)*(x[1]^(nu-1)/theta^nu)*(1+(x[1]/theta)^nu)^(1/gamma -1)*exp(1 - (1 + (x[1]/theta)^nu)^(1/gamma)))^x[2]*(exp(1 - (1 + (x[1]/theta)^nu)^(1/gamma)))^(1-x[2])
  }
  -sum(log(apply(dat,1,dens.survi)))
}

###############################################
#### pgW test function for multiple alphas ####
###############################################

# function performs a test on both shape parameters of the power generalized Weibull distribution
# Each test contains a test for nu and a test for gamma
# and is performed 10 times for 10 different significance levels, alphas = 0.01,...,0.1

# Input: eo = estimation output

# output: vector of length 10 containing the test results for each alpha
# (alpha increases from left to right)


test.pgW = function(eo){
  if(is.vector(eo) == F){
    rej.H0 = rep(NA, 10)
  } else if(is.numeric(eo$hessian) == F){
    rej.H0 = rep(NA, 10)
  } else{
    varmatrix = try(solve(eo$hessian))
    # print(varmatrix)
    if(is.matrix(varmatrix) == F){
      rej.H0 = rep(NA, 10)
    }
    else if(sum(is.na(varmatrix)) > 0){ #added
      rej.H0 = rep(NA, 10)
    }
    else if(sum(is.nan(varmatrix)) > 0){ #added
      rej.H0 = rep(NA, 10)
    }
    else if(varmatrix[2,2] < 0 | varmatrix[3,3] < 0){
      rej.H0 = rep(NA, 10)
    }
    else{
      test.pgW.inside = function(alpha){
        CI.nu = eo$estimate[1] + c(-1,1)*qnorm(1-alpha/2)*sqrt(varmatrix[2,2])
        rej.nu = as.numeric(0 <= CI.nu[1] | CI.nu[2] <= 0) #here: 1 if rejected, 0 if not rejected
        CI.gamma = eo$estimate[2] + c(-1,1)*qnorm(1-alpha/2)*sqrt(varmatrix[3,3])
        rej.gamma = as.numeric(0 <= CI.gamma[1] | CI.gamma[2] <= 0)
        return(c(rej.nu, rej.gamma))
      }
      alphas = (1:10)/100
      rej = sapply(alphas, test.pgW.inside)
      rej.nu = rej[1,]
      rej.gamma = rej[2,]
      rej.H0 = rej.nu*rej.gamma #here: 1 if H0 rejected, 0 if H0 not rejected
    }
    return(pgWtest = rej.H0)
  }
}



#################################################
#### Weibull test for multiple alphas ####
#################################################

# function performs a test on the shape parameter of the weibull function
# H0: shape parameter does not differ from 1

# input: summary of a survreg with Weibull distribution, no predictors

# output: vector of length 10 containing the test result for different significance
# levels reaching from 0.01 to 0.1 (equidistant steps
# if H0 is rejected, then we get a 1, else a 0

test.W = function(est){
  test.W.inside = function(alpha){
    rej.W = as.numeric(try(isTRUE(est[[9]][2,4] < alpha)))
    return(rej.W)
  }
  alphas = 1:10/100
  rej.H0 = sapply(alphas, test.W.inside)
  return(Wtest = rej.H0)
}

##################################
#### data generating function ####
##################################

# input parameters (vector with 5 entries):
#   ##sample size (n)
#   ##backround rate (br)
#   ##ADR rate as proportion of the number of backround events (adr)
#     -> mean of the ADR occuring timepoint is generated randomly
#        via a uniform distribution
#   ##relative standard deviation of the adr generating process (rel.sd)
#     -> sd for the generating process = rel.sd*censor
#   ##end of the observation period (censor)
#     (start assumed to be at 0)

#output is a list containing:
#   ##data = simulated dataframe containing status and time

datagen <- function(genpar){
  n = genpar[1]
  br = genpar[2]
  adr = genpar[3]
  rel.sd = genpar[4]
  censor = genpar[5]
  # Number of br & adr cases in the study
  n.br = rbinom(1,n, prob = br)
  n.adr = rbinom(1,n, prob = br*adr)
  # event time for backround event candidates
  t.br = runif(n.br, min = 0, max = censor)
  #t.br=rexp(n.br,ln(br/censor))
  # mean timepoint for adr events
  m.adr = runif(1,0,censor)
  sd.adr = rel.sd*censor
  # event time for adr event candidates
  t.adr = rnorm(n.adr, mean = m.adr, sd = sd.adr)
  #print(c(m.adr,sd.adr))
  # adjustment for negativ timepoints
  if(sum(t.adr<=0)>0){
    t.adr[t.adr<=0] = 1
  }
  # number of events in the observation period
  n.events = length(c(t.br, t.adr))
  # status vector for whole sample
  status = c(rep(1,n.events), rep(0, n-n.events))
  # time vector (continuous values) for whole sample
  time = c(t.br, t.adr, rep(censor,n-n.events))
  # time vector (daily scale) for whole sample
  time = ceiling(time)
  dat = data.frame(time, status)
  # print(c(length(t.br), length(t.adr),table(status)))
  return(dat)
}


naming = function(result){
  out = rep(NA, length(result))
  out[which(result == 1)] = "signal"
  out[which(result == 0)] = "no signal"
  return(out)
}



# Testing -----------------------------------------------------------------

n = c(100,200,300,400)
br = c(0.01, 0.05, 0.1)
adr = c(1, 0.5)
rel.sd = c( 0.05)
censor = 365
pars = as.matrix(expand.grid(n=n,br=br,adr=adr,rel.sd=rel.sd,censor=censor))

data = datagen(pars[12,])
test_result = weibull_test(data)

