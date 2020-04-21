#' Generating function for Weibull Data.
#'
#' \code{genweibulldata()} function used mainly internally by
#' \code{weibulltrialsimulator()} and \code{weibulltrialsimulatornohist()} functions
#'  to generate data for a two-arm clinical trial, experimental and control groups.
#'  Can be used to generate random trial data.
#'
#' @param sample_size  Number of subjects per arm.
#' @param scale1       Scale parameter used in call to \code{rweibull()}.
#'   Used only in control arm.
#' @param hazard_ratio Desired Hazard Ratio between experimental and control groups.
#' @param common_shape Shape parameter used in call to \code{rweibull()}.
#'   Used in both arms.
#' @param censor_value Value at which time-to-event data are right censored.
#'
#' @return \code{genweibulldata()} returns a data frame with columns: 'id', 'treatment',
#'   'event_time', and 'status'.
#'
#' @examples
#' SampleHistData <- genweibulldata(sample_size=60, scale1=2.82487,
#'                                  hazard_ratio=0.6, common_shape=3,
#'                                  censor_value=3)
#' SampleHistData
#' @export
genweibulldata <- function(sample_size, scale1, hazard_ratio, common_shape, censor_value) {

  # --------------------------------------------------------------- #
	# The function genweibulldata simulates a balanced clinical trial
	# with 'sample_size' subjects per arm using a weibull distribution.
	# 'scale1' is the weibull scale parameter, and 'common_shape' is the
	# weibull common shape parameter for both arms, where scale and
	# shape parameters are defined according to the rweibull() function
  # in R.  'censor_value' is the value when right censoring occurs. As
  # of 9/5/2016, genweibulldata only generates data with right
	# censoring.  Random right censoring is not incorporated.
	# 'hazard_ratio is the ratio of group hazards (experimental group
	# over control group).
	#
	# In the code below time1, scale1, test1, status1, etc. are data
	# for the control goup.
	# In the code below time2, scale2, test2, status2, etc. are data
	# for the experimental group.
  # --------------------------------------------------------------- #

	#Define experimental group scale parameter given control group
	# scale parameter, commom shape parameter in both groups, and
	# the user specified hazard ratio.
    scale2 <- exp(log(scale1) - (1/common_shape) * log(hazard_ratio))

    # Create event times for both groups
    time1 <- stats::rweibull(sample_size, shape = common_shape, scale = scale1)
    time2 <- stats::rweibull(sample_size, shape = common_shape, scale = scale2)

    # Create variables needed for simulation, if censor_value is specified.
    if (!is.null(censor_value) == TRUE) {
        test1 <- (time1 > censor_value) #Identify which times need right censoring.
        test2 <- (time2 > censor_value)
        status1 <- rep(1, sample_size) #Initialize the Status variable.
        status2 <- rep(1, sample_size)
		#For all observations that need to be right censored, set the time
		# value to the right censor value and set the status value to 0
		# (indicating right censoring).  Status=1 implies observed event.
		# Status=0 implies right censored event.
        time1[test1] <- censor_value
        time2[test2] <- censor_value
        status1[test1] <- 0
        status2[test2] <- 0
    }
	#Create status variable if censor_value is not specified (in such a case
	# status = 1 for all observeations.)
    if (is.null(censor_value) == TRUE) {
        status1 <- rep(1, sample_size)
        status2 <- rep(1, sample_size)
    }

	#Take all data created above and put into a data frame that contains
	# the required variables.
    subjid <- seq(from = 1, to = 2 * sample_size)
    trt <- c(rep(0, sample_size), rep(1, sample_size))
    time <- c(time1, time2)
    status <- c(status1, status2)

    gendata <- data.frame(subjid, trt, time, status)
    colnames(gendata) <- c("id", "treatment", "event_time", "status")

    return(gendata)
}

#' Log-likelihood function for two-arm trial with historical data using Weibull distribution.
#'
#' \code{weibullloglike()} function only used internally by
#' \code{weibulltrialsimulator()} function to estimate Weibull model parameters
#' when clinical trial involves experimental and control groups as well as historical
#' control data.
#'
#' The Weibull log-likelihood is calculated by using \code{dweibull()}.  The following
#' derivation was used to reparametrize the model parameters so that \code{dweibull()}
#' could be used in the log-likelihood function and still estimate the hazard ratio as a
#' parameter.
#'
#' If we let \eqn{scale = exp(-beta0/shape - beta1*treatment/shape)}, then that for the control
#' group: \eqn{scale0 = exp(-beta0/shape)}.  Similarly, for the experimental group
#' \eqn{scale1 = exp(-beta0/shape - beta1/shape)}.
#'
#' Now for the Weibull distribution \deqn{hazard(t) = (shape)*(scale^(-shape))*(t^(shape-1))}
#' From this equation, we can derive the Weibull hazard ratio, HR, (experimental over control) as a
#' function of \eqn{scale0}, \eqn{scale1}, and the shape parameter \eqn{HR = (scale0/scale1)^shape}.
#'
#' Substituting for \eqn{scale0} and eqn{scale1} using functions of \eqn{beta0}, \eqn{beta1}, and
#' \eqn{shape}, we can express the hazard ratio in terms of \eqn{beta0} and \eqn{beta1}. When
#' we do this we have \deqn{HR = (exp(-beta0/shape)/exp(-beta0/shape - beta1/shape))^shape}
#' After reducing terms we have a simply formula for the hazard ratio,
#' \eqn{HR = (exp(beta1/shape))^shape}, which reduces further to \eqn{exp(beta1)}.
#'
#' Therefore, \eqn{beta1} is the log hazard ratio of experimental over control.  Similarly
#' \eqn{log(scale0) = -beta0/shape} so the nuisance parameter \eqn{beta0} is equal to
#' \eqn{(-shape)log(scale0)}.
#'
#' \code{weibullloglike()} should not be called directly by user.
#'
#' @param params  Three element vector of Weibull parameters.  The third element is
#'   the shape parameter used in \code{dweibull()}.  The first and second elements
#'   are the intercept (beta0), and treatment effect (beta1), parameters as defined in
#'   details section.  The beta1 parameter is the log hazard ratio.
#' @param randdata  Dataset of randomly generated trial data.  Randomized trial datasets
#'   must have 4 columns: id, treatment, event_time, and status.  The value of treatment
#'   must be 0 (control) or 1 (experimental).  The values of event_time must be positive.
#'   The values of status must be 0 (right censored event) or 1 (observed event).
#' @param histdata Dataset of historical control data.  Historical datasets must have 4 columns:
#'   id, treatment, event_time, and status.  The value of treatment should be 0.  The
#'   values of event_time must be positive.  The values of status must be 0 (right
#'   censored event) or 1 (observed event).
#' @param a0 Power prior parameter where 0 implies historical control data is ignored and 1 implies
#'   all information in historical control data is used.  A value between 0 and 1 partially includes
#'   the historical control data.
#'
#' @return \code{weibullloglike()} returns a value of the loglikelihood function
#'   given a set of Weibull parameters, randomly generated trial data, and observed
#'   historical control data.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
weibullloglike <- function(params, randdata, histdata, a0) {

  # --------------------------------------------------------------- #
	#  This function calculates the Weibull log-likelihood given
	#  a vector of parameter values, a dataset of randomized trial
	#  data (two arms, no covariates beyond treatment), and a dataset
	#  of historical control data.
	#  The weibull shape parameter is common in both randomized groups
	#  and the historical control group.  The scale parameter is
	#  assumed to be the same in both control groups.  The scale
  #  parameter for the randomized experimental group is a function
  #  of the control scale and the treatment effect.  The parameters
  #  are beta0, beta1, and v.  v is the common shape parameter.
  #  beta0 and beta1 are regression parameters that are linked to
	#  the Weibull scale parameter via the exp() function.
	#  beta1 is the log hazard ratio (experimental group over control
	#  group), while beta0 = -v*log(scale0).
  #  Note, scale0 is the scale parameter for controls.
  # --------------------------------------------------------------- #

    # Get params
    beta0 <- params[1]
    beta1 <- params[2]
    v <- params[3]

	# Calculate the scale parameter vector for all randomized observations.
    b_i <- exp((-1 * beta0/v) + (-1 * beta1/v) * randdata$treatment)
	# Calculate the log-likelihood values for all randomized observations.
    ll_r <- randdata$status * stats::dweibull(randdata$event_time, shape = v, scale = b_i, log = TRUE) + (1 - randdata$status) *
        stats::pweibull(randdata$event_time, shape = v, scale = b_i, log.p = TRUE, lower.tail = FALSE)

	# Calculate the scale parameter vector for all historical controls.  All values in this vector are the same.
	# Note that bh_i is the same as the randomized control scale parameter for randomized observations.
    bh_i <- exp(-1 * beta0/v)
	# Calculate the log-likelihood values for all historical control observations.
    ll_h <- histdata$status * stats::dweibull(histdata$event_time, shape = v, scale = bh_i, log = TRUE) + (1 - histdata$status) *
        stats::pweibull(histdata$event_time, shape = v, scale = bh_i, log.p = TRUE, lower.tail = FALSE)

	# Calculate the overall log-likelihood by adding the randomized log-likelihood to the historical control
	# log-likelihood by a0, where a0 is the power prior parameter.  This a0 value is defined by the
	# user and not estimated via object function optimization.
    ll <- sum(ll_r) + a0 * sum(ll_h)

	# Return the sum of all individual elements to the negative log-likelihood
    return(-ll)
}

#' Log-likelihood function for two-arm trial with no historical data using Weibull distribution.
#'
#' \code{weibullloglikenohist()} function only used internally by
#' \code{weibulltrialsimulatornohist()} function to estimate Weibull model parameters
#' when clinical trial involves experimental and control groups but no historical control
#' data.
#'
#' The Weibull log-likelihood is calculated by using \code{dweibull()}.  The following
#' derivation was used to reparametrize the model parameters so that \code{dweibull()}
#' could be used in the log-likelihood function and still estimate the hazard ratio as a
#' parameter.
#'
#' If we let \eqn{scale = exp(-beta0/shape - beta1*treatment/shape)}, then that for the control
#' group: \eqn{scale0 = exp(-beta0/shape)}.  Similarly, for the experimental group
#' \eqn{scale1 = exp(-beta0/shape - beta1/shape)}.
#'
#' Now for the Weibull distribution \deqn{hazard(t) = (shape)*(scale^(-shape))*(t^(shape-1))}
#' From this equation, we can derive the Weibull hazard ratio, HR, (experimental over control) as a
#' function of \eqn{scale0}, \eqn{scale1}, and the shape parameter \eqn{HR = (scale0/scale1)^shape}.
#'
#' Substituting for \eqn{scale0} and eqn{scale1} using functions of \eqn{beta0}, \eqn{beta1}, and
#' \eqn{shape}, we can express the hazard ratio in terms of \eqn{beta0} and \eqn{beta1}. When
#' we do this we have \deqn{HR = (exp(-beta0/shape)/exp(-beta0/shape - beta1/shape))^shape}
#' After reducing terms we have a simply formula for the hazard ratio,
#' \eqn{HR = (exp(beta1/shape))^shape}, which reduces further to \eqn{exp(beta1)}.
#'
#' Therefore, \eqn{beta1} is the log hazard ratio of experimental over control.  Similarly
#' \eqn{log(scale0) = -beta0/shape} so the nuisance parameter \eqn{beta0} is equal to
#' \eqn{(-shape)log(scale0)}.
#'
#' \code{weibullloglike()} should not be called directly by user.
#'
#' @param params  Three element vector of Weibull parameters.  The third element is
#'   the shape parameter used in \code{dweibull()}.  The first and second elements
#'   are the intercept (beta0), and treatment effect (beta1), parameters as defined in
#'   details section.  The beta1 parameter is the log hazard ratio.
#' @param randdata  Dataset of randomly generated trial data. Randomized trial datasets
#'   must have 4 columns: id, treatment, event_time, and status.  The value of treatment
#'   must be 0 (control) or 1 (experimental).  The values of event_time must be positive.
#'   The values of status must be 0 (right censored event) or 1 (observed event).
#'
#' @return \code{weibullloglikenohist()} returns a value of the loglikelihood function
#'   given a set of Weibull parameters and randomly generated trial data.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
weibullloglikenohist <- function(params, randdata) {

  # --------------------------------------------------------------- #
  #  This function calculates the Weibull log-likelihood given
  #  a vector of parameter values, and a dataset of randomized trial
  #  data (two arms, no covariates beyond treatment).  Historical
  #  data is not utilized in this log-likelihood function.
  #  The weibull shape parameter is common in both randomized groups.
  #  The scale parameter for the randomized experimental group is a
  #  function of the control scale and the treatment effect.
  #  The parameters are beta0, beta1, and v.  v is the common shape
  #  parameter.  beta0 and beta1 are regression parameters that are
  #  linked to the Weibull scale parameter via the exp() function.
  #  beta1 is the log hazard ratio (experimental group over control
  #  group), while beta0 = -v*log(scale0).
  #  Note, scale0 is the scale parameter for controls.
  # --------------------------------------------------------------- #

    # Get params
    beta0 <- params[1]
    beta1 <- params[2]
    v <- params[3]

    # Calculate the scale parameter vector for all randomized observations.
    b_i <- exp((-1 * beta0/v) + (-1 * beta1/v) * randdata$treatment)
    # Calculate the log-likelihood values for all randomized observations.
    ll_r <- randdata$status * stats::dweibull(randdata$event_time, shape = v, scale = b_i, log = TRUE) + (1 - randdata$status) *
        stats::pweibull(randdata$event_time, shape = v, scale = b_i, log.p = TRUE, lower.tail = FALSE)

	# Return the sum of all individual elements to the negative log-likelihood
    return(sum(-ll_r))
}


#' Simulate a single randomized trial using a Weibull outcome and information from
#' historical controls.
#'
#' \code{weibulltrialsimulator()} function only used internally by
#' \code{weibull_sim()} function to run a single trial simulation involving historical
#' control data and a Weibull outcome.
#'
#' The simulation of a trial with a Weibull outcome involving historical control data returns
#' an estimate of the hazard ratio as well as an estimate of the log hazard ratio variance.
#' Finally the simulation returns an indication of whether or not the simulated trial led to
#' a rejection of the null hypothesis (1) or not (0).
#'
#' \code{weibulltrialsimulator()} should not be called directly by user.
#'
#' @param sample_size_val Number of subjects per arm.
#' @param histdata Dataset of historical data.  Historical datasets must have 4 columns:
#'   id, treatment, event_time, and status.  The value of treatment should be 0.  The
#'   values of event_time must be positive.  The values of status must be 0 (right
#'   censored event) or 1 (observed event).
#' @param scale1_val Randomized control arm scale parameter used in call to \code{rweibull()}.
#' @param hazard_ratio_val Desired hazard ratio between randomized experimental and control arms.
#' @param common_shape_val Randomized shape parameter used in call to \code{rweibull()}.
#'   Used in both randomized arms.
#' @param censor_value Value at which time-to-event data are right censored.
#' @param a0_val A power prior parameter ranging from 0 to 1, where 0
#'   implies no information from historical data should be used, 1 implies all of
#'   the information from historical data should be used.  A value between 0 and 1
#'   implies that a proportion of the information from historical data will be used.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#'
#' @return \code{weibulltrialsimulator()} returns a vector of simulation results. The
#'   first element is an estimated hazard ratio, the second element is the estimated
#'   variance of the log hazard ratio, and the third element is a 0/1 variable indicator
#'   whether or not the trial rejected the null hypothesis (1) or failed to reject
#'   the null hypothesis (0).
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
weibulltrialsimulator <- function(sample_size_val, histdata, scale1_val, hazard_ratio_val, common_shape_val, censor_value,
    a0_val, alpha) {

    # --------------------------------------------------------------- #
    #  This function simulates a two-arm Bayesian trial where
    #  historical data is utilized in the parameter estimation.
    # --------------------------------------------------------------- #

    # First, Generate weibull trial data given the user defined trial characteristics.
    sampleranddata <- genweibulldata(sample_size = sample_size_val, scale1 = scale1_val, hazard_ratio = hazard_ratio_val,
        common_shape = common_shape_val, censor_value = censor_value)
    # Make sure the trial data has at least one not right censored observation.
    if (sum(sampleranddata$event_time == censor_value) == dim(sampleranddata)[1]) {
        stop("Simulated trial data must have at least one observation that is not right censored.")
    }
    # Generate initial values for your call to optim()
    initializemodel <- survival::survreg(survival::Surv(event_time, status) ~ treatment, dist = "weibull", data = sampleranddata)

    initialbeta0 <- -1 * initializemodel$coefficients[1]/initializemodel$scale
    initialbeta1 <- -1 * initializemodel$coefficients[2]/initializemodel$scale
    initialv <- 1/initializemodel$scale

    # Generate the Bayesian CLT based parameter estimates needed for inference on hazard ratio.
    fitmod <- stats::optim(c(initialbeta0, initialbeta1, initialv), weibullloglike, randdata = sampleranddata, histdata = histdata,
        a0 = a0_val, method = "Nelder-Mead", hessian = TRUE)

    #Extract model parameters and statistics
    modparm <- fitmod$par
    covarmat <- solve(fitmod$hessian)

    weibullloghazard_ratio <- modparm[2]

    weibullhazard_ratio <- exp(weibullloghazard_ratio)
    lower_weibullhazard_ratio <- exp(weibullloghazard_ratio - stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2]))
    upper_weibullhazard_ratio <- exp(weibullloghazard_ratio + stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2]))

    #Make a decision about the simulated trial, reject or fail to reject null hypothesis.
    reject <- ifelse(((lower_weibullhazard_ratio > 1) | (upper_weibullhazard_ratio < 1)), 1, 0)
    output <- c(weibullhazard_ratio, covarmat[2, 2], reject)

    #Return the hazard ratio, the estimated variance of the log hazard ratio, and the trial decision.
    names(output) <- c("hazard_ratio", "loghazard_ratio_var", "reject")
    return(output)

}


#' Simulate a single randomized trial using a Weibull outcome but not including any information from
#' historical controls.
#'
#' \code{weibulltrialsimulatornohist()} function only used internally by
#' \code{simple_weibull_sim()} function to run a single trial simulation involving
#' a Weibull outcome but no historical control data.
#'
#' The simulation of a trial with a Weibull outcome involving no historical control data returns
#' an estimate of the hazard ratio as well as an estimate of the log hazard ratio variance.
#' Finally the simulation returns an indication of whether or not the simulated trial led to
#' a rejection of the null hypothesis (1) or not (0).
#'
#' \code{weibulltrialsimulatornohist()} should not be called directly by user.
#'
#' @param sample_size_val Number of subjects per arm.
#' @param scale1_val Scale parameter used in call to \code{rweibull()}.
#'   Used only in randomized control arm.
#' @param hazard_ratio_val Desired hazard ratio between experimental and control arms.
#' @param common_shape_val Shape parameter used in call to \code{rweibull()}.
#'   Used in both randomized arms.
#' @param censor_value Value at which time-to-event data are right censored.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#'
#' @return \code{weibulltrialsimulatornohist()} returns a vector of simulation results. The
#'   first element is an estimated hazard ratio, the second element is the estimated
#'   variance of the log hazard ratio, and the third element is a 0/1 variable indicator
#'   whether or not the trial rejected the null hypothesis (1) or failed to reject
#'   the null hypothesis (0).
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
weibulltrialsimulatornohist <- function(sample_size_val, scale1_val, hazard_ratio_val, common_shape_val, censor_value,
    alpha) {

  # --------------------------------------------------------------- #
  #  This function simulates a two-arm Bayesian trial where
  #  historical data is not utilized in the parameter estimation.
  # --------------------------------------------------------------- #

  # First, Generate weibull trial data given the user defined trial characteristics.
    sampleranddata <- genweibulldata(sample_size = sample_size_val, scale1 = scale1_val, hazard_ratio = hazard_ratio_val,
        common_shape = common_shape_val, censor_value = censor_value)
    # Make sure the trial data has at least one not right censored observation.
    if (sum(sampleranddata$event_time == censor_value) == dim(sampleranddata)[1]) {
        stop("Simulated trial data must have at least one observation that is not right censored.")
    }

    #Unlike Bernoulli, Poisson, Gaussian, and Lognormal, I cannot use survreg
    #directly, because survreg does not use a parameterization where the hazard
    #ratio is a model parameter.
    # Generate initial values for your call to optim()
    initializemodel <- survival::survreg(survival::Surv(event_time, status) ~ treatment, dist = "weibull", data = sampleranddata)

    initialbeta0 <- -1 * initializemodel$coefficients[1]/initializemodel$scale
    initialbeta1 <- -1 * initializemodel$coefficients[2]/initializemodel$scale
    initialv <- 1/initializemodel$scale

    # Generate the Bayesian CLT based parameter estimates needed for inference on hazard ratio.
    fitmod <- stats::optim(c(initialbeta0, initialbeta1, initialv), weibullloglikenohist, randdata = sampleranddata, method = "Nelder-Mead",
        hessian = TRUE)

    #Extract model parameters and statistics
    modparm <- fitmod$par
    covarmat <- solve(fitmod$hessian)

    weibullloghazard_ratio <- modparm[2]

    weibullhazard_ratio <- exp(weibullloghazard_ratio)
    lower_weibullhazard_ratio <- exp(weibullloghazard_ratio - stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2]))
    upper_weibullhazard_ratio <- exp(weibullloghazard_ratio + stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2]))

    #Make a decision about the simulated trial, reject or fail to reject null hypothesis.
    reject <- ifelse(((lower_weibullhazard_ratio > 1) | (upper_weibullhazard_ratio < 1)), 1, 0)
    output <- c(weibullhazard_ratio, covarmat[2, 2], reject)

    #Return the hazard ratio, the estimated variance of the log hazard ratio, and the trial decision.
    names(output) <- c("hazard_ratio", "loghazard_ratio_var", "reject")

    return(output)

}


#' Repeated Two Arm Bayesian Clinical Trial Simulation with Historical Data and
#' Weibull Outcome.
#'
#' \code{weibull_sim()} function only used internally by \code{historic_sim()}
#' function to run a set of trial simulations involving historical
#' control data and a Weibull outcome.  User defined simulation parameters are
#' used to generate a set of trial scenarios.  Each scenario is simulated multiple
#' times and then means are taken to calculate estimates of power, hazard ratios,
#' and other user requested study summary statistics like variance of hazard
#' ratio, bias (on hazard ratio scale), and mse (on hazard ratio scale).
#' The number of repeated simulations is defined by the user.
#'
#' \code{weibull_sim()} should not be called directly by user.
#'
#' @param trial_reps Number of trials to replicate within each combination of
#'   a0_val, subj_per_arm, effect_vals, and rand_control_diff.  As the number
#'   of trials increases, the precision of the estimate will increase. Default is
#'   100.
#' @param subj_per_arm A vector of sample sizes, all of which must be positive
#'   integers.
#' @param a0_vals A vector of power prior parameters ranging from 0 to 1, where 0
#'   implies no information from historical data should be used, 1 implies all of
#'   the information from historical data should be used.  A value between 0 and 1
#'   implies that a proportion of the information from historical data will be used.
#' @param effect_vals A vector of hazard ratios between randomized arms (randomized
#'   experimental over control), all of which must be positive.
#' @param rand_control_diff A vector of hazard ratios (randomized controls over
#'   historical controls) representing differences between historical and randomized
#'   controls.
#' @param hist_control_data A dataset of historical data.  Default is \code{NULL}.
#'   Historical datasets must have 4 columns: id, treatment, event_time, and
#'   status.  The value of treatment should be 0.  The values of event_time must
#'   be positive.  The values of status must be 0 (right censored event) or
#'   1 (observed event).
#' @param censor_value A single value at which right censoring occurs when
#'   simulating randomized subject outcomes.  Default is \code{NULL}, where
#'   \code{NULL} implies no right censoring.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#' @param get_var A TRUE/FALSE indicator of whether an array of variance
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param get_bias A TRUE/FALSE indicator of whether an array of bias
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param get_mse A TRUE/FALSE indicator of whether an array of MSE
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param quietly A TRUE/FALSE indicator of whether notes are printed
#'   to output about simulation progress as the simulation runs.  If
#'   running interactively in RStudio or running in the R console,
#'   \code{quietly} can be set to FALSE.  If running in a Notebook or
#'   knitr document, \code{quietly} needs to be set to TRUE.  Otherwise
#'   each note will be printed on a separate line and it will take up
#'   a lot of output space.  Default is \code{TRUE}.
#'
#' @return \code{weibull_sim()} returns an S3 object of class bayes_ctd_array.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
weibull_sim <- function(trial_reps=100, subj_per_arm, a0_vals, effect_vals,
                        rand_control_diff, hist_control_data, censor_value,
                        alpha=0.05, get_var=FALSE, get_bias=FALSE,
                        get_mse=FALSE, quietly=TRUE) {

    # --------------------------------------------------------------- #
    # For a set of user specified scenarios (defined by combinations
    # of user specified parameters), simulate "trial_reps" trials
    # and estimate power, hazard ratio estimate, and if requested by user:
    # variance of hazard ratio, bias , and mse.  Using a Weibull oucome
    # and incorporating data from historical controls.
    # --------------------------------------------------------------- #

    # Need to take the historical data and generate distributional parameter estimates
    histdata = hist_control_data
    hist_model <- survival::survreg(survival::Surv(event_time, status) ~ 1, dist = "weibull", data = histdata)
    bparm_histc <- exp(hist_model$coefficients)
    aparm_histc <- 1/hist_model$scale


    # Initialize arrays to hold power, var, mse, and bias estimate results as requested.
    len_val <- length(rand_control_diff) * length(effect_vals) * length(a0_vals) * length(subj_per_arm)
    power_results <- array(rep(0, len_val), c(length(subj_per_arm), length(a0_vals), length(effect_vals), length(rand_control_diff)))
    est_results <- array(rep(0, len_val), c(length(subj_per_arm), length(a0_vals), length(effect_vals), length(rand_control_diff)))
    if (get_mse == TRUE) {
        mse_results <- array(rep(0, len_val), c(length(subj_per_arm), length(a0_vals), length(effect_vals), length(rand_control_diff)))
    }
    if (get_bias == TRUE) {
        bias_results <- array(rep(0, len_val), c(length(subj_per_arm), length(a0_vals), length(effect_vals), length(rand_control_diff)))
    }
    if (get_var == TRUE) {
        var_results <- array(rep(0, len_val), c(length(subj_per_arm), length(a0_vals), length(effect_vals), length(rand_control_diff)))
    }

    # Cycle through the passed values in rand_control_diff, effect_vals, a0_val, and subj_per_arm to generate the
    # requested trial characteristics.
    for (diffs in 1:length(rand_control_diff)) {
        # Need to adjust the randomized control scale parameter given the historical control scale parameters and the hazard
        # ratios given in rand_control_diff
        bparm_randc <- exp(log(bparm_histc) - (1/aparm_histc) * log(rand_control_diff[diffs]))

        for (effvals in 1:length(effect_vals)) {
            for (a0vals in 1:length(a0_vals)) {
                for (sizes in 1:length(subj_per_arm)) {
                  if (!quietly){
                    cat("\r", c(subj_per_arm[sizes], a0_vals[a0vals], effect_vals[effvals], rand_control_diff[diffs]))
                  }
                  # For each combination of rand_control_diff, effect_vals, a0_vals, and subj_per_arm, simulate the trial
                  #trial_reps times and then calculate the mean reject rate to estimate power.  For bias, work on the
                  #hazard ratio scale and take the mean of all differences between estimated hazard ratios and the
                  #hazard ratio.  For mse, calculate the mean of squared differences between the
                  #estimated hazard ratios and the true hazard ratio value.
                  collect <- matrix(rep(0, 3 * trial_reps), ncol = 3)
                  for (k in 1:trial_reps) {
                    # sample_size_val will be equal to both arms
                    collect[k, ] <- weibulltrialsimulator(sample_size_val = subj_per_arm[sizes], histdata, scale1_val = bparm_randc,
                      hazard_ratio_val = effect_vals[effvals], common_shape_val = aparm_histc, censor_value = censor_value,
                      a0_val = a0_vals[a0vals], alpha = alpha)
                  }
                  #collect is a matrix of data, hazard ratio in 1st column, log hazard ratio variance
                  # in second column, and a vector of 0/1s in third column indicating whether or
                  # not trial represented by row led to a rejection of null hypothesis (1) or not (0).
                  # Note that collect gets rewritten for each scenario.
                  colnames(collect) <- c("hazard_ratio", "log_hazard_ratio_var", "reject")
                  #Start calculating means for each scenarios and placing the means in the proper
                  # array.  Every simulation will contain an array of power results and hazard
                  # ratio estimates.
                  power_results[sizes, a0vals, effvals, diffs] <- mean(collect[, 3])
                  est_results[sizes, a0vals, effvals, diffs] <- mean(collect[, 1])
                  if (get_bias == TRUE) {
                    bias_results[sizes, a0vals, effvals, diffs] <- mean(collect[, 1] - effect_vals[effvals])
                  }
                  if (get_var == TRUE) {
                    var_results[sizes, a0vals, effvals, diffs] <- mean((collect[, 1]*sqrt(collect[, 2]))^2)
                  }
                  if (get_mse == TRUE) {
                    mse_results[sizes, a0vals, effvals, diffs] <- mean((collect[, 1] - effect_vals[effvals])^2)
                  }
                  if (!quietly){
                    cat("\r", "                                                                                    ")
                  }
                }
            }
        }
    }
    cat("\n")

    #Lines 576 through 879 simply apply names to the dimensions of array created by the
    # simulation depending on values get_bias, get_var, and get_mse.
    if (get_bias == FALSE & get_var == FALSE & get_mse == FALSE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results)
        names(output) <- c("power", "est")
    }
    if (get_bias == FALSE & get_var == FALSE & get_mse == TRUE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(mse_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(mse_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(mse_results)[[2]] <- as.character(a0_vals)
        dimnames(mse_results)[[3]] <- as.character(effect_vals)
        dimnames(mse_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results, mse_results)
        names(output) <- c("power", "est", "mse")
    }
    if (get_bias == TRUE & get_var == FALSE & get_mse == FALSE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(bias_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(bias_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(bias_results)[[2]] <- as.character(a0_vals)
        dimnames(bias_results)[[3]] <- as.character(effect_vals)
        dimnames(bias_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results, bias_results)
        names(output) <- c("power", "est", "bias")
    }
    if (get_bias == TRUE & get_var == FALSE & get_mse == TRUE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(bias_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(bias_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(bias_results)[[2]] <- as.character(a0_vals)
        dimnames(bias_results)[[3]] <- as.character(effect_vals)
        dimnames(bias_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(mse_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(mse_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(mse_results)[[2]] <- as.character(a0_vals)
        dimnames(mse_results)[[3]] <- as.character(effect_vals)
        dimnames(mse_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results, bias_results, mse_results)
        names(output) <- c("power", "est", "bias", "mse")
    }
    if (get_bias == FALSE & get_var == TRUE & get_mse == FALSE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(var_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(var_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(var_results)[[2]] <- as.character(a0_vals)
        dimnames(var_results)[[3]] <- as.character(effect_vals)
        dimnames(var_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results, var_results)
        names(output) <- c("power", "est", "var")
    }
    if (get_bias == FALSE & get_var == TRUE & get_mse == TRUE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(var_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(var_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(var_results)[[2]] <- as.character(a0_vals)
        dimnames(var_results)[[3]] <- as.character(effect_vals)
        dimnames(var_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(mse_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(mse_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(mse_results)[[2]] <- as.character(a0_vals)
        dimnames(mse_results)[[3]] <- as.character(effect_vals)
        dimnames(mse_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results, var_results, mse_results)
        names(output) <- c("power", "est", "var", "mse")
    }
    if (get_bias == TRUE & get_var == TRUE & get_mse == FALSE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(var_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(var_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(var_results)[[2]] <- as.character(a0_vals)
        dimnames(var_results)[[3]] <- as.character(effect_vals)
        dimnames(var_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(bias_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(bias_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(bias_results)[[2]] <- as.character(a0_vals)
        dimnames(bias_results)[[3]] <- as.character(effect_vals)
        dimnames(bias_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results, var_results, bias_results)
        names(output) <- c("power", "est", "var", "bias")
    }
    if (get_bias == TRUE & get_var == TRUE & get_mse == TRUE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(bias_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(bias_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(bias_results)[[2]] <- as.character(a0_vals)
        dimnames(bias_results)[[3]] <- as.character(effect_vals)
        dimnames(bias_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(var_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(var_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(var_results)[[2]] <- as.character(a0_vals)
        dimnames(var_results)[[3]] <- as.character(effect_vals)
        dimnames(var_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(mse_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(mse_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(mse_results)[[2]] <- as.character(a0_vals)
        dimnames(mse_results)[[3]] <- as.character(effect_vals)
        dimnames(mse_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results, var_results, bias_results, mse_results)
        names(output) <- c("power", "est", "var", "bias", "mse")
    }

    #Create an list of results and apply the bayes_ctd_array class to the list, then
    # return the output object.
    class_out <- list(data = output, subj_per_arm = subj_per_arm, a0_vals = a0_vals, effect_vals = effect_vals, rand_control_diff = rand_control_diff, objtype = 'historic')
    class(class_out) <- append("bayes_ctd_array", class(class_out))
    return(class_out)

}

#' Two Arm Bayesian Clinical Trial Simulation with no Historical Data and
#' Weibull Outcome.
#'
#' \code{simple_weibull_sim()} function only used internally by
#' \code{simple_sim()} function to run a set of trial simulations involving no
#' historical control data and a Weibull outcome.  User defined simulation
#' parameters are used to generate a set of trial scenarios.  Each scenario is
#' simulated multiple times and then means are taken to calculate estimates
#' of power, hazard ratios, and other user requested study summary statistics
#' like variance of hazard ratio, bias (on hazard ratio scale), and
#' mse (on hazard ratio scale).  The number of repeated simulations is
#' defined by the user.
#'
#' \code{simple_weibull_sim()} should not be called directly by user.
#'
#' @param trial_reps Number of trials to replicate within each combination of
#'   subj_per_arm and effect_vals.  As the number of trials increases, the
#'   precision of the estimate will increase. Default is 100.
#' @param subj_per_arm A vector of sample sizes, all of which must be positive
#'   integers.
#' @param effect_vals A vector of hazard ratios (randomized experimental over control),
#'   all of which must be positive.
#' @param scale1_value scale parameter value for randomized controls. Used in call
#'   to \code{rweibull()}.
#' @param common_shape_value shape parameter value assumed common in both arms.
#'   Used in call to \code{rweibull()}.
#' @param censor_value A single value at which right censoring occurs when
#'   simulating randomized subject outcomes.  Default is \code{NULL}, where
#'   \code{NULL} implies no right censoring.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#' @param get_var A TRUE/FALSE indicator of whether an array of variance
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param get_bias A TRUE/FALSE indicator of whether an array of bias
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param get_mse A TRUE/FALSE indicator of whether an array of MSE
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param quietly A TRUE/FALSE indicator of whether notes are printed
#'   to output about simulation progress as the simulation runs.  If
#'   running interactively in RStudio or running in the R console,
#'   \code{quietly} can be set to FALSE.  If running in a Notebook or
#'   knitr document, \code{quietly} needs to be set to TRUE.  Otherwise
#'   each note will be printed on a separate line and it will take up
#'   a lot of output space.  Default is \code{TRUE}.
#'
#' @return \code{simple_weibull_sim()} returns an S3 object of class bayes_ctd_array.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
simple_weibull_sim <- function(trial_reps=100, subj_per_arm, effect_vals, scale1_value,
                               common_shape_value, censor_value, alpha=0.05,
                               get_var=FALSE, get_bias=FALSE, get_mse=FALSE,
                               quietly=TRUE) {

  # --------------------------------------------------------------- #
  # For a set of user specified scenarios (defined by combinations
  # of user specified parameters), simulate "trial_reps" trials
  # and estimate power, hazard ratio estimate, and if requested by user:
  # variance of hazard ratio, bias, and mse.  Using a Weibull oucome
  # but historical control data is not used.
  # --------------------------------------------------------------- #

    #The rand_control_diff and a0_val dimensions will be set to 1, and the value for
    # rand_control_diff will be 1 and a0_val will be set to 0.  All summaries will
    # be set up to ignore these dimensions for simple (no historical data) simulations.
    rand_control_diff <- 1
    a0_vals <- 0
    # Initialize arrays to hold power, var, mse, and bias estimate results as requested.
    len_val <- length(rand_control_diff) * length(effect_vals) * length(a0_vals) * length(subj_per_arm)
    power_results <- array(rep(0, len_val), c(length(subj_per_arm), length(a0_vals), length(effect_vals), length(rand_control_diff)))
    est_results <- array(rep(0, len_val), c(length(subj_per_arm), length(a0_vals), length(effect_vals), length(rand_control_diff)))
    if (get_mse == TRUE) {
        mse_results <- array(rep(0, len_val), c(length(subj_per_arm), length(a0_vals), length(effect_vals), length(rand_control_diff)))
    }
    if (get_bias == TRUE) {
        bias_results <- array(rep(0, len_val), c(length(subj_per_arm), length(a0_vals), length(effect_vals), length(rand_control_diff)))
    }
    if (get_var == TRUE) {
        var_results <- array(rep(0, len_val), c(length(subj_per_arm), length(a0_vals), length(effect_vals), length(rand_control_diff)))
    }

    # Cycle through the passed values in rand_control_diff, effect_vals, a0_val, and subj_per_arm to generate the
    # requested trial characteristics.  Note that rand_control_diff is set to 1 and a0_val is set to 0.
    for (diffs in 1:length(rand_control_diff)) {
        for (effvals in 1:length(effect_vals)) {
            for (a0vals in 1:length(a0_vals)) {
                for (sizes in 1:length(subj_per_arm)) {
                  if (!quietly){
                    cat("\r", c(subj_per_arm[sizes], a0_vals[a0vals], effect_vals[effvals], rand_control_diff[diffs]))
                  }
                  # For each combination of rand_control_diff, effect_vals, a0_val, and subj_per_arm, simulate the trial
                  # trial_reps times and then calculate the mean reject rate to estimate power.  For bias, work on the
                  #hazard ratio scale and take the mean of all differences between estimated hazard ratios and the
                  #true hazard ratio.  For mse, calculate the mean of squared differences between the estimated
                  #hazard ratios and the true hazard ratio value.  Note that rand_control_diff is set to 1 and
                  #a0_val is set to 0.
                  collect <- matrix(rep(0, 3 * trial_reps), ncol = 3)
                  for (k in 1:trial_reps) {
                    # sample_size_val will be equal in both arms
                    collect[k, ] <- weibulltrialsimulatornohist(sample_size_val = subj_per_arm[sizes], scale1_val = scale1_value,
                      hazard_ratio_val = effect_vals[effvals], common_shape_val = common_shape_value, censor_value = censor_value,
                      alpha = alpha)
                  }
                  #collect is a matrix of data, hazard ratio in 1st column, log hazard ratio variance
                  # in second column, and a vector of 0/1s in third column indicating whether or
                  # not trial represented by row led to a rejection of null hypothesis (1) or not (0).
                  # Note that collect gets rewritten for each scenario.
                  colnames(collect) <- c("hazard_ratio", "log_hazard_ratio_var", "reject")
                  #Start calculating means for each scenarios and placing the means in the proper
                  # array.  Every simulation will contain an array of power results and hazard
                  # ratio estimates.
                  power_results[sizes, a0vals, effvals, diffs] <- mean(collect[, 3])
                  est_results[sizes, a0vals, effvals, diffs] <- mean(collect[, 1])
                  if (get_bias == TRUE) {
                    bias_results[sizes, a0vals, effvals, diffs] <- mean(collect[, 1] - effect_vals[effvals])
                  }
                  if (get_var == TRUE) {
                    var_results[sizes, a0vals, effvals, diffs] <- mean((collect[, 1]*sqrt(collect[, 2]))^2)
                  }
                  if (get_mse == TRUE) {
                    mse_results[sizes, a0vals, effvals, diffs] <- mean((collect[, 1] - effect_vals[effvals])^2)
                  }
                  if (!quietly){
                    cat("\r", "                                                                                    ")
                  }
                }
            }
        }
    }
    cat("\n")

    #Lines 1011 through 1314 simply apply names to the dimensions of array created by the
    # simulation depending on values get_bias, get_var, and get_mse.
    if (get_bias == FALSE & get_var == FALSE & get_mse == FALSE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results)
        names(output) <- c("power", "est")
    }
    if (get_bias == FALSE & get_var == FALSE & get_mse == TRUE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(mse_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(mse_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(mse_results)[[2]] <- as.character(a0_vals)
        dimnames(mse_results)[[3]] <- as.character(effect_vals)
        dimnames(mse_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results, mse_results)
        names(output) <- c("power", "est", "mse")
    }
    if (get_bias == TRUE & get_var == FALSE & get_mse == FALSE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(bias_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(bias_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(bias_results)[[2]] <- as.character(a0_vals)
        dimnames(bias_results)[[3]] <- as.character(effect_vals)
        dimnames(bias_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results, bias_results)
        names(output) <- c("power", "est", "bias")
    }
    if (get_bias == TRUE & get_var == FALSE & get_mse == TRUE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(bias_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(bias_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(bias_results)[[2]] <- as.character(a0_vals)
        dimnames(bias_results)[[3]] <- as.character(effect_vals)
        dimnames(bias_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(mse_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(mse_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(mse_results)[[2]] <- as.character(a0_vals)
        dimnames(mse_results)[[3]] <- as.character(effect_vals)
        dimnames(mse_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results, bias_results, mse_results)
        names(output) <- c("power", "est", "bias", "mse")
    }
    if (get_bias == FALSE & get_var == TRUE & get_mse == FALSE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(var_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(var_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(var_results)[[2]] <- as.character(a0_vals)
        dimnames(var_results)[[3]] <- as.character(effect_vals)
        dimnames(var_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results, var_results)
        names(output) <- c("power", "est", "var")
    }
    if (get_bias == FALSE & get_var == TRUE & get_mse == TRUE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(var_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(var_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(var_results)[[2]] <- as.character(a0_vals)
        dimnames(var_results)[[3]] <- as.character(effect_vals)
        dimnames(var_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(mse_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(mse_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(mse_results)[[2]] <- as.character(a0_vals)
        dimnames(mse_results)[[3]] <- as.character(effect_vals)
        dimnames(mse_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results, var_results, mse_results)
        names(output) <- c("power", "est", "var", "mse")
    }
    if (get_bias == TRUE & get_var == TRUE & get_mse == FALSE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(var_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(var_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(var_results)[[2]] <- as.character(a0_vals)
        dimnames(var_results)[[3]] <- as.character(effect_vals)
        dimnames(var_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(bias_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(bias_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(bias_results)[[2]] <- as.character(a0_vals)
        dimnames(bias_results)[[3]] <- as.character(effect_vals)
        dimnames(bias_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results, var_results, bias_results)
        names(output) <- c("power", "est", "var", "bias")
    }
    if (get_bias == TRUE & get_var == TRUE & get_mse == TRUE) {
        if (length(subj_per_arm) == 1) {
            dimnames(power_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(power_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(power_results)[[2]] <- as.character(a0_vals)
        dimnames(power_results)[[3]] <- as.character(effect_vals)
        dimnames(power_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(est_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(est_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(est_results)[[2]] <- as.character(a0_vals)
        dimnames(est_results)[[3]] <- as.character(effect_vals)
        dimnames(est_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(bias_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(bias_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(bias_results)[[2]] <- as.character(a0_vals)
        dimnames(bias_results)[[3]] <- as.character(effect_vals)
        dimnames(bias_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(var_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(var_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(var_results)[[2]] <- as.character(a0_vals)
        dimnames(var_results)[[3]] <- as.character(effect_vals)
        dimnames(var_results)[[4]] <- as.character(rand_control_diff)

        if (length(subj_per_arm) == 1) {
            dimnames(mse_results)[[1]] <- list(as.character(subj_per_arm))
        }
        if (length(subj_per_arm) > 1) {
            dimnames(mse_results)[[1]] <- as.character(subj_per_arm)
        }
        dimnames(mse_results)[[2]] <- as.character(a0_vals)
        dimnames(mse_results)[[3]] <- as.character(effect_vals)
        dimnames(mse_results)[[4]] <- as.character(rand_control_diff)
        output <- list(power_results, est_results, var_results, bias_results, mse_results)
        names(output) <- c("power", "est", "var", "bias", "mse")
    }

    #Create an list of results and apply the bayes_ctd_array class to the list, then
    # return the output object.
    class_out <- list(data = output, subj_per_arm = subj_per_arm, a0_vals = 0, effect_vals = effect_vals, rand_control_diff = 1, objtype = 'simple')
    class(class_out) <- append("bayes_ctd_array", class(class_out))
    return(class_out)

}

