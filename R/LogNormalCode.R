#' Generating function for Lognormal Data.
#'
#' \code{genlognormaldata()} function used mainly internally by
#' \code{lognormaltrialsimulator()} and \code{lognormaltrialsimulatornohist()} functions
#' to generate data for a two-arm clinical trial, experimental and control groups.
#' Can be used to generate random trial data.
#'
#' @param sample_size  Number of subjects per arm.
#' @param mu1 meanlog parameter used in call to \code{rlnorm()}.
#'   Used only in control arm.
#' @param mean_ratio Desired Mean Ratio between experimental and control groups.
#' @param common_sd sdlog parameter used in call to \code{rlnorm()}.
#'   Used in both arms.
#' @param censor_value Value at which time-to-event data are right censored.
#'
#' @return \code{genlognormaldata()} returns a data frame with columns: 'id', 'treatment',
#'   'event_time', and 'status'.
#'
#' @examples
#' samplehistdata <- genlognormaldata(sample_size=60, mu1=1.06, mean_ratio=0.6,
#'                                    common_sd=1.25, censor_value=3)
#' samplehistdata
#' @export
genlognormaldata <- function(sample_size, mu1, mean_ratio, common_sd, censor_value) {

  # --------------------------------------------------------------- #
	# The function genlognormaldata simulates a balanced clinical trial
	# with 'sample_size' subjects per arm using a lognormal distribution.
  # 'm1' is the lognormal location parameter, and 'common_sd' is the
	# lognormal scale parameter for both arms.  'censor_value' is the
	# value when right censoring occurs.  As of 9/5/2016,
	# genlognormaldata only generates data with right censoring.
	# Random right censoring is not incorporated.  'mean_ratio is the
	# ratio of group means (experimental group over control group).
	#
	# In the code below time1, mu1, test1, status1, etc. are data
	# for the control goup.
	# In the code below time2, mu2, test2, status2, etc. are data
	# for the experimental group.
	# --------------------------------------------------------------- #

    # mu1 is the lognormal distribution mu parameter for the control group.
	  # The lognormal mean for control group is:
    lnmean1 <- exp(mu1 + 0.5 * common_sd^2)
    # given user specified mean_ratio, the lognormal mean for experimental
	  # group is:
    lnmean2 <- mean_ratio * lnmean1
    # Now, I need to transform the lognormal mean in the experimental group
	  # to the mu parameter for experimental group.
    mu2 <- log(lnmean2) - 0.5 * common_sd^2

	  # Create event times for both groups
    time1 <- stats::rlnorm(sample_size, meanlog = mu1, sdlog = common_sd)
    time2 <- stats::rlnorm(sample_size, meanlog = mu2, sdlog = common_sd)

    # Create variables needed for simulation when censor_value is specified.
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

#' Log-likelihood function for two-arm trial with historical data using Lognormal
#' distribution.
#'
#' \code{lognormalloglike()} function only used internally by
#' \code{lognormaltrialsimulator()} function to estimate Lognormal model parameters
#' when clinical trial involves experimental and control groups as well as historical
#' control data.  The lognormal log-likelihood is calculated by modeling \code{log(data)}
#' as a Gaussian random variable. Not to be called directly by user.
#'
#' @param params  Three element vector of Lognormal parameters.  Third element is log(sd),
#'   where sd is a parameter required by dnorm().  The first and second elements
#'   are the intercept (beta0) and treatment effect parameter (beta1), where the treatment effect is
#'   a log mean ratio (experimental group over control group).  The mu parameter required by
#'   dnorm() is equal to params[1] + params[2]*treatment.  It is assumed that the log(sd)
#'   parameter is the same in both randomized and historical data.  It is assumed that
#'   the mu parameter in the randomized and historical control data is equal to params[1].
#' @param randdata  Dataset of randomly generated trial data.  Randomized trial datasets
#'   must have 4 columns: id, treatment, event_time, and status.  The value of treatment
#'   must be 0 (control) or 1 (experimental).  The values of event_time must be positive.
#'   The values of status must be 0 (right censored event) or 1 (observed event).
#' @param histdata Dataset of historical data.  Historical datasets must have 4 columns:
#'   id, treatment, event_time, and status.  The value of treatment should be 0.  The
#'   values of event_time must be positive.  The values of status must be 0 (right
#'   censored event) or 1 (observed event).
#' @param a0 Power prior parameter: 0 implies historical data is ignored and 1 implies
#'   all information in historical data is used.
#'
#' @return \code{lognormalloglike()} returns a value of the loglikelihood function
#'   given a set of Lognormal parameters, randomly generated trial data, and observed
#'   historical data.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
lognormalloglike <- function(params, randdata, histdata, a0) {
  # --------------------------------------------------------------- #
	#  This function calculates the Lognormal log-likelihood given
	#  a vector of parameter values, a dataset of randomized trial
	#  data (two arms, no covariates beyond treatment), and a dataset
	#  of historical control data.
	#  The lognormal scale parameter is common in both randomized
	#  groups and the historical control group.  The lognormal location
	#  parameter is assumed to be the same in both control groups.
	#  The location parameter for the randomized experimental group is
	#  a linear function of the control location and the treatment
	#  effect.  The parameters are beta0, beta1, and logs.  logs is the
	#  log of the common scale parameter.  beta0 and beta1 are regression
	#  parameters that are linked to the Lognoraml location parameter,
	#  mu, via the identity link.  beta1 is the log mean ratio
	#  (experimental group over control group), while beta0 is the
	#  location parameter for controls.
  # --------------------------------------------------------------- #

    # Get params.  Note that the scale parameter is on the log scale.
	  # To calculate the lognormal standard deviation value one must use exp(logs).
    beta0 <- params[1]
    beta1 <- params[2]
    logs <- params[3]

    # Calculate the lognormal mu parameter vector for all randomized observations.
    mu_i <- beta0 + beta1 * randdata$treatment
	  # Calculate the log-likelihood values for all randomized observations.
	  # Note that I am modeling the lognormal by taking the log of event_time and then
	  # using the Gaussian distribution for the log transformed values.  Could have used
	  # lnorm(), but did it this way just because it is what I am used to doing.
    ll_R <- randdata$status * stats::dnorm(log(randdata$event_time), mean = mu_i, sd = exp(logs), log = TRUE) + (1 - randdata$status) *
      stats::pnorm(log(randdata$event_time), mean = mu_i, sd = exp(logs), log.p = TRUE, lower.tail = FALSE)

	  # Calculate the loglikelihood values for all historical control observations.
	  # Note that the lognormal mu parameter is equal to beta0.
    ll_H <- histdata$status * stats::dnorm(log(histdata$event_time), mean = beta0, sd = exp(logs), log = TRUE) + (1 - histdata$status) *
      stats::pnorm(log(histdata$event_time), mean = beta0, sd = exp(logs), log.p = TRUE, lower.tail = FALSE)

	  # Calculate the overall log likelihood by adding the randomized log-likelihood to the historical control
	  # log-likelihood by a0, where a0 is the power prior parameter.  This a0 value is defined by the
	  # user and not estimated via object function optimization.
    ll <- sum(ll_R) + a0 * sum(ll_H)

	  # Return the sum of all individual elements to the negative log-likelihood
    return(-ll)
}

#' Log-likelihood function for two-arm trial with no historical data using Lognormal distribution.
#'
#' \code{lognormalloglikenohist()} function only used internally by
#' \code{lognormaltrialsimulatornohist()} function to estimate Lognormal model parameters
#' when clinical trial involves experimental and control groups but no historical control
#' data.  The lognormal log-likelihood is calculated by modeling \code{log(data)} as
#' a Gaussian random variable. Not to be called directly by user.
#'
#' @param params  Three element vector of Lognormal parameters.  Third element is log(sd),
#'   where sd is a parameter required by dnorm().  The first and second elements
#'   are the intercept (beta0) and treatment effect parameter (beta1), where the treatment effect is
#'   a log mean ratio (experimental group over control group).  The mu parameter required by
#'   dnorm() is equal to params[1] + params[2]*treatment.  It is assumed that
#'   the mu parameter in the randomized control data is equal to params[1].
#' @param randdata  Dataset of randomly generated trial data.  Randomized trial datasets
#'   must have 4 columns: id, treatment, event_time, and status.  The value of treatment
#'   must be 0 (control) or 1 (experimental).  The values of event_time must be positive.
#'   The values of status must be 0 (right censored event) or 1 (observed event).
#'
#' @return \code{lognormalloglikenohist()} returns a value of the loglikelihood function
#'   given a set of Lognormal parameters and randomly generated trial data.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
lognormalloglikenohist <- function(params, randdata) {
    # --------------------------------------------------------------- #
	#  This function calculates the Lognormal log-likelihood given
	#  a vector of parameter values, and a dataset of randomized trial
	#  data (two arms, no covariates beyond treatment).  Historical
    #  data is not utilized in this log-likelihood function.
	#  The lognormal scale parameter is common in both randomized
	#  groups.  The location parameter for the randomized experimental
	#  group is a linear function of the control location and the
	#  treatment effect.  The parameters are beta0, beta1, and logs.
	#  logs is the log of the common scale parameter.  beta0 and beta1
	#  are regression parameters that are linked to the Lognoraml
	#  location parameter, mu, via the identity link.  beta1 is the
	#  log mean ratio (experimental group over control group), while
	#  beta0 is the location parameter for controls.
    # --------------------------------------------------------------- #

    # Get params.  Note that the scale parameter is on the log scale.
	# To calculate the lognormal standard deviation value one must use exp(logs).
    beta0 <- params[1]
    beta1 <- params[2]
    logs <- params[3]

    # Calculate the lognormal mu parameter vector for all randomized observations.
    mu_i <- beta0 + beta1 * randdata$treatment
	# Calculate the log-likelihood values for all randomized observations.
	# Note that I am modeling the lognormal by taking the log of event_time and then
	# using the Gaussian distribution for the log transformed values.  Could have used
	# lnorm(), but did it this way just because it is what I am used to doing.
    ll_R <- randdata$status * stats::dnorm(log(randdata$event_time), mean = mu_i, sd = exp(logs), log = TRUE) + (1 - randdata$status) *
      stats::pnorm(log(randdata$event_time), mean = mu_i, sd = exp(logs), log.p = TRUE, lower.tail = FALSE)

	# Return the sum of all individual elements to the negative log-likelihood
    return(sum(-ll_R))
}


#' Simulate a single randomized trial using a Lognormal outcome and information from
#' historical controls.
#'
#' \code{lognormaltrialsimulator()} function only used internally by
#' \code{lognormal_sim()} function to run a single trial simulation involving historical
#' control data and a Lognormal outcome.
#'
#' The simulation of a trial with a Lognormal outcome involving historical control data returns
#' an estimate of the mean ratio as well as an estimate of the log mean ratio variance.
#' Finally the simulation returns an indication of whether or not the simulated trial led to
#' a rejection of the null hypothesis (1) or not (0).
#'
#' \code{lognormaltrialsimulator()} should not be called directly by user.
#'
#' @param sample_size_val Number of subjects per arm.
#' @param histdata Dataset of historical data.  Historical datasets must have 4 columns:
#'   id, treatment, event_time, and status.  The value of treatment should be 0.  The
#'   values of event_time must be positive.  The values of status must be 0 (right
#'   censored event) or 1 (observed event).
#' @param mu1_val Randomized control arm meanlog parameter used in call to \code{rlnorm()}.
#' @param mean_ratio_val Desired mean ratio between randomized experimental and control arms.
#' @param common_sd_val Randomized sdlog parameter used in call to \code{rlnorm()}.
#'   Used in both randomized arms.
#' @param censor_value Value at which time-to-event data are right censored.
#' @param a0_val A power prior parameter ranging from 0 to 1, where 0
#'   implies no information from historical data should be used, 1 implies all of
#'   the information from historical data should be used.  A value between 0 and 1
#'   implies that a proportion of the information from historical data will be used.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#'
#' @return \code{lognormaltrialsimulator()} returns a vector of simulation results. The
#'   first element is an estimated mean ratio, the second element is the estimated
#'   variance of the log mean ratio, and the third element is a 0/1 variable indicator
#'   whether or not the trial rejected the null hypothesis (1) or failed to reject
#'   the null hypothesis (0).
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
lognormaltrialsimulator <- function(sample_size_val, histdata, mu1_val, mean_ratio_val, common_sd_val, censor_value,
    a0_val, alpha) {
    # --------------------------------------------------------------- #
    #  This function simulates a two-arm Bayesian trial where
    #  historical data is utilized in the parameter estimation.
    # --------------------------------------------------------------- #

    # First, Generate lognormal trial data given the user defined trial characteristics.
    sampleranddata <- genlognormaldata(sample_size = sample_size_val, mu1 = mu1_val, mean_ratio = mean_ratio_val, common_sd = common_sd_val,
        censor_value = censor_value)
	# Make sure the trial data has at least one not right censored observation.
    if (sum(sampleranddata$event_time == censor_value) == dim(sampleranddata)[1]) {
        stop("Simulated trial data must have at least one observation that is not right censored.")
    }
    # Generate initial values for your call to optim()
    initializemodel <- survival::survreg(survival::Surv(event_time, status) ~ treatment, dist = "lognormal", data = sampleranddata)

    initialbeta0 <- initializemodel$coefficients[1]
    initialbeta1 <- initializemodel$coefficients[2]
    initiallogs <- log(initializemodel$scale)

	# Generate the Bayesian CLT based parameter estimates needed for inference on mean ratio.
    fitmod <- stats::optim(c(initialbeta0, initialbeta1, initiallogs), lognormalloglike, randdata = sampleranddata, histdata = histdata,
        a0 = a0_val, method = "Nelder-Mead", hessian = TRUE)

    #Extract model parameters and statistics
    modparm <- fitmod$par
    covarmat <- solve(fitmod$hessian)

    lognormallogmeanratio <- modparm[2]

    lognormalmeanratio <- exp(lognormallogmeanratio)
    lower_lognormalmeanratio <- exp(lognormallogmeanratio - stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2]))
    upper_lognormalmeanratio <- exp(lognormallogmeanratio + stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2]))

	#Make a decision about the simulated trial, reject or fail to reject null hypothesis.
    reject <- ifelse(((lower_lognormalmeanratio > 1) | (upper_lognormalmeanratio < 1)), 1, 0)
    output <- c(lognormalmeanratio, covarmat[2, 2], reject)

	#Return the mean ratio, the estimated variance of the log mean ratio, and the trial decision.
    names(output) <- c("mean_ratio", "log_mean_ratio_var", "reject")
    return(output)

}

#' Simulate a single randomized trial using a Lognormal outcome but not including any information from
#' historical controls.
#'
#' \code{lognormaltrialsimulatornohist()} function only used internally by
#' \code{simple_lognormal_sim()} function to estimate Lognormal model parameters
#' when clinical trial involves experimental and control groups but no historical control
#' data.
#'
#' The simulation of a trial with a Lognormal outcome involving no historical control data returns
#' an estimate of the mean ratio as well as an estimate of the log mean ratio variance.
#' Finally the simulation returns an indication of whether or not the simulated trial led to
#' a rejection of the null hypothesis (1) or not (0).
#'
#' \code{lognormaltrialsimulatornohist()} should not be called directly by user.
#'
#' @param sample_size_val Number of subjects per arm.
#' @param mu1_val meanlog parameter used in call to \code{rlnorm()}.
#'   Used only in randomized control arm.
#' @param mean_ratio_val Desired mean ratio between experimental and control arms.
#' @param common_sd_val sdlog parameter used in call to \code{rlnorm()}.
#'   Used in both arms.
#' @param censor_value Value at which time-to-event data are right censored.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#'
#' @return \code{lognormaltrialsimulatornohist()} returns a vector of simulation results. The
#'   first element is an estimated mean ratio, the second element is the estimated
#'   variance of the log mean ratio, and the third element is a 0/1 variable indicator
#'   whether or not the trial rejected the null hypothesis (1) or failed to reject
#'   the null hypothesis (0).
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
lognormaltrialsimulatornohist <- function(sample_size_val, mu1_val, mean_ratio_val, common_sd_val, censor_value, alpha) {

  # --------------------------------------------------------------- #
  #  This function simulates a two-arm Bayesian trial where
  #  historical data is not utilized in the parameter estimation.
  # --------------------------------------------------------------- #

    # First, Generate lognormal trial data given the user defined trial characteristics.
    sampleranddata <- genlognormaldata(sample_size = sample_size_val, mu1 = mu1_val, mean_ratio = mean_ratio_val, common_sd = common_sd_val,
        censor_value = censor_value)
	# Make sure the trial data has at least one not right censored observation.
    if (sum(sampleranddata$event_time == censor_value) == dim(sampleranddata)[1]) {
        stop("Simulated trial data must have at least one observation that is not right censored.")
    }
    # Generate the Bayesian CLT based parameter estimates needed for inference on mean ratio.
    initializemodel <- survival::survreg(survival::Surv(event_time, status) ~ treatment, dist = "lognormal", data = sampleranddata)

    #Extract model parameters and statistics
    modparm <- initializemodel$coefficients
    covarmat <- stats::vcov(initializemodel)

    lognormallogmeanratio <- modparm[2]

    lognormalmeanratio <- exp(lognormallogmeanratio)
    lower_lognormalmeanratio <- exp(lognormallogmeanratio - stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2]))
    upper_lognormalmeanratio <- exp(lognormallogmeanratio + stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2]))

	#Make a decision about the simulated trial, reject or fail to reject null hypothesis.
    reject <- ifelse(((lower_lognormalmeanratio > 1) | (upper_lognormalmeanratio < 1)), 1, 0)
    output <- c(lognormalmeanratio, covarmat[2, 2], reject)

	#Return the hazard ratio, the estimated variance of the log hazard ratio, and the trial decision.
    names(output) <- c("mean_ratio", "log_mean_ratio_var", "reject")
    return(output)

}



#' Repeated Two Arm Bayesian Clinical Trial Simulation with Historical Data and
#' Lognormal Outcome.
#'
#' \code{lognormal_sim()} function only used internally by \code{historic_sim()}
#' function to run a set of trial simulations involving historical
#' control data and a Lognormal outcome.  User defined simulation parameters are
#' used to generate a set of trial scenarios.  Each scenario is simulated multiple
#' times and then means are taken to calculate estimates of power, mean ratios,
#' and other user requested study summary statistics like variance of mean
#' ratio, bias (on mean ratio scale), and mse (on mean ratio scale).
#' The number of repeated simulations is defined by the user.
#'
#' \code{lognormal_sim()} should not be called directly by user.
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
#' @param effect_vals A vector of mean ratios (randomized experimental over control),
#'   all of which must be positive.
#' @param rand_control_diff For Lognormal outcomes this is a vector of mean ratios
#'   (randomized controls over historical controls) that represent differences
#'   between randomized and historical controls.
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
#' @return \code{lognormal_sim()} returns an S3 object of class bayes_ctd_array.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
lognormal_sim <- function(trial_reps=100, subj_per_arm, a0_vals, effect_vals,
                          rand_control_diff, hist_control_data, censor_value,
                          alpha=0.05, get_var=FALSE, get_bias=FALSE, get_mse=FALSE,
                          quietly=TRUE) {

    # --------------------------------------------------------------- #
    # For a set of user specified scenarios (defined by combinations
    # of user specified parameters), simulate "trial_reps" trials
    # and estimate power, mean ratio estimate, and if requested by user:
    # variance of mean ratio, bias, and mse.  Using a Lognormal oucome
    # and incorporating data from historical controls.
    # --------------------------------------------------------------- #

    # Need to take the historical data and generate distributional parameter estimates
    histdata <- hist_control_data
    hist_model <- survival::survreg(survival::Surv(event_time, status) ~ 1, dist = "lognormal", data = histdata)
    initialmu1 <- hist_model$coefficients[1]
    initialsd1 <- hist_model$scale


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
        # Need to adjust the randomized control mean given the historical control mean and the mean ratios given in
        # rand_control_diff.  Note, mean ratios in rand_control_diff are ratios of lognormal means
        # and not mu parameters of lognormal distributions (means of variables on log scale), so
        # you first need to calculate the lognormal mean for the historical control group.

        # lognormal mean for historical control group is:
        lnmean1 <- exp(initialmu1 + 0.5 * initialsd1^2)
        # given mean_ratio, lognormal for randomized control group is:
        lnmean2 <- rand_control_diff[diffs] * lnmean1
        # Now that we have the randomized control lognormal mean, we need to transform this
        # back to a lognormal mu parameter:
        # lognormal mu parameter for randomized control group is:
        adjmu1 <- log(lnmean2) - 0.5 * initialsd1^2

        for (effvals in 1:length(effect_vals)) {
            for (a0vals in 1:length(a0_vals)) {
                for (sizes in 1:length(subj_per_arm)) {
                  if (!quietly){
                    cat("\r", c(subj_per_arm[sizes], a0_vals[a0vals], effect_vals[effvals], rand_control_diff[diffs]))
                  }
                  # For each combination of rand_control_diff, effect_vals, a0_vals, and subj_per_arm, simulate the trial
                  #trial_reps times and then calculate the mean reject rate to estimate power.  For bias, work on the
                  #mean ratio scale and take the mean of all differences between estimated mean ratios and the
                  #true mean ratio.  For mse, calculate the mean of squared differences between the
                  #estimated mean ratios and the true mean ratio value.
                  collect <- matrix(rep(0, 3 * trial_reps), ncol = 3)
                  for (k in 1:trial_reps) {
                    # sample_size_val will be equal to both arms
                    collect[k, ] <- lognormaltrialsimulator(sample_size_val = subj_per_arm[sizes], histdata, mu1_val = adjmu1,
                      mean_ratio_val = effect_vals[effvals], common_sd_val = initialsd1, censor_value = censor_value,
                      a0_val = a0_vals[a0vals], alpha = alpha)
                  }
                  #collect is a matrix of data, mean ratio in 1st column, log mean ratio variance
                  # in second column, and a vector of 0/1s in third column indicating whether or
                  # not trial represented by row led to a rejection of null hypothesis (1) or not (0).
                  # Note that collect gets rewritten for each scenario.
                  colnames(collect) <- c("mean_ratio", "log_mean_ratio_var", "reject")
                  #Start calculating means for each scenarios and placing the means in the proper
                  # array.  Every simulation will contain an array of power results and mean
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

    #Lines 534 through 837 simply apply names to the dimensions of array created by the
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
    class_out <- list(data = output, subj_per_arm = subj_per_arm, a0_vals = a0_vals, effect_vals = effect_vals, rand_control_diff = rand_control_diff, objtype= 'historic')
    class(class_out) <- append("bayes_ctd_array", class(class_out))
    return(class_out)

}


#' Repeated Two Arm Bayesian Clinical Trial Simulation with no Historical Data and
#' Lognormal Outcome.
#'
#' \code{simple_lognormal_sim()} function only used internally by \code{simple_sim()}
#' function to run a set of trial simulations involving no
#' historical control data and a Lognormal outcome.  User defined simulation
#' parameters are used to generate a set of trial scenarios.  Each scenario is
#' simulated multiple times and then means are taken to calculate estimates
#' of power, mean ratios, and other user requested study summary statistics
#' like variance of mean ratio, bias (on mean ratio scale), and
#' mse (on mean ratio scale).  The number of repeated simulations is
#' defined by the user.
#'
#' \code{simple_lognormal_sim()} should not be called directly by user.
#'
#' @param trial_reps Number of trials to replicate within each combination of
#'   subj_per_arm and effect_vals.  As the number of trials increases, the
#'   precision of the estimate will increase. Default is 100.
#' @param subj_per_arm A vector of sample sizes, all of which must be positive
#'   integers.
#' @param effect_vals A vector of mean ratios (randomized experimental over control),
#'   all of which must be positive.
#' @param mu1_val meanlog parameter value for randomized control arm. Used in call
#'   to \code{rlnorm()}.
#' @param common_sd_val sdlog parameter value used in both randomized arms. Used in call to
#'   \code{rlnorm()}.
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
#' @return \code{simple_lognormal_sim()} returns an S3 object of class bayes_ctd_array.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
simple_lognormal_sim <- function(trial_reps, subj_per_arm, effect_vals, mu1_val, common_sd_val, censor_value, alpha,
    get_var, get_bias, get_mse, quietly=TRUE) {

  # --------------------------------------------------------------- #
  # For a set of user specified scenarios (defined by combinations
  # of user specified parameters), simulate "trial_reps" trials
  # and estimate power, mean ratio estimate, and if requested by user:
  # variance of mean ratio, bias, and mse.  Using a Lognormal oucome
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
                  #mean ratio scale and take the mean of all differences between estimated mean ratios and the
                  #true mean ratio.  For mse, calculate the mean of squared differences between the estimated
                  #mean ratios and the true mean ratio value.  Note that rand_control_diff is set to 1 and
                  #a0_val is set to 0.
                  collect <- matrix(rep(0, 3 * trial_reps), ncol = 3)
                  for (k in 1:trial_reps) {
                    # sample_size_val will be equal to both arms
                    collect[k, ] <- lognormaltrialsimulatornohist(sample_size_val = subj_per_arm[sizes], mu1_val = mu1_val,
                      mean_ratio_val = effect_vals[effvals], common_sd_val = common_sd_val, censor_value = censor_value,
                      alpha = alpha)
                  }
                  #collect is a matrix of data, mean ratio in 1st column, log mean ratio variance
                  # in second column, and a vector of 0/1s in third column indicating whether or
                  # not trial represented by row led to a rejection of null hypothesis (1) or not (0).
                  # Note that collect gets rewritten for each scenario.
                  colnames(collect) <- c("mean_ratio", "log_mean_ratio_var", "reject")
                  #Start calculating means for each scenarios and placing the means in the proper
                  # array.  Every simulation will contain an array of power results and mean
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

    #Lines 970 through 1273 simply apply names to the dimensions of array created by the
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
    class_out <- list(data = output, subj_per_arm = subj_per_arm, a0_vals = 0, effect_vals = effect_vals, rand_control_diff = 1, objtype= 'simple')
    class(class_out) <- append("bayes_ctd_array", class(class_out))
    return(class_out)

}
