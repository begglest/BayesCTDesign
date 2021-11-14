#' Generating function for Bernoulli Data.
#'
#' \code{genlogisticdata()} function used mainly internally by
#' \code{logistictrialsimulator()} function to generate data for a two-arm
#' clinical trial, experimental and control groups.  Can be used to generate
#' random trial data.
#'
#' @param sample_size  Number of subjects per arm.
#' @param prob1 prob parameter used in call to \code{rbinom()}.
#'   Used only in control arm.
#' @param odds_ratio Desired Odds Ratio between experimental and control groups.
#'
#' @return \code{genlogisticdata()} returns a data frame with columns: 'id', 'treatment',
#'   and 'y'.
#'
#' @examples
#' samplehistdata <- genbernoullidata(sample_size=60, prob1=0.6, odds_ratio=0.6)
#' samplehistdata
#' @export
#'
genbernoullidata <- function(sample_size, prob1, odds_ratio) {

    # --------------------------------------------------------------- #
	# The function genbernoullidata simulates a balanced clinical trial
	# with 'sample_size' subjects per arm using a binomial distribution.
	# 'prob1' is the proportion for events in controls.  'odds_ratio is
    # the ratio of group odds (experimental group over control group).
	#
	# In the code below y1 and prob1 are data for the control goup.
	# In the code below y2 and prob1 are data for the experimental
	# group.
	# --------------------------------------------------------------- #

    # prob1 is the bernoulli distribution event probability parameter
	# for the control group. Given a user specified odds ratio, I
	# need to calculate the event probability for the experimental
	# group.
    A <- odds_ratio * (prob1 / (1 - prob1))
    prob2 <- A / (1 + A)

	# Create outcomes for both groups.
    y1 <- stats::rbinom(sample_size, size = 1, prob = prob1)
    y2 <- stats::rbinom(sample_size, size = 1, prob = prob2)

	#Take all data created above and put into a data frame that contains
	# the required variables.
    subjid <- seq(from = 1, to = 2 * sample_size)
    trt <- c(rep(0, sample_size), rep(1, sample_size))
    y <- c(y1, y2)

    gendata <- data.frame(subjid, trt, y)
    colnames(gendata) <- c("id", "treatment", "y")

    return(gendata)
}


#' Log-likelihood function for two-arm trial with historical data using Bernoulli
#' distribution.
#'
#' \code{bernoulliloglike()} function only used internally by
#' \code{bernoullitrialsimulator()} function to estimate Bernoulli model parameters
#' when clinical trial involves experimental and control groups as well as historical
#' control data.  The Bernoulli log-likelihood is calculated by modeling \code{data}
#' as a Bernoulli random variable. Not to be called directly by user.
#'
#' @param params  Two element vector of Bernoulli parameters.  The first and second elements
#'   are the intercept (beta0) and treatment effect parameter (beta1), where the treatment effect is
#'   a log odds ratio (experimental group over control group).  The prob parameter required by
#'   dbinom() is equal to exp(params[1] + params[2]*treatment) / (1 + exp(params[1] + params[2]*treatment)).
#'   It is assumed that the params[1] parameter is the same in both randomized and historical data.
#'   It is assumed that the prob parameter for dbinom() in the randomized and historical control data is
#'   equal to exp(params[1]) / (1 + exp(params[1])).
#' @param randdata  Dataset of randomly generated trial data.  Randomized trial datasets
#'   must have 3 columns: id, treatment, and y.  The value of treatment must be 0 (control)
#'   or 1 (experimental).  The values of y must be 0 or 1.
#' @param histdata Dataset of historical data.  Historical datasets must have 3 columns: id,
#'   treatment, and y.  The value of treatment should be 0.  The values of y must be 0 or 1.
#' @param a0 Power prior parameter: 0 implies historical data is ignored and 1 implies
#'   all information in historical data is used.
#'
#' @return \code{bernoulliloglike()} returns a value of the loglikelihood function
#'   given a set of Bernoulli parameters, randomly generated trial data, and observed
#'   historical data.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
bernoulliloglike <- function(params, randdata, histdata, a0) {

    # --------------------------------------------------------------- #
	#  This function calculates the Bernoulli log-likelihood given
	#  a vector of parameter values, a dataset of randomized trial
	#  data (two arms, no covariates beyond treatment), and a dataset
	#  of historical control data.
	#  This function has two parameters, beta0 and beta1.
	#  beta0 is the control group log odds, and beta 1 is the log
	#  odds ratio (experimental group over control group).
    # --------------------------------------------------------------- #

    # Get params
    beta0 <- params[1]
    beta1 <- params[2]

    # Calculate the logit vector for all randomized observations.
	# beta0 and beta1 are regression parameters linked to the event
	# probabilities via the logit function.  beta1 is the log odds
	# ratio (experimental group over conrol group).  beta0 is the log
	# odds among of the control group.
    logit_i <- beta0 + beta1 * randdata$treatment
	# Using the logit vector, calculate event probabilities for all
	# randomized observations.
    prob_i <- exp(logit_i) / (1 + exp(logit_i))

	# Calculate the log-likelihood values for all randomized observations.
    ll_R <- stats::dbinom(randdata$y, size = 1, prob = prob_i, log = TRUE)

	# Calculate the event probability vector for all historical control
	# observations.  Note that this event probability is the same
	# as the event probability for all randomized control observations.
    probH_i <- exp(beta0) / (1 + exp(beta0))
	# Calculate the loglikelihood values for all historical control
	# observations.
    ll_H <- stats::dbinom(histdata$y, size = 1, prob = probH_i, log = TRUE)

	# Calculate the overall log likelihood by adding the randomized log-likelihood to the historical control
	# log-likelihood by a0, where a0 is the power prior parameter.  This a0 value is defined by the
	# user and not estimated via object function optimization.
    ll <- sum(ll_R) + a0 * sum(ll_H)

	# Return the sum of all individual elements to the negative log-likelihood
    return(-ll)
}


#' Simulate a single randomized trial using a Bernoulli outcome and information from
#' historical controls (Logistic regression model).
#'
#' \code{bernoullitrialsimulator()} function only used internally by
#' \code{bernoulli_sim()} function to run a single trial simulation involving historical
#' control data and a Bernoulli (0/1) outcome.
#'
#' The simulation of a trial with a Bernoulli outcome involving historical control data returns
#' an estimate of the odds ratio as well as an estimate of the log odds ratio variance.
#' Finally the simulation returns an indication of whether or not the simulated trial led to
#' a rejection of the null hypothesis (1) or not (0).
#'
#' \code{bernoullitrialsimulator()} should not be called directly by user.
#'
#'
#' @param sample_size_val Number of subjects per arm.
#' @param histdata Dataset of historical data.  Historical datasets must have 3 columns: id,
#'   treatment, and y.  The value of treatment should be 0.  The values of y must be 0 or 1.
#' @param prob1_val prob parameter value for randomized control arm.  Used in call to \code{rbinom()}.
#' @param odds_ratio_val Desired odds ratio between randomized experimental and control groups.
#' @param a0_val A power prior parameter ranging from 0 to 1, where 0
#'   implies no information from historical data should be used, 1 implies all of
#'   the information from historical data should be used.  A value between 0 and 1
#'   implies that a proportion of the information from historical data will be used.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#'
#' @return \code{bernoullitrialsimulator()} returns a vector of simulation results. The
#'   first element is an estimated odds ratio, the second element is the estimated
#'   variance of the log odds ratio, and the third element is a 0/1 variable indicator
#'   whether or not the trial rejected the null hypothesis (1) or failed to reject
#'   the null hypothesis (0).
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
bernoullitrialsimulator <- function(sample_size_val, histdata, prob1_val, odds_ratio_val, a0_val, alpha) {

    # --------------------------------------------------------------- #
    #  This function simulates a two-arm Bayesian trial where
    #  historical data is utilized in the parameter estimation.
    # --------------------------------------------------------------- #

    # First, Generate Bernoulli trial data given the user defined trial characteristics.
    sampleranddata <- genbernoullidata(sample_size = sample_size_val, prob1 = prob1_val, odds_ratio = odds_ratio_val)

    # Generate initial values for your call to optim()
    initializemodel <- stats::glm(y ~ treatment, family = stats::binomial(link = "logit"), data = sampleranddata)

    initialbeta0 <- initializemodel$coefficients[1]
    initialbeta1 <- initializemodel$coefficients[2]

	# Generate the Bayesian CLT based parameter estimates needed for inference on odds ratio.
    fitmod <- stats::optim(c(initialbeta0, initialbeta1), bernoulliloglike, randdata = sampleranddata, histdata = histdata,
        a0 = a0_val, method = "Nelder-Mead", hessian = TRUE)

	#Extract model parameters and statistics
    modparm <- fitmod$par
    covarmat <- solve(fitmod$hessian)

    logoddsratio <- modparm[2]

    odds_ratio <- exp(logoddsratio)
    lower_oddsratio <- exp(logoddsratio - stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2]))
    upper_oddsratio <- exp(logoddsratio + stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2]))

	#Make a decision about the simulated trial, reject or fail to reject null hypothesis.
    reject <- ifelse(((lower_oddsratio > 1) | (upper_oddsratio < 1)), 1, 0)
    output <- c(odds_ratio, covarmat[2, 2], reject)

	#Return the odds ratio, the estimated variance of the log odds ratio, and the trial decision.
    names(output) <- c("odds_ratio", "log_or_var", "reject")
    return(output)

}


#' Simulate a single randomized trial using a Bernoulli outcome but not including any information from
#' historical controls (Logistic regression model).
#'
#' \code{bernoullitrialsimulator()} function only used internally by
#' \code{simple_bernoulli_sim()} function to run a single trial simulation involving historical
#' control data and a Bernoulli (0/1) outcome.
#'
#' The simulation of a trial with a Bernoulli outcome without historical control data returns
#' an estimate of the odds ratio as well as an estimate of the log odds ratio variance.
#' Finally the simulation returns an indication of whether or not the simulated trial led to
#' a rejection of the null hypothesis (1) or not (0).
#'
#' \code{bernoullitrialsimulatornohist()} should not be called directly by user.
#'
#' @param sample_size_val Number of subjects per arm.
#' @param prob1_val prob parameter value for randomized control arm.  Used in call to \code{rbinom()}.
#' @param odds_ratio_val Desired odds ratio between randomized experimental and control groups.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#'
#' @return \code{bernoullitrialsimulatornohist()} returns a vector of simulation results. The
#'   first element is an estimated odds ratio, the second element is the estimated
#'   variance of the log odds ratio, and the third element is a 0/1 variable indicator
#'   whether or not the trial rejected the null hypothesis (1) or failed to reject
#'   the null hypothesis (0).
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
bernoullitrialsimulatornohist <- function(sample_size_val, prob1_val, odds_ratio_val, alpha) {

    # --------------------------------------------------------------- #
    #  This function simulates a two-arm Bayesian trial where
    #  historical data is not utilized in the parameter estimation.
	#  No call to optim is necessary since using flat priors and
	#  canonical link (expected and observed information are identical)
    # --------------------------------------------------------------- #

    # First, Generate Bernoulli trial data given the user defined trial characteristics.
    sampleranddata <- genbernoullidata(sample_size = sample_size_val, prob1 = prob1_val, odds_ratio = odds_ratio_val)

    # Generate the Bayesian CLT based parameter estimates needed for inference on odds ratio.
    initializemodel <- stats::glm(y ~ treatment, family = stats::binomial(link = "logit"), data = sampleranddata)

    modparm <- initializemodel$coefficients
    covarmat <- stats::vcov(initializemodel)

    logoddsratio <- modparm[2]

    odds_ratio <- exp(logoddsratio)
    lower_oddsratio <- exp(logoddsratio - stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2]))
    upper_oddsratio <- exp(logoddsratio + stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2]))

	#Make a decision about the simulated trial, reject or fail to reject null hypothesis.
    reject <- ifelse(((lower_oddsratio > 1) | (upper_oddsratio < 1)), 1, 0)

	#Return the odds ratio, the estimated variance of the log odds ratio, and the trial decision.
    output <- c(odds_ratio, covarmat[2, 2], reject)
    names(output) <- c("odds_ratio", "log_or_var", "reject")
    return(output)

}


#' Repeated Two Arm Bayesian Clinical Trial Simulation with Historical Data and
#' Bernoulli Outcome (Logistic regression model).
#'
#' \code{bernoulli_sim()} function only used internally by \code{historic_sim()}
#' function to run a set of trial simulations involving historical
#' control data and a Bernoulli outcome.  User defined simulation parameters are
#' used to generate a set of trial scenarios.  Each scenario is simulated multiple
#' times and then means are taken to calculate estimates of power, odds ratios,
#' and other user requested study summary statistics like variance of odds
#' ratio, bias (on odds ratio scale), and mse (on odds ratio scale).
#' The number of repeated simulations is defined by the user.
#'
#' \code{bernoulli_sim()} should not be called directly by user.
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
#' @param effect_vals A vector of odds ratios (randomized experimental over control),
#'   all of which must be positive.
#' @param rand_control_diff For Bernoulli outcomes this is a vector of odds ratios
#'   (randomized controls over historical controls) that represent odds ratios
#'   between randomized and historical controls.
#' @param hist_control_data A dataset of historical data.  Default is \code{NULL}.
#'   Historical datasets must have 3 columns: id, treatment, and y.  The value of
#'   treatment should be 0.  The values of y must be 0 or 1.
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
#' @return \code{bernoulli_sim()} returns an S3 object of class bayes_ctd_array.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
bernoulli_sim <- function(trial_reps=100, subj_per_arm, a0_vals, effect_vals,
                          rand_control_diff, hist_control_data, alpha=0.05,
                          get_var=FALSE, get_bias=FALSE, get_mse=FALSE,
                          quietly=TRUE) {

    # --------------------------------------------------------------- #
    # For a set of user specified scenarios (defined by combinations
    # of user specified parameters), simulate "trial_reps" trials
    # and estimate power, odds ratio estimate, and if requested
    # by user: variance of odds ratio, bias, and mse.  Using a
    # Bernoulli oucome and incorporating data from historical controls.
    # --------------------------------------------------------------- #

  # Need to take the historical data and generate distributional parameter estimates
    histdata = hist_control_data
    hist_model <- stats::glm(y ~ 1, family = stats::binomial(link = "logit"), data = histdata)
    initialprob1 <- exp(hist_model$coefficients[1])/(1 + exp(hist_model$coefficients[1]))


    # Initialize arrays to hold power, mse, and bias estimate results as requested.
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
        # Need to adjust the randomized control odds given the historical control odds and the odds ratios given in
        # rand_control_diff, then calculate the event prob for randomized controls.
        rand_cont_odds <- (initialprob1 / (1 - initialprob1)) * rand_control_diff[diffs]
        adjprob1 <- rand_cont_odds / (1 + rand_cont_odds)

        for (effvals in 1:length(effect_vals)) {
            for (a0vals in 1:length(a0_vals)) {
                for (sizes in 1:length(subj_per_arm)) {
                  if (!quietly){
                    cat("\r", c(subj_per_arm[sizes], a0_vals[a0vals], effect_vals[effvals], rand_control_diff[diffs]))
                  }
                  # For each combination of rand_control_diff, effect_vals, a0_vals, and subj_per_arm, simulate the trial
                  # trial_reps times and then calculate the mean reject rate to estimate power.  For bias, work on the
                  #odds ratio scale and take mean of differences between estimated odds ratios and the true
                  #odds ratio.  For mse, calculate the mean of squared differences between the estimated odds ratios
                  #and the true odds ratio.  Note that rand_control_diff is set to 1 and a0_val is set to 0.
                  collect <- matrix(rep(0, 3 * trial_reps), ncol = 3)
                  for (k in 1:trial_reps) {
                    # sample_size_val will be equal to both arms
                    collect[k, ] <- bernoullitrialsimulator(sample_size_val = subj_per_arm[sizes], histdata, prob1_val = adjprob1,
                      odds_ratio_val = effect_vals[effvals], a0_val = a0_vals[a0vals], alpha = alpha)
                  }
                  #collect is a matrix of data, odds ratio in 1st column, log odds ratio variance
                  # in second column, and a vector of 0/1s in third column indicating whether or
                  # not trial represented by row led to a rejection of null hypothesis (1) or not (0).
                  # Note that collect gets rewritten for each scenario.
                  colnames(collect) <- c("odds_ratio", "log_or_var", "reject")
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

    #Lines 404 through 707 simply apply names to the dimensions of array created by the
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
#' Bernoulli Outcome (Logistic regression model).
#'
#' \code{simple_bernoulli_sim()} function only used internally by \code{simple_sim()}
#' function to run a set of trial simulations involving no historical
#' control data and a Bernoulli outcome.  User defined simulation parameters are
#' used to generate a set of trial scenarios.  Each scenario is simulated multiple
#' times and then means are taken to calculate estimates of power, odds ratios,
#' and other user requested study summary statistics like variance of odds
#' ratio, bias (on odds ratio scale), and mse (on odds ratio scale).
#' The number of repeated simulations is defined by the user.
#'
#' \code{simple_bernoulli_sim()} should not be called directly by user.
#'
#' @param trial_reps Number of trials to replicate within each combination of
#'   subj_per_arm and effect_vals.  As the number of trials increases, the
#'   precision of the estimate will increase. Default is 100.
#' @param subj_per_arm A vector of sample sizes, all of which must be positive
#'   integers.
#' @param effect_vals A vector of odds ratios (randomized experimental over control),
#'   all of which must be positive.
#' @param prob1_val prob parameter value for randomized control arm. Used in call to \code{rbinom()}.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#' @param get_var A TRUE/FALSE indicator of whether or not an array of variance
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param get_bias A TRUE/FALSE indicator of whether or not an array of bias
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param get_mse A TRUE/FALSE indicator of whether or not an array of MSE
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param quietly A TRUE/FALSE indicator of whether notes are printed
#'   to output about simulation progress as the simulation runs.  If
#'   running interactively in RStudio or running in the R console,
#'   \code{quietly} can be set to FALSE.  If running in a Notebook or
#'   knitr document, \code{quietly} needs to be set to TRUE.  Otherwise
#'   each note will be printed on a separate line and it will take up
#'   a lot of output space.  Default is \code{TRUE}.
#'
#' @return \code{simple_bernoulli_sim()} returns an S3 object of class bayes_ctd_array.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
simple_bernoulli_sim <- function(trial_reps=100, subj_per_arm, effect_vals, prob1_val,
                                 alpha=0.05, get_var=FALSE, get_bias=FALSE,
                                 get_mse=FALSE, quietly=TRUE) {

  # --------------------------------------------------------------- #
  # For a set of user specified scenarios (defined by combinations
  # of user specified parameters), simulate "trial_reps" trials
  # and estimate power, odds ratio estimate, and if requested by user:
  # variance of odds ratio, bias, and mse.  Using a Bernoulli oucome
  # but historical control data is not used.
  # --------------------------------------------------------------- #

  #The rand_control_diff and a0_val dimensions will be set to 1, and the value for
  # rand_control_diff will be 1 and a0_val will be set to 0.  All summaries will
  # be set up to ignore these dimensions for simple (no historical data) simulations.
    rand_control_diff <- 1
    a0_vals <- 0
    # Initialize arrays to hold power, mse, and bias estimate results as requested.
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
    # requested trial characteristics. Note that rand_control_diff is set to 1 and a0_val is set to 0.
    for (diffs in 1:length(rand_control_diff)) {
        for (effvals in 1:length(effect_vals)) {
            for (a0vals in 1:length(a0_vals)) {
                for (sizes in 1:length(subj_per_arm)) {
                  if (!quietly){
                    cat("\r", c(subj_per_arm[sizes], a0_vals[a0vals], effect_vals[effvals], rand_control_diff[diffs]))
                  }
                  # For each combination of rand_control_diff, effect_vals, a0_val, and subj_per_arm, simulate the trial
                  # trial_reps times and then calculate the mean reject rate to estimate power.  For bias, work on the
                  #odds ratio scale and take mean of differences between estimated odds ratios and the true
                  #odds ratio.  For mse, calculate the mean of squared differences between the estimated odds ratios
                  #and the true odds ratio.  Note that rand_control_diff is set to 1 and a0_val is set to 0.
                  collect <- matrix(rep(0, 3 * trial_reps), ncol = 3)
                  for (k in 1:trial_reps) {
                    # sample_size_val will be equal to both arms
                    collect[k, ] <- bernoullitrialsimulatornohist(sample_size_val = subj_per_arm[sizes], prob1_val = prob1_val,
                      odds_ratio_val = effect_vals[effvals], alpha = alpha)
                  }
                  #collect is a matrix of data, odds ratio in 1st column, log odds ratio variance
                  # in second column, and a vector of 0/1s in third column indicating whether or
                  # not trial represented by row led to a rejection of null hypothesis (1) or not (0).
                  # Note that collect gets rewritten for each scenario.
                  colnames(collect) <- c("odds_ratio", "log_or_var", "reject")
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

    #Lines 830 through 1133 simply apply names to the dimensions of array created by the
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
