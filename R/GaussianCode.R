#' Generating function for Gaussian Data.
#'
#' \code{gengaussiandata()} function used mainly internally by
#' \code{gaussiantrialsimulator()} function to generate data for a two-arm
#' clinical trial, experimental and control groups.  Can be used to generate
#' random trial data.
#'
#' @param sample_size  Number of subjects per arm.
#' @param mu1 mean parameter used in call to \code{rnorm()}.
#'   Used only in control arm.
#' @param mean_diff Desired Mean Difference between experimental and control groups.
#' @param common_sd sd parameter used in call to \code{rnorm()}.
#'   Used in both arms.
#'
#' @return \code{gengaussiandata()} returns a data frame with columns: 'id', 'treatment',
#'   and 'y'.
#'
#' @examples
#' samplehistdata <- gengaussiandata(sample_size=60, mu1=25, mean_diff=0, common_sd=3)
#' samplehistdata
#' @export
gengaussiandata <- function(sample_size, mu1, mean_diff, common_sd) {

    # --------------------------------------------------------------- #
	# The function gengaussiandata simulates a balanced clinical trial
	# with 'sample_size' subjects per arm using a normal distribution.
	# 'm'1 is the gaussian location parameter, and 'common_sd' is the
	# gaussian scale parameter for both arms.  'mean_diff' is the
	# difference in group means (experimental group minus control group).
	#
	# In the code below y1 and mu1 are data for the control goup.
	# In the code below y2 and mu2 are data for the experimental group.
    # --------------------------------------------------------------- #

    # mu1 is the normal distribution mu parameter for the control group.
	# given mean_diff, normal mean for experimental group is:
    mu2 <- mu1 + mean_diff

	# Create outcomes for both groups.
    y1 <- stats::rnorm(sample_size, mean = mu1, sd = common_sd)
    y2 <- stats::rnorm(sample_size, mean = mu2, sd = common_sd)

	#Take all data created above and put into a data frame that contains
	# the required variables.
    subjid <- seq(from = 1, to = 2 * sample_size)
    trt <- c(rep(0, sample_size), rep(1, sample_size))
    y <- c(y1, y2)

    gendata <- data.frame(subjid, trt, y)
    colnames(gendata) <- c("id", "treatment", "y")

    return(gendata)
}


#' Log-likelihood function for two-arm trial with historical data using Gaussian
#' distribution.
#'
#' \code{gaussianloglike()} function used only used internally by
#' \code{gaussiantrialsimulator()} function to estimate Gaussian model parameters
#' when clinical trial involves experimental and control groups as well as historical
#' control data.  The Gaussian log-likelihood is calculated by modeling \code{data}
#' as a Gaussian random variable. Not to be called directly by user.
#'
#' @param params  Three element vector of Gaussian parameters.  Third element is log(sd),
#'   where sd is a parameter required by dnorm().  The first and second elements
#'   are the intercept and treatment effect parameter, where the treatment effect is
#'   a mean difference (experimental group minus control group).  The mu parameter required by
#'   dnorm() is equal to params[1] + params[2]*treatment.  It is assumed that the log(sd)
#'   parameter is the same in both randomized and historical data.  It is assumed that
#'   the mu parameter in the randomized and historical control data is equal to params[1].
#' @param randdata  Dataset of randomly generated trial data.  Randomized trial datasets
#'   must have 3 columns: id, treatment, and y.  The value of treatment must be 0 (control)
#'   or 1 (experimental).  The values of y must be numeric.
#' @param histdata Dataset of historical data.  Historical datasets must have 3 columns: id,
#'   treatment, and y.  The value of treatment should be 0.  The values of y must be
#'   numeric.
#' @param a0 Power prior parameter: 0 implies historical data is ignored and 1 implies
#'   all information in historical data is used.
#'
#' @return \code{gaussianloglike()} returns a value of the loglikelihood function
#'   given a set of Gaussian parameters, randomly generated trial data, and observed
#'   historical data.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
gaussianloglike <- function(params, randdata, histdata, a0) {

    # --------------------------------------------------------------- #
	#  This function calculates the Gaussian log-likelihood given
	#  a vector of parameter values, a dataset of randomized trial
	#  data (two arms, no covariates beyond treatment), and a dataset
	#  of historical control data.
	#  The standard deviation parameter is common in both randomized
	#  groups and the historical control group.  The mean parameter is
	#  assumed to be the same in both control groups.  The mean
    #  parameter for the randomized experimental group is the sum
    #  of the control mean and the treatment effect.  The parameters
    #  are beta0, beta1, and logs.  logs is the common log standard
	#  deviation in all three groups.  It follows that the common
    #  standard deviation is equal to exp(logs).  beta0 and beta1
	#  are regression parameters that are linked to the gaussian mean
	#  via the identity function.
	#  beta0 is the control group mean, and beta 1 is the experimental
	#  group mean minus the control group mean.
    # --------------------------------------------------------------- #

    # Get params
    beta0 <- params[1]
    beta1 <- params[2]
    logs <- params[3]

    # Calculate the mean parameter vector for all randomized observations.
	# Note that beta0 is simply the control mean and beta1 is the
	# experimental group mean minus the control group mean.
    mu_i <- beta0 + beta1 * randdata$treatment
	# Calculate the log-likelihood values for all randomized observations.
    ll_R <- stats::dnorm(randdata$y, mean = mu_i, sd = exp(logs), log = TRUE)

	# Calculate the log-likelihood values for all historical control observations.
	# Note that it is assumed the mean among randomized and historical controls
	# is the same.
    ll_H <- stats::dnorm(histdata$y, mean = beta0, sd = exp(logs), log = TRUE)

	# Calculate the overall log-likelihood by adding the randomized log-likelihood to the historical control
	# log-likelihood by a0, where a0 is the power prior parameter.  This a0 value is defined by the
	# user and not estimated via object function optimization.
    ll <- sum(ll_R) + a0 * sum(ll_H)
    return(-ll)
}


#' Simulate a single randomized trial using a Gaussian outcome and information from
#' historical controls.
#'
#' \code{gaussiantrialsimulator()} function only used internally by
#' \code{gaussian_sim()} function to run a single trial simulation involving historical
#' control data and a Gaussian outcome.
#'
#' The simulation of a trial with a Gaussian outcome involving historical control data returns
#' an estimate of the mean difference as well as an estimate of the mean difference variance.
#' Finally the simulation returns an indication of whether or not the simulated trial led to
#' a rejection of the null hypothesis (1) or not (0).
#'
#' \code{gaussiantrialsimulator()} should not be called directly by user.
#'
#' @param sample_size_val Number of subjects per arm.
#' @param histdata Dataset of historical data.  Historical datasets must have 3 columns: id,
#'   treatment, and y.  The value of treatment should be 0.  The values of y must be
#'   numeric.
#' @param mu1_val Randomized control arm mean parameter used in call to \code{rnorm()}.
#' @param mean_diff_val Desired mean difference between randomized experimental and control groups.
#' @param common_sd_val Randomized sd parameter used in call to \code{rnorm()}.
#'   Used in both randomized arms.
#' @param a0_val A power prior parameter ranging from 0 to 1, where 0
#'   implies no information from historical data should be used, 1 implies all of
#'   the information from historical data should be used.  A value between 0 and 1
#'   implies that a proportion of the information from historical data will be used.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#'
#' @return \code{gaussiantrialsimulator()} returns a vector of simulation results. The
#'   first element is an estimated mean difference, the second element is the estimated
#'   variance of the mean difference, and the third element is a 0/1 variable indicator
#'   whether or not the trial rejected the null hypothesis (1) or failed to reject
#'   the null hypothesis (0).
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
gaussiantrialsimulator <- function(sample_size_val, histdata, mu1_val, mean_diff_val, common_sd_val, a0_val, alpha) {

    # --------------------------------------------------------------- #
    #  This function simulates a two-arm Bayesian trial where
    #  historical data is utilized in the parameter estimation.
    # --------------------------------------------------------------- #

	# First, Generate Gaussian trial data given the user defined trial characteristics.
    sampleranddata <- gengaussiandata(sample_size = sample_size_val, mu1 = mu1_val, mean_diff = mean_diff_val, common_sd = common_sd_val)

    # Generate initial values for your call to optim()
    initializemodel <- stats::lm(y ~ treatment, data = sampleranddata)

    initialbeta0 <- initializemodel$coefficients[1]
    initialbeta1 <- initializemodel$coefficients[2]
    initiallogs <- log(summary(initializemodel)$sigma)

	# Generate the Bayesian CLT based parameter estimates needed for inference on hazard ratio.
    fitmod <- stats::optim(c(initialbeta0, initialbeta1, initiallogs), gaussianloglike, randdata = sampleranddata, histdata = histdata,
        a0 = a0_val, method = "Nelder-Mead", hessian = TRUE)

	#Extract model parameters and statistics
    modparm <- fitmod$par
    covarmat <- solve(fitmod$hessian)

    gaussianmeandiff <- modparm[2]

    lower_gaussianmeandiff <- gaussianmeandiff - stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2])
    upper_gaussianmeandiff <- gaussianmeandiff + stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2])

	#Make a decision about the simulated trial, reject or fail to reject null hypothesis.
    reject <- ifelse(((lower_gaussianmeandiff > 0) | (upper_gaussianmeandiff < 0)), 1, 0)
    output <- c(gaussianmeandiff, covarmat[2, 2], reject)

	#Return the mean, the estimated variance of the mean difference, and the trial decision.
    names(output) <- c("mean_diff", "mean_diff_var", "reject")
    return(output)

}


#' Simulate a single randomized trial using a Gaussian outcome but not including any information from
#' historical controls.
#'
#' \code{gaussiantrialsimulatornohist()} function only used internally by
#' \code{simple_gaussian_sim()} function to estimate Gaussian model parameters
#' when clinical trial involves experimental and control groups but no historical control
#' data.
#'
#' The simulation of a trial with a Gaussian outcome involving no historical control data returns
#' an estimate of the mean difference as well as an estimate of the mean difference variance.
#' Finally the simulation returns an indication of whether or not the simulated trial led to
#' a rejection of the null hypothesis (1) or not (0).
#'
#' \code{gaussiantrialsimulatornohist()} should not be called directly by user.
#'
#' @param sample_size_val Number of subjects per arm.
#' @param mu1_val Randomized control arm mean parameter used in call to \code{rnorm()}.
#' @param mean_diff_val Desired mean difference between randomized experimental and control groups.
#' @param common_sd_val Randomized sd parameter used in call to \code{rnorm()}.
#'   Used in both randomized arms.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#'
#' @return \code{gaussiantrialsimulatornohist()} returns a vector of simulation results. The
#'   first element is an estimated mean difference, the second element is the estimated
#'   variance of the mean difference, and the third element is a 0/1 variable indicator
#'   whether or not the trial rejected the null hypothesis (1) or failed to reject
#'   the null hypothesis (0).
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
gaussiantrialsimulatornohist <- function(sample_size_val, mu1_val, mean_diff_val, common_sd_val, alpha) {

    # --------------------------------------------------------------- #
    #  This function simulates a two-arm Bayesian trial where
    #  historical data is not utilized in the parameter estimation.
	  #  No call to optim is necessary since using flat priors and
	  #  canonical link (expected and observed information are identical)
    # --------------------------------------------------------------- #

    # First, Generate Gaussian trial data given the user defined trial characteristics.
    sampleranddata <- gengaussiandata(sample_size = sample_size_val, mu1 = mu1_val, mean_diff = mean_diff_val, common_sd = common_sd_val)

    # Generate the Bayesian CLT based parameter estimates needed for inference on mean difference.
    X = matrix(c(rep(1,length(sampleranddata$y)),sampleranddata$treatment),ncol=2)
    Beta_hat    = solve(t(X)%*%X)%*%t(X)%*%sampleranddata$y
    mleVal = X%*%Beta_hat
    Sigma_hat =  sqrt(sum((sampleranddata$y-mleVal)^2))/sqrt(2*sample_size_val)
    Covarmatinv <- matrix(rep(0,9),ncol=3)
    Covarmatinv[1:2,1:2] = t(X)%*%X/(Sigma_hat^2)
    Covarmatinv[1:2,3] = c(0,0)
    Covarmatinv[3,1:2] = c(0,0)
    Covarmatinv[3,3] = (Sigma_hat^(-6))*sum((sampleranddata$y-mleVal)^2)-(sample_size_val/(Sigma_hat^4))
    Covarmatinv = solve(Covarmatinv)


    #Extract model parameters and statistics
    mle <- Beta_hat
    covarmat <- Covarmatinv

    gaussianmeandiff <- mle[2]

    lower_gaussianmeandiff <- gaussianmeandiff - stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2])
    upper_gaussianmeandiff <- gaussianmeandiff + stats::qnorm(1 - alpha/2) * sqrt(covarmat[2, 2])

	#Make a decision about the simulated trial, reject or fail to reject null hypothesis.
    reject <- ifelse(((lower_gaussianmeandiff > 0) | (upper_gaussianmeandiff < 0)), 1, 0)
    output <- c(gaussianmeandiff, covarmat[2, 2], reject)

	#Return the mean difference, the estimated variance of the mean difference, and the trial decision.
    names(output) <- c("mean_diff", "mean_diff_var", "reject")
    return(output)

}


#' Repeated Two Arm Bayesian Clinical Trial Simulation with Historical Data and
#' Gaussian Outcome.
#'
#' \code{gaussian_sim()} function only used internally by \code{historic_sim()}
#' function to run a set of trial simulations involving historical
#' control data and a Gaussian outcome.  User defined simulation parameters are
#' used to generate a set of trial scenarios.  Each scenario is simulated multiple
#' times and then means are taken to calculate estimates of power, mean difference,
#' and other user requested study summary statistics like variance of mean
#' difference, bias, and mse.  The number of repeated simulations is defined
#' by the user.
#'
#' \code{gaussian_sim()} should not be called directly by user.
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
#' @param effect_vals A vector of mean differences (randomized experimental minus control).
#' @param rand_control_diff For Gaussian outcomes this is a vector of mean differences
#'   (randomized controls minus historical controls) that represent differences
#'   between randomized and historical controls.
#' @param hist_control_data A dataset of historical data.  Default is \code{NULL}.
#'   Historical datasets must have 3 columns: id, treatment, and y.  The value of
#'   treatment should be 0.  The values of y must be numeric.
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
#' @return \code{gaussian_sim()} returns an S3 object of class bayes_ctd_array.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
gaussian_sim <- function(trial_reps=100, subj_per_arm, a0_vals, effect_vals,
                         rand_control_diff, hist_control_data, alpha=0.05,
                         get_var=FALSE, get_bias=FALSE, get_mse=FALSE,
                         quietly=TRUE) {

    # --------------------------------------------------------------- #
    # For a set of user specified scenarios (defined by combinations
    # of user specified parameters), simulate "trial_reps" trials
    # and estimate power, mean difference estimate, and if requested
    # by user: variance of mean difference, bias, and mse.  Using a
    # Gaussian oucome and incorporating data from historical controls.
    # --------------------------------------------------------------- #

    # Need to take the historical data and generate distributional parameter estimates
    histdata <- hist_control_data
    initializemodel <- stats::lm(y ~ 1, data = histdata)
    initialmu1 <- initializemodel$coefficients[1]
    initialsig <- summary(initializemodel)$sigma


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
        # Need to adjust the randomized control mean given the historical control mean and the mean differences given in
        # rand_control_diff
        adjmu1 <- initialmu1 + rand_control_diff[diffs]

        for (effvals in 1:length(effect_vals)) {
            for (a0vals in 1:length(a0_vals)) {
                for (sizes in 1:length(subj_per_arm)) {
                  if (!quietly){
                    cat("\r", c(subj_per_arm[sizes], a0_vals[a0vals], effect_vals[effvals], rand_control_diff[diffs]))
                  }
                  # For each combination of rand_control_diff, effect_vals, a0_vals, and subj_per_arm, simulate the trial
                  #trial_reps times and then calculate the mean reject rate to estimate power.  For bias, work on the
                  #untransformed scale and take mean of differences between estimated mean differences and the true mean
                  #mean difference.  For mse, calculate the mean of squared differences between the estimated mean
                  #difference and the true mean difference.
                  collect <- matrix(rep(0, 3 * trial_reps), ncol = 3)
                  for (k in 1:trial_reps) {
                    # sample_size_val will be equal to both arms
                    collect[k, ] <- gaussiantrialsimulator(sample_size_val = subj_per_arm[sizes], histdata, mu1_val = adjmu1,
                      mean_diff_val = effect_vals[effvals], common_sd_val = initialsig, a0_val = a0_vals[a0vals], alpha = alpha)
                  }
                  #collect is a matrix of data, mean difference in 1st column, mean difference variance
                  # in second column, and a vector of 0/1s in third column indicating whether or
                  # not trial represented by row led to a rejection of null hypothesis (1) or not (0).
                  # Note that collect gets rewritten for each scenario.
                  colnames(collect) <- c("mean_diff", "mean_diff_var", "reject")
                  #Start calculating means for each scenarios and placing the means in the proper
                  # array.  Every simulation will contain an array of power results and mean
                  # difference estimates.
                  power_results[sizes, a0vals, effvals, diffs] <- mean(collect[, 3])
                  est_results[sizes, a0vals, effvals, diffs] <- mean(collect[, 1])
                  if (get_bias == TRUE) {
                    bias_results[sizes, a0vals, effvals, diffs] <- mean(collect[, 1] - effect_vals[effvals])
                  }
                  if (get_var == TRUE) {
                    var_results[sizes, a0vals, effvals, diffs] <- mean(collect[, 2])
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

    #Lines 416 through 719 simply apply names to the dimensions of array created by the
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
#' Gaussian Outcome.
#'
#' \code{simple_gaussian_sim()} function only used internally by \code{simple_sim()}
#' function to run a set of trial simulations involving no historical
#' control data and a Gaussian outcome.  User defined simulation parameters are
#' used to generate a set of trial scenarios.  Each scenario is simulated multiple
#' times and then means are taken to calculate estimates of power, mean difference,
#' and other user requested study summary statistics like variance of mean
#' difference, bias, and mse.  The number of repeated simulations is defined
#' by the user.
#'
#' \code{simple_gaussian_sim()} should not be called directly by user.
#'
#' @param trial_reps Number of trials to replicate within each combination of
#'   subj_per_arm and effect_vals.  As the number of trials increases, the
#'   precision of the estimate will increase. Default is 100.
#' @param subj_per_arm A vector of sample sizes, all of which must be positive
#'   integers.
#' @param effect_vals A vector of mean differences (randomized experimental minus
#'   controls).
#' @param mu1_val mean parameter value for randomized control arm. Used in call to
#'   \code{rnorm()}.
#' @param common_sd_val sd parameter value used in both randomized arms.  Used in call
#'   to \code{rnorm()}.
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
#' @return \code{simple_gaussian_sim()} returns an S3 object of class bayes_ctd_array.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
simple_gaussian_sim <- function(trial_reps=100, subj_per_arm, effect_vals, mu1_val,
                                common_sd_val, alpha=0.05, get_var=FALSE,
                                get_bias=FALSE, get_mse=FALSE, quietly=TRUE) {

    # --------------------------------------------------------------- #
    # For a set of user specified scenarios (defined by combinations
    # of user specified parameters), simulate "trial_reps" trials
    # and estimate power, mean difference estimate, and if requested
    # by user: variance of mean difference, bias, and mse.  Using a
    # Gaussian oucome but historical control data is not used.
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
                  #untransformed scale and take mean of differences between estimated mean differences and the true mean
                  #mean difference.  For mse, calculate the mean of squared differences between the estimated mean
                  #difference and the true mean difference.  Note that rand_control_diff is set to 1 and
                  #a0_val is set to 0.
                  collect <- matrix(rep(0, 3 * trial_reps), ncol = 3)
                  for (k in 1:trial_reps) {
                    # sample_size_val will be equal to both arms
                    collect[k, ] <- gaussiantrialsimulatornohist(sample_size_val = subj_per_arm[sizes], mu1_val = mu1_val,
                      mean_diff_val = effect_vals[effvals], common_sd_val = common_sd_val, alpha = alpha)
                  }
                  #collect is a matrix of data, mean difference in 1st column, mean difference variance
                  # in second column, and a vector of 0/1s in third column indicating whether or
                  # not trial represented by row led to a rejection of null hypothesis (1) or not (0).
                  # Note that collect gets rewritten for each scenario.
                  colnames(collect) <- c("mean_diff", "mean_diff_var", "reject")
                  #Start calculating means for each scenarios and placing the means in the proper
                  # array.  Every simulation will contain an array of power results and mean
                  # difference estimates.
                  power_results[sizes, a0vals, effvals, diffs] <- mean(collect[, 3])
                  est_results[sizes, a0vals, effvals, diffs] <- mean(collect[, 1])
                  if (get_bias == TRUE) {
                    bias_results[sizes, a0vals, effvals, diffs] <- mean(collect[, 1] - effect_vals[effvals])
                  }
                  if (get_var == TRUE) {
                    var_results[sizes, a0vals, effvals, diffs] <- mean(collect[, 2])
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

    #Lines 846 through 1149 simply apply names to the dimensions of array created by the
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
