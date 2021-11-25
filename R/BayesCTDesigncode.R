
#' Two Arm Bayesian Clinical Trial Simulation with Historical Data
#'
#' \code{historic_sim()} returns an S3 object of class \code{bayes_ctd_array}, which
#' will contain simulation results for power, statistic estimation, bias,
#' variance, and mse as requested by user.
#'
#' The object \code{bayes_ctd_array} has 6 elements: a list containing simulation
#' results (\code{data}), copies of the 4 function arguments \code{subj_per_arm},
#' \code{a0_vals}, \code{effect_vals}, and \code{rand_control_diff}, and finally
#' a \code{objtype} value indicating that \code{historic_sim()} was used. Each element of
#' \code{data} is a four-dimensional array, where each dimension is determined by the
#' length of parameters \code{subj_per_arm}, \code{a0_vals}, \code{effect_vals}, and
#' \code{rand_control_diff}. The size of \code{data} depends on which results are
#' requested by the user. At a minimum, at least one of \code{subj_per_arm},
#' \code{a0_vals}, \code{effect_vals}, or \code{rand_control_diff} must contain at
#' least 2 values, while the other three must contain at least 1 value. The \code{data}
#' list will always contain two elements: an array of power results (\code{power}) and
#' an array of estimation results (\code{est}). In addition to \code{power} and
#' \code{est}, data may also contain elements \code{var}, \code{bias}, or \code{mse},
#' depending on the values of \code{get_var}, \code{get_bias}, and \code{get_mse}. The
#' values returned in \code{est} are in the form of hazard ratios, mean ratios, odds
#' ratios, or mean differences depending on the value of \code{outcome_type}.  For a
#' Gaussian outcome, the estimation results are differences in group means (experimental
#' group minus control group). For a logistic outcome, the estimation results are odds
#' ratios (experimental group over control group). For lognormal and Poisson outcomes,
#' the estimation results are mean ratios (experimental group over control group). For a
#' piecewise exponential or a Weibull outcome, the estimation results are hazard
#' ratios (experimental group over control group).  The values returned in \code{bias},
#' \code{var}, and \code{mse} are on the scale of the values returned in
#' \code{est}.
#'
#' The object \code{bayes_ctd_array} has two primary methods, \code{print()} and
#' \code{plot()}, for printing and plotting slices of the arrays contained in
#' \code{bayes_ctd_array$data}.
#'
#' As dimensions of the four dimensional array increases, the time required to complete
#' the simulation will increase; however, it will be faster than a similar simulation
#' based on repeated calls to MCMC routines to analyze each simulated trial.
#'
#' The meaning of the estimation results, and the test used to generate power results,
#' depends on the outcome used. In all cases, power is based on a two-sided test
#' involving a (1-alpha)100\% credible interval, where the interval is used to determine
#' if the null hypothesis should be rejected (null value outside of the interval) or
#' not rejected (null value inside the interval). For a Gaussian outcome, the 95\%
#' credible interval is an interval for the difference in group means
#' (experimental group minus control group), and the test determines if 0 is in or
#' outside of the interval. For a Bernoulli outcome, the 95\% credible interval
#' is an interval for the odds ratio (experimental group over control group),
#' and the test determines if 1 is in or outside of the interval. For a lognormal or
#' a Poisson outcome, the 95\% credible interval is an interval for the mean ratio
#' (experimental group over control group), and the test determines if 1 is in or
#' outside of the interval. Finally, for a piecewise exponential or a Weibull outcome,
#' the 95\% credible interval is an interval for the hazard ratio (experimental group
#' over control group), and the test determines if 1 is in or outside of the interval.
#'
#' Please refer to the examples for illustration of package use.
#'
#' @param trial_reps Number of trials to replicate within each combination of
#'   \code{a0_vals}, \code{subj_per_arm}, \code{effect_vals}, and \code{rand_control_parms}.
#'   As the number of trials increases, the precision of the estimate will increase.
#'   Default is 100.
#' @param outcome_type Outcome distribution. Must be equal to \code{weibull},
#'   \code{lognormal}, \code{pwe} (Piecewise Exponential), \code{gaussian},
#'   \code{bernoulli}, or \code{poisson}.  Default is \code{weibull}.
#' @param subj_per_arm A vector of sample sizes, all of which must be positive
#'   integers.  Default is \code{c(50, 100, 150, 200, 250)}.
#' @param a0_vals A vector of power prior parameters ranging from 0 to 1, where 0
#'   implies no information from historical data should be used, and 1 implies all of
#'   the information from historical data should be used.  A value between 0 and 1
#'   implies that a proportion of the information from historical data will be used.
#'   Default is \code{c(0, 0.33, 0.67, 1)}.
#' @param effect_vals A vector of effects that should be reasonable for the
#'   outcome_type being studied, hazard ratios for Weibull, odds ratios for
#'   Bernoulli, mean ratios for Poisson, etc..  When \code{effect_vals} contain
#'   the null effect for a given \code{outcome_type}, the \code{power} component
#'   of \code{data} will contain an estimate of Type One Error.  In order to
#'   have a good set of Type One Error estimates, \code{trial_reps} need to be
#'   at least 10,000.  In such a case, if the total number of combinations
#'   made up from \code{subj_per_arm}, \code{a0_vals}, \code{effect_vals}, and
#'   \code{rand_control_diff} is very large, the time to complete the simulation
#'   can be substantial.  Default is \code{c(0.6, 1, 1.4)}.
#' @param rand_control_diff For piecewise exponential and Weibull outcomes, this is
#'   a vector of hazard ratios (randomized controls over historical controls)
#'   representing differences between historical and randomized controls.  For
#'   lognormal and Poisson outcomes, this is a vector of mean ratios (randomized
#'   controls over historical controls).  For a Bernoulli outcome, this is a vector
#'   of odds ratios (randomized controls over historical controls).  For a Gaussian
#'   outcome, this is a vector of mean differences (randomized minus historical
#'   controls). Default is \code{c(0.8, 1, 1.2)}.
#' @param hist_control_data A dataset of historical data.  Default is \code{NULL}.
#'   For survival outcomes, historical datasets must have 4 columns: id, treatment,
#'   event_time, and status.  The value of treatment should be 0.  For other
#'   outcomes, historical datasets must have columns: id, treatment, and y.
#' @param time_vec  A vector of time values which are used to create time periods
#'   within which the exponential hazard is constant.  Only used for piecewise
#'   exponential models.  Default is \code{NULL}.
#' @param censor_value A single value at which right censoring occurs when
#'   simulating randomized subject outcomes.  Used with survival outcomes.
#'   Default is \code{NULL}, where \code{NULL} implies no right censoring.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#' @param get_var A TRUE/FALSE indicator of whether an array of variance
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param get_bias A TRUE/FALSE indicator of whether an array of bias
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param get_mse A TRUE/FALSE indicator of whether an array of MSE
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param seedval A seed value for pseudo-random number generation.
#' @param quietly A TRUE/FALSE indicator of whether notes are printed
#'   to output about simulation progress as the simulation runs.  If
#'   running interactively in RStudio or running in the R console,
#'   \code{quietly} can be set to FALSE.  If running in a Notebook or
#'   knitr document, \code{quietly} needs to be set to TRUE.  Otherwise
#'   each note will be printed on a separate line and it will take up
#'   a lot of output space.  Default is \code{TRUE}.
#'
#' @return \code{historic_sim()} returns an S3 object of class \code{bayes_ctd_array}.
#'   As noted in details, an object of class \code{bayes_ctd_array }has 6 elements: a
#'   list of simulation results (\code{data}), copies of the 4 function arguments
#'   \code{subj_per_arm}, \code{a0_vals}, \code{effect_vals}, and
#'   \code{rand_control_diff}, and finally \code{objtype} indicating that \code{historic_sim()}
#'   was used. See details for a discussion about the contents of
#'   \code{data}. Results from the simulation contained in the \code{bayes_ctd_array}
#'   object can be printed or plotted using the \code{print()} and
#'   \code{plot()} methods. The results can also be accessed using basic list
#'   element identification and array slicing. For example, to get the 4-dimensional
#'   array of power results from a simulation, one could use the code
#'   \code{bayes_ctd_array$data$power}, where \code{bayes_ctd_array} is replaced
#'   with the name of the variable containing the \code{bayes_ctd_array} object. If
#'   one wanted a table of power for sample size by a0, while holding effect equal to
#'   the first considered value and control differences equal to the second considered
#'   value, then the code is \code{bayes_ctd_array$data$power[,,1,2]}, where
#'   \code{bayes_ctd_array} is replaced with the name of the variable containing the
#'   \code{bayes_ctd_array} object.
#'
#' @examples
#' #Generate a sample of historical data for use in example.
#' set.seed(2250)
#' SampleHistData <- genweibulldata(sample_size=60, scale1=2.82487,
#'                                  hazard_ratio=0.6, common_shape=3,
#'                                  censor_value=3)
#' histdata <- subset(SampleHistData, subset=(treatment==0))
#' histdata$id <- histdata$id+10000
#'
#' #Run a Weibull simulation, using historic_sim().
#' #For meaningful results, trial_reps needs to be much larger than 2.
#' weibull_test <- historic_sim(trial_reps = 2, outcome_type = "weibull",
#'                              subj_per_arm = c(50, 100, 150),
#'                              a0_vals = c(0, 0.50, 1),
#'                              effect_vals = c(0.6, 1),
#'                              rand_control_diff = c(0.8, 1),
#'                              hist_control_data = histdata, time_vec = NULL,
#'                              censor_value = 3, alpha = 0.05, get_var = TRUE,
#'                              get_bias = TRUE, get_mse = TRUE, seedval=123,
#'                              quietly=TRUE)
#'
#' #Tabulate the simulation results for power.
#' test_table <- print(x=weibull_test, measure="power",
#'                     tab_type="WX|YZ", effect_val=0.6,
#'                     rand_control_diff_val=1.0)
#' print(test_table)
#'
#' \donttest{
#' #Create a plot of the power simulation results.
#' plot(x=weibull_test, measure="power", tab_type="WX|YZ",
#'      smooth=FALSE, plot_out=TRUE, effect_val=0.6,
#'      rand_control_diff_val=1.0)
#' #Create a plot of the estimated hazard ratio simulation results.
#' plot(x=weibull_test, measure="est", tab_type="WX|YZ",
#'      smooth=FALSE, plot_out=TRUE, effect_val=0.6,
#'      rand_control_diff_val=1.0)
#' #Create a plot of the hazard ratio variance simulation results.
#' plot(x=weibull_test, measure="var", tab_type="WX|YZ",
#'      smooth=FALSE, plot_out=TRUE, effect_val=0.6,
#'      rand_control_diff_val=1.0)
#' #Create a plot of the hazard ratio bias simulation results.
#' plot(x=weibull_test, measure="bias", tab_type="WX|YZ",
#'      smooth=FALSE, plot_out=TRUE, effect_val=0.6,
#'      rand_control_diff_val=1.0)
#' #Create a plot of the hazard ratio mse simulation results.
#' plot(x=weibull_test, measure="mse", tab_type="WX|YZ",
#'      smooth=FALSE, plot_out=TRUE, effect_val=0.6,
#'      rand_control_diff_val=1.0)
#'
#' #Create other power plots using different values for tab_type
#' plot(x=weibull_test, measure="power", tab_type="XY|WZ",
#'      smooth=FALSE, plot_out=TRUE, subj_per_arm_val=150,
#'      rand_control_diff_val=1.0)
#'
#' plot(x=weibull_test, measure="power", tab_type="XZ|WY",
#'      smooth=FALSE, plot_out=TRUE, subj_per_arm_val=150, effect_val=0.6)
#'
#' plot(x=weibull_test, measure="power", tab_type="YZ|WX",
#'      smooth=FALSE, plot_out=TRUE, subj_per_arm_val=150, a0_val=0.5)
#'
#' plot(x=weibull_test, measure="power", tab_type="WY|XZ",
#'      smooth=FALSE, plot_out=TRUE, rand_control_diff_val=1, a0_val=0.5)
#'
#' plot(x=weibull_test, measure="power", tab_type="WZ|XY",
#'      smooth=FALSE, plot_out=TRUE, effect_val=0.6, a0_val=0.5)
#' }
#'
#' \donttest{
#' #Run Poisson simulation, using historic_sim(), but set two design characteristic
#' # parameters to only 1 value.
#' #Note: historic_sim() can take a while to run.
#' #Generate a sample of historical poisson data for use in example.
#' set.seed(2250)
#' samplehistdata <- genpoissondata(sample_size=60, mu1=1, mean_ratio=1.0)
#' histdata <- subset(samplehistdata, subset=(treatment==0))
#' histdata$id <- histdata$id+10000
#'
#' #For meaningful results, trial_reps needs to be larger than 100.
#' poisson_test <- historic_sim(trial_reps = 100, outcome_type = "poisson",
#'                               subj_per_arm = c(50, 75, 100, 125, 150, 175, 200, 225, 250),
#'                               a0_vals = c(1),
#'                               effect_vals = c(0.6),
#'                               rand_control_diff = c(0.6, 1, 1.6),
#'                               hist_control_data = histdata, time_vec = NULL,
#'                               censor_value = 3, alpha = 0.05, get_var = TRUE,
#'                               get_bias = TRUE, get_mse = TRUE, seedval=123,
#'                               quietly=TRUE)
#'
#' #Tabulate the simulation results for power.
#' test_table <- print(x=poisson_test, measure="power",
#'                     tab_type=NULL)
#' print(test_table)
#'
#' #Create a plot of the power simulation results.
#' plot(x=poisson_test, measure="power", tab_type=NULL,
#'      smooth=FALSE, plot_out=TRUE)
#' }
#'
#' \donttest{
#' #At least one of subj_per_arm, a0_vals, effect_vals, or rand_control_diff
#' #must contain at least 2 values.
#' #Generate a sample of historical lognormal data for use in example.
#' set.seed(2250)
#' samplehistdata <- genlognormaldata(sample_size=60, mu1=1.06, mean_ratio=0.6, common_sd=1.25,
#'                                    censor_value=3)
#' histdata <- subset(samplehistdata, subset=(treatment==0))
#' histdata$id <- histdata$id+10000
#'
#' #Run a Lognormal simulation, using historic_sim().
#' #For meaningful results, trial_reps needs to be larger than 100.
#' lognormal_test <- historic_sim(trial_reps = 100, outcome_type = "lognormal",
#'                                subj_per_arm = c(25,50,75,100,125,150,175,200,225,250),
#'                                a0_vals = c(1.0),
#'                                effect_vals = c(0.6),
#'                                rand_control_diff = c(1.8),
#'                                hist_control_data = histdata, time_vec = NULL,
#'                                censor_value = 3, alpha = 0.05, get_var = TRUE,
#'                                get_bias = TRUE, get_mse = TRUE, seedval=123,
#'                                quietly=TRUE)
#'
#' test_table <- print(x=lognormal_test, measure="power",
#'                     tab_type=NULL)
#' print(test_table)
#' #Create a plot of the power simulation results.
#' plot(x=lognormal_test, measure="power", tab_type=NULL,
#'      smooth=TRUE, plot_out=TRUE)
#' }
#'
#' @export
historic_sim <- function(trial_reps = 100, outcome_type = "weibull", subj_per_arm = c(50, 100, 150, 200, 250), a0_vals = c(0,
    0.33, 0.67, 1), effect_vals = c(0.6, 1, 1.4), rand_control_diff = c(0.8, 1, 1.2), hist_control_data = NULL, time_vec = NULL,
    censor_value = NULL, alpha = 0.05, get_var = FALSE, get_bias = FALSE, get_mse = FALSE, seedval=NULL, quietly=TRUE){

	#set random seed
	set.seed(seedval)
    #------------- Go through all the high level checks for proper input. -------------#
    global_error_checks(outcome_type, a0_vals, subj_per_arm, hist_control_data, rand_control_diff, get_var, get_bias,
        get_mse)

    #-- Given outcome type go through low level checks for proper input and run simulation --#
    if (tolower(outcome_type) == "weibull") {
        weibull_error_checks(effect_vals, hist_control_data, rand_control_diff, censor_value, alpha)
        results <- weibull_sim(trial_reps = trial_reps, subj_per_arm = subj_per_arm, a0_vals = a0_vals, effect_vals = effect_vals,
            rand_control_diff = rand_control_diff, hist_control_data = hist_control_data, censor_value = censor_value,
            alpha = alpha, get_var = get_var, get_bias = get_bias, get_mse = get_mse, quietly=quietly)

    } else if (tolower(outcome_type) == "lognormal") {
        lognormal_error_checks(effect_vals, hist_control_data, rand_control_diff, censor_value, alpha)
        results <- lognormal_sim(trial_reps = trial_reps, subj_per_arm = subj_per_arm, a0_vals = a0_vals, effect_vals = effect_vals,
            rand_control_diff = rand_control_diff, hist_control_data = hist_control_data, censor_value = censor_value,
            alpha = alpha, get_var = get_var, get_bias = get_bias, get_mse = get_mse, quietly=quietly)

    } else if (tolower(outcome_type) == "pwe") {
        pwe_error_checks(effect_vals, hist_control_data, rand_control_diff, time_vec, censor_value, alpha)
        results <- pwe_sim(trial_reps = trial_reps, subj_per_arm = subj_per_arm, a0_vals = a0_vals, effect_vals = effect_vals,
            rand_control_diff = rand_control_diff, hist_control_data = hist_control_data, time_vec_val = time_vec, censor_value = censor_value,
            alpha = alpha, get_var = get_var, get_bias = get_bias, get_mse = get_mse, quietly=quietly)

    } else if (tolower(outcome_type) == "gaussian") {
        gaussian_error_checks(effect_vals, hist_control_data, rand_control_diff, alpha)
        results <- gaussian_sim(trial_reps = trial_reps, subj_per_arm = subj_per_arm, a0_vals = a0_vals, effect_vals = effect_vals,
            rand_control_diff = rand_control_diff, hist_control_data = hist_control_data, alpha = alpha, get_var = get_var,
            get_bias = get_bias, get_mse = get_mse, quietly=quietly)

    } else if (tolower(outcome_type) == "bernoulli") {
        bernoulli_error_checks(effect_vals, hist_control_data, rand_control_diff, alpha)
        results <- bernoulli_sim(trial_reps = trial_reps, subj_per_arm = subj_per_arm, a0_vals = a0_vals, effect_vals = effect_vals,
            rand_control_diff = rand_control_diff, hist_control_data = hist_control_data, alpha = alpha, get_var = get_var,
            get_bias = get_bias, get_mse = get_mse, quietly=quietly)

    } else if (tolower(outcome_type) == "poisson") {
        poisson_error_checks(effect_vals, hist_control_data, rand_control_diff, alpha)
        results <- poisson_sim(trial_reps = trial_reps, subj_per_arm = subj_per_arm, a0_vals = a0_vals, effect_vals = effect_vals,
            rand_control_diff = rand_control_diff, hist_control_data = hist_control_data, alpha = alpha, get_var = get_var,
            get_bias = get_bias, get_mse = get_mse, quietly=quietly)
    }

    if ((length(subj_per_arm) == 1 & length(a0_vals) == 1 & length(effect_vals) == 1) |
       (length(subj_per_arm) == 1 & length(a0_vals) == 1 & length(rand_control_diff) == 1) |
       (length(subj_per_arm) == 1 & length(effect_vals) == 1 & length(rand_control_diff) == 1) |
       (length(a0_vals) == 1 & length(effect_vals) == 1 & length(rand_control_diff) == 1)){
         results$objtype <- 'realsimple'
    }
    results
}

#' Two Arm Bayesian Clinical Trial Simulation without Historical Data
#'
#' \code{simple_sim()} returns an S3 object of class \code{bayes_ctd_array}, which
#' will contain simulation results for power, statistic estimation, bias, variance,
#' and mse as requested by user.
#'
#' The object \code{bayes_ctd_array} has 6 elements: a list containing simulation
#' results (\code{data}), copies of the 4 function arguments \code{subj_per_arm},
#' \code{a0_vals}, \code{effect_vals}, and \code{rand_control_diff}, and finally
#' a \code{objtype} value indicating that \code{simple_sim()} was used. Each element of
#' \code{data} is a four-dimensional array, where each dimension is determined by the
#' length of parameters \code{subj_per_arm}, \code{a0_vals}, \code{effect_vals}, and
#' \code{rand_control_diff}. The size of \code{data} depends on which results are
#' requested by the user. At a minimum, at least one of \code{subj_per_arm},
#' \code{a0_vals}, \code{effect_vals}, or \code{rand_control_diff} must contain at
#' least 2 values, while the other three must contain at least 1 value.  The \code{data}
#' list will always contain two elements: an array of power results (\code{power}) and
#' an array of estimation results (\code{est}).  In addition to \code{power} and
#' \code{est}, \code{data} may also contain elements \code{var}, \code{bias}, or
#' \code{mse}, depending on the values of \code{get_var}, \code{get_bias}, and
#' \code{get_mse}. The values returned in \code{est} are in the form of hazard ratios,
#' mean ratios, odds ratios, or mean differences depending on the value of
#' \code{outcome_type}.   For a Gaussian outcome, the estimation results are
#' differences in group means (experimental group minus control group). For a
#' logistic outcome, the estimation results are odds ratios (experimental group over
#' control group). For lognormal and Poisson outcomes, the estimation results are mean
#' ratios (experimental group over control group). For a piecewise exponential or a
#' Weibull outcome, the estimation results are hazard ratios (experimental group over
#' control group).  The values returned in \code{bias}, \code{var}, and \code{mse} are
#' on the scale of the values returned in \code{est}.
#'
#' The object \code{bayes_ctd_array} has two primary methods, \code{print()} and
#' \code{plot()}, for printing and plotting slices of the arrays contained in
#' \code{bayes_ctd_array$data}.
#'
#' As dimensions of the four dimensional array increases, the time required to complete
#' the simulation will increase; however, it will be faster than a similar simulation
#' based on repeated calls to MCMC routines to analyze each simulated trial.
#'
#' The meaning of the estimation results, and the test used to generate power results,
#' depends on the outcome used. In all cases, power is based on a two-sided test
#' involving a (1-alpha)100\% credible interval, where the interval is used to determine
#' if the null hypothesis should be rejected (null value outside of the interval) or
#' not rejected (null value inside the interval). For a Gaussian outcome, the 95\%
#' credible interval is an interval for the difference in group means
#' (experimental group minus control group), and the test determines if 0 is in or
#' outside of the interval. For a Bernoulli outcome, the 95\% credible interval
#' is an interval for the odds ratio (experimental group over control group),
#' and the test determines if 1 is in or outside of the interval. For a lognormal or
#' a Poisson outcome, the 95\% credible interval is an interval for the mean ratio
#' (experimental group over control group), and the test determines if 1 is in or
#' outside of the interval. Finally, for a piecewise exponential or a Weibull outcome,
#' the 95\% credible interval is an interval for the hazard ratio (experimental group
#' over control group), and the test determines if 1 is in or outside of the interval.
#'
#' For a Gaussian outcome, the \code{control_parms} values should be \code{(mean, sd)},
#' where mean is the mean parameter for the control group used in a call to \code{rnorm()},
#' and sd is the common sd parameter for both groups used in a call to\code{rlnorm()}.
#'
#' For a Bernoulli outcome, the \code{control_parms} values should be \code{(prob)}, where
#' prob is the event probability for the control group used in a call to \code{rbinom()}.
#'
#' For a lognormal outcome, the \code{control_parms} values should be \code{(meanlog, sdlog)},
#' where meanlog is the meanlog parameter for the control group used in a call to
#' \code{rlnorm()}, and sdlog is the common sdlog parameter for both groups used in
#' a call to \code{rlnorm()}.
#'
#' For a Poisson outcome, the \code{control_parms} value should be \code{(lambda)}, where
#' lambda is the lambda parameter for the control group used in a call to \code{rpois()} and
#' is equal to the mean of a Poisson distribution.
#'
#' For a Weibull outcome, the \code{control_parms} values should be \code{(scale, shape)},
#' where scale is the scale parameter for the control group used in a call to
#' \code{rweibull()}, and shape is the common shape parameter for both groups used in
#' a call to \code{rweibull()}.
#'
#' For a piecewise exponential outcome, the \code{control_parms} values should be a vector
#' of lambdas used in a call to \code{eha::rpch()}.  Each element in \code{control_parms}
#' is a hazard for an interval defined by the \code{time_vec} parameter.
#'
#' Please refer to the examples for illustration of package use.
#'
#' @param trial_reps Number of trials to replicate within each combination of
#'   \code{a0_vals}, \code{subj_per_arm}, \code{effect_vals}, and \code{rand_control_parms}.
#'   As the number of trials increases, the precision of the estimate will increase.
#'   Default is 100.
#' @param outcome_type Outcome distribution. Must be equal to \code{weibull},
#'   \code{lognormal}, \code{pwe} (Piecewise Exponential), \code{gaussian},
#'   \code{bernoulli}, or \code{poisson}.  Default is \code{weibull}.
#' @param subj_per_arm A vector of sample sizes, all of which must be positive
#'   integers.  Default is \code{c(50, 100, 150, 200, 250)}.
#' @param effect_vals A vector of effects that should be reasonable for the
#'   outcome_type being studied, hazard ratios for Weibull, odds ratios for
#'   Bernoulli, mean ratios for Poisson, etc..  When \code{effect_vals} contain
#'   the null effect for a given \code{outcome_type}, the \code{power} component
#'   of \code{data} will contain an estimate of Type One Error.  In order to
#'   have a good set of Type One Error estimates, \code{trial_reps} need to be
#'   at least 10,000.  In such a case, if the total number of combinations
#'   made up from \code{subj_per_arm}, \code{a0_vals}, \code{effect_vals}, and
#'   \code{rand_control_diff} is very large, the time to complete the simulation
#'   can be substantial. Default is \code{c(0.6, 1, 1.4)}.
#' @param control_parms A vector of parameter values defining the outcome
#'   distribution for randomized controls. See Details for what is required for
#'   each \code{outcome_type}.
#' @param time_vec  A vector of time values that are used to create time periods
#'   within which the exponential hazard is constant.  Only used for piecewise
#'   exponential models.  Default is \code{NULL}.
#' @param censor_value A single value at which right censoring occurs when
#'   simulating randomized subject outcomes.  Used with survival outcomes.
#'   Default is \code{NULL}, where \code{NULL} implies no right censoring.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#' @param get_var A TRUE/FALSE indicator of whether an array of variance
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param get_bias A TRUE/FALSE indicator of whether an array of bias
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param get_mse A TRUE/FALSE indicator of whether an array of MSE
#'   estimates will be returned.  Default is \code{FALSE}.
#' @param seedval A seed value for pseudo-random number generation.
#' @param quietly A TRUE/FALSE indicator of whether notes are printed
#'   to output about simulation progress as the simulation runs.  If
#'   running interactively in RStudio or running in the R console,
#'   \code{quietly} can be set to FALSE.  If running in a Notebook or
#'   knitr document, \code{quietly} needs to be set to TRUE.  Otherwise
#'   each note will be printed on a separate line and it will take up
#'   a lot of output space.  Default is \code{TRUE}.
#'
#' @return \code{simple_sim()} returns an S3 object of class \code{bayes_ctd_array}.
#'   As noted in Details, an object of class \code{bayes_ctd_array} has 6 elements: a
#'   list containing simulation results (\code{data}), copies of the 4 function
#'   arguments \code{subj_per_arm}, \code{a0_vals}, \code{effect_vals}, and
#'   \code{rand_control_diff}, and finally \code{objtype} indicating that \code{simple_sim()}
#'   was used. See Details for a discussion about the contents of
#'   \code{data}. Results from the simulation contained in the \code{bayes_ctd_array}
#'   object can be printed or plotted using the \code{print()} and
#'   \code{plot()} methods. The results can also be accessed using basic list
#'   element identification and array slicing. For example, to get the power results
#'   from a simulation, one could use the code \code{bayes_ctd_array$data$power}, where
#'   \code{bayes_ctd_array} is replaced with the name of the variable containing the
#'   \code{bayes_ctd_array} object. Even though this is a 4-dimensional array, the power
#'   results only occupy a single 2-dimensional table. To print this 2-dimensional table,
#'   one would use the code \code{bayes_ctd_array$data$power[,1,,1]}, where
#'   \code{bayes_ctd_array} is replaced with the name of the variable containing the
#'   \code{bayes_ctd_array} object.
#'
#' @examples
#' #Run a Weibull simulation, using simple_sim().
#' #For meaningful results, trial_reps needs to be much larger than 2.
#' weibull_test <- simple_sim(trial_reps = 2, outcome_type = "weibull",
#'                            subj_per_arm = c(50, 100, 150, 200),
#'                            effect_vals = c(0.6, 1, 1.4),
#'                            control_parms = c(2.82487,3), time_vec = NULL,
#'                            censor_value = NULL, alpha = 0.05,
#'                            get_var = TRUE, get_bias = TRUE, get_mse = TRUE,
#'                            seedval=123, quietly=TRUE)
#'
#' #Tabulate the simulation results for power.
#' test_table <- print(x=weibull_test, measure="power",
#'                     tab_type=NULL, subj_per_arm_val=NULL, a0_val=NULL,
#'                     effect_val=NULL, rand_control_diff_val=NULL)
#' print(test_table)
#'
#' #Create a plot of the power simulation results.
#' plot(x=weibull_test, measure="power", tab_type=NULL,
#'      smooth=FALSE, plot_out=TRUE)
#' #Create a plot of the estimated hazard ratio simulation results.
#' plot(x=weibull_test, measure="est", tab_type=NULL,
#'      smooth=FALSE, plot_out=TRUE)
#' #Create a plot of the hazard ratio variance simulation results.
#' plot(x=weibull_test, measure="var", tab_type=NULL,
#'      smooth=FALSE, plot_out=TRUE)
#' #Create a plot of the hazard ratio bias simulation results.
#' plot(x=weibull_test, measure="bias", tab_type=NULL,
#'      smooth=FALSE, plot_out=TRUE)
#' #Create a plot of the hazard ratio mse simulation results.
#' plot(x=weibull_test, measure="mse", tab_type=NULL,
#'      smooth=FALSE, plot_out=TRUE)
#' @export
simple_sim <- function(trial_reps = 100, outcome_type = "weibull", subj_per_arm = c(50, 100, 150, 200, 250), effect_vals = c(0.6,1,
    1.4), control_parms = NULL, time_vec = NULL, censor_value = NULL, alpha = 0.05, get_var = FALSE, get_bias = FALSE,
    get_mse = FALSE, seedval=NULL, quietly=TRUE) {

	#set random seed
	set.seed(seedval)
    #------------- Go through all the high level checks for proper input. -------------#
    global_error_checks_simple(outcome_type, subj_per_arm, get_var, get_bias, get_mse)

    #-- Given outcome type go through low level checks for proper input and run simulation --#
    if (tolower(outcome_type) == "weibull") {
        weibull_error_checks_simple(effect_vals, control_parms, censor_value, alpha)
        results <- simple_weibull_sim(trial_reps = trial_reps, subj_per_arm = subj_per_arm, effect_vals = effect_vals,
            scale1_value = control_parms[1], common_shape_value = control_parms[2], censor_value = censor_value, alpha = alpha,
            get_var = get_var, get_bias = get_bias, get_mse = get_mse, quietly=quietly)

    } else if (tolower(outcome_type) == "lognormal") {
        lognormal_error_checks_simple(effect_vals, control_parms, censor_value, alpha)
        results <- simple_lognormal_sim(trial_reps = trial_reps, subj_per_arm = subj_per_arm, effect_vals = effect_vals,
            mu1_val = control_parms[1], common_sd_val = control_parms[2], censor_value = censor_value, alpha = alpha,
            get_var = get_var, get_bias = get_bias, get_mse = get_mse, quietly=quietly)

    } else if (tolower(outcome_type) == "pwe") {
        pwe_error_checks_simple(effect_vals, time_vec, control_parms, censor_value, alpha)
        results <- simple_pwe_sim(trial_reps = trial_reps, subj_per_arm = subj_per_arm, effect_vals = effect_vals, time_vec_val = time_vec,
            rc_hazards = control_parms, censor_value = censor_value, alpha = alpha, get_var = get_var, get_bias = get_bias,
            get_mse = get_mse, quietly=quietly)

    } else if (tolower(outcome_type) == "gaussian") {
        gaussian_error_checks_simple(effect_vals, control_parms, alpha)
        results <- simple_gaussian_sim(trial_reps = trial_reps, subj_per_arm = subj_per_arm, effect_vals = effect_vals,
            mu1_val = control_parms[1], common_sd_val = control_parms[2], alpha = alpha, get_var = get_var, get_bias = get_bias,
            get_mse = get_mse, quietly=quietly)

    } else if (tolower(outcome_type) == "bernoulli") {
        bernoulli_error_checks_simple(effect_vals, control_parms, alpha)
        results <- simple_bernoulli_sim(trial_reps = trial_reps, subj_per_arm = subj_per_arm, effect_vals = effect_vals,
            prob1_val = control_parms[1], alpha = alpha, get_var = get_var, get_bias = get_bias, get_mse = get_mse, quietly=quietly)

    } else if (tolower(outcome_type) == "poisson") {
        poisson_error_checks_simple(effect_vals, control_parms, alpha)
        results <- simple_poisson_sim(trial_reps = trial_reps, subj_per_arm = subj_per_arm, effect_vals = effect_vals,
            mu1_val = control_parms[1], alpha = alpha, get_var = get_var, get_bias = get_bias, get_mse = get_mse, quietly=quietly)
    }
}



#' Print Data from Two Arm Bayesian Clinical Trial Simulation.
#'
#' \code{print.bayes_ctd_array()} takes an S3 object of class \code{bayes_ctd_array}, and
#' prints a two dimensional slice from the data generated by a clinical trial simulation
#' using \code{historic_sim()} or \code{simple_sim()}.
#'
#' If the object of class \code{bayes_ctd_array} is created by \code{historic_sim()},
#' then the function \code{print()} allows the user to print user-specified 1- and 2-
#' dimensional slices of the simulation results based on slicing code described
#' below.  If the object of class \code{bayes_ctd_array} is created by
#' \code{simple_sim()}, a basic table of characteristic by sample size and effect is created.
#'
#' If the object of class \code{bayes_ctd_array} is created by \code{simple_sim()}, then
#' all four trial characteristics (\code{subj_per_arm_val}, \code{a0_vals},
#' \code{effect_val}, and \code{rand_control_diff_val}) can be ignored, as can the
#' parameter defining what type of table to print, \code{tab_type}.  A call to
#' \code{print()} will require the user to specify a measure (power, est, var, bias,
#' or mse).
#'
#' If the object of class \code{bayes_ctd_array} is created by \code{historic_sim()},
#' a call to \code{print()} will require the user to specify a measure
#' (power, est, var, bias, or mse) and may require the user to specify a table type.
#' A table type, \code{tab_type}, will be required if 3 of the 4 trial characteristics
#' are equal to a vector of 2 or more values.  The table type specification
#' uses the letters W, X, Y, and Z.  The letter W represents the subject per arm
#' dimension.  The letter X represents the a0 dimension.  The letter Y represents
#' the effect dimension.  The letter Z represents the control difference dimension.
#' To define a slice of the 4-dimensional array, these letters are put into an AB|CD
#' pattern.  The two letters to the right of the vertical bar define which variables
#' are held constant.  The two letters to the left of the vertical bar define which
#' variables are going to show up in the rows (first letter) and in the columns (second
#' letter).  For example if tab_type equals \code{WX|YZ}, then effect and control
#' differences will be held constant, while sample size will be represented by the rows
#' in the generated table and a0 values will be represented by the columns.  The actual
#' values that are printed in the tables depend on what measure is requested in the
#' parameter \code{measure}.
#'
#' \itemize{
#'   \item \code{tab_type='WX|YZ'}, Sample Size by a0
#'   \item \code{tab_type='WY|XZ'}, Sample Size by Effect
#'   \item \code{tab_type='WZ|XY'}, Sample Size by Control Differences
#'   \item \code{tab_type='XY|WZ'}, a0 by Effect
#'   \item \code{tab_type='XZ|WY'}, a0 by Control Differences
#'   \item \code{tab_type='YZ|WX'}, Effect by Control Differences
#'   \item \code{tab_type='ZX|WY'}, Control Differences by a0
#'   \item \code{tab_type='XW|YZ'}, a0 by Sample Size
#'   \item \code{tab_type='YW|XZ'}, Effect by Sample Size
#'   \item \code{tab_type='YX|WZ'}, Effect by a0
#'   \item \code{tab_type='ZW|XY'}, Control Differences by Sample Size
#'   \item \code{tab_type='ZY|WX'}, Control Differences by Effect
#' }
#'
#' It is very important to populate the values of \code{subj_per_arm_val},
#' \code{a0_vals}, \code{effect_val}, and \code{rand_control_diff_val} correctly given
#' the value of tab_type, when the object of class \code{bayes_ctd_array} is created by
#' \code{historic_sim()} and at least 3 of the four parameters have more than one
#' value.  On the other hand, if 2 or more of the four parameters have only one value,
#' then \code{subj_per_arm_val}, \code{a0_vals}, \code{effect_val},
#' \code{rand_control_diff_val}, as well as \code{tab_type} can be ignored.  If the last
#' two letters are \code{YZ}, then \code{effect_val} and \code{rand_control_diff_val}
#' must be populated.  If the last two letters are \code{XZ}, then \code{a0_vals} and
#' \code{rand_control_diff_val} must be populated.  If the last two letters are
#' \code{XY}, then \code{a0_vals} and \code{effect_val} must be populated.  If the last
#' two letters are \code{WZ}, then \code{sample_val} and \code{rand_control_diff_val}
#' must be populated.  If the last two letters are \code{WY}, then \code{sample_size_val}
#' and \code{effect_val} must be populated.  If the last two letters are \code{WX}, then
#' \code{sample_size_val} and \code{a0_vals} must be populated.
#'
#' If the object of class \code{bayes_ctd_array} is created by \code{simple_sim()}, the
#' parameters \code{tab_type}, \code{subj_per_arm_val}, \code{a0_vals}, \code{effect_val},
#' and \code{rand_control_diff_val} are ignored.
#'
#' @param x Name of object of class \code{bayes_ctd_array} containing
#'   data from clinical trial simulation.
#' @param measure Must be equal to \code{power}, \code{est}, \code{var}, \code{bias},
#'   or \code{mse}.  Default is \code{power}.  Case does not matter.
#' @param tab_type A character string that must equal \code{WX|YZ}, \code{WY|XZ},
#'   \code{WZ|XY}, \code{XY|WZ}, \code{XZ|WY}, \code{YZ|WX}, \code{ZX|WY}, \code{XW|YZ},
#'   \code{YW|XZ}, \code{YX|WZ}, \code{ZW|XY}, \code{ZX|WY}, \code{ZY|WX} when
#'   \code{x} is generated by \code{historic_sim()}.  Default is
#'   \code{WX|YZ}.  When \code{x} is generated by \code{simple_sim()},
#'   \code{tab_type} is ignored.
#' @param subj_per_arm_val Must be non-missing, if \code{x} is generated
#'   by \code{historic_sim()} and sample size is being held constant.
#'   If \code{x} is generated by \code{historic_sim()} and sample size
#'   is being held constant, \code{subj_per_arm_val} must equal a value submitted
#'   to \code{historic_sim()} within the \code{subj_per_arm} parameter.  When
#'   \code{x} is generated by \code{simple_sim()}, \code{subj_per_arm_val}
#'   is ignored.
#' @param a0_val Must be non-missing, if \code{x} is generated
#'   by \code{historic_sim()} and a0, the power prior parameter, is being held
#'   constant.  If \code{x} is generated by \code{historic_sim()} and
#'   a0 is being held constant, \code{a0_val} must equal a value submitted
#'   to \code{historic_sim()} within the \code{a0_vals} parameter.  When
#'   \code{x} is generated by \code{simple_sim()}, \code{a0_val} is
#'   ignored.
#' @param effect_val  Must be non-missing, if \code{x} is generated
#'   by \code{historic_sim()} and effect is being held constant.  If
#'   \code{x} is generated by \code{historic_sim()} and effect is being
#'   held constant, \code{effect_val} must equal a value submitted to
#'   \code{historic_sim()} within the \code{effect_vals} parameter.  When
#'   \code{x} is generated by \code{simple_sim()}, \code{effect_val} is
#'   ignored.
#' @param rand_control_diff_val Must be non-missing, if \code{x} is
#'   generated by \code{historic_sim()} and differences between randomized
#'   and historical controls are being held constant.  If \code{x}
#'   is generated by \code{historic_sim()} and control differences are being
#'   held constant, \code{rand_control_diff_val} must equal a value submitted to
#'   \code{historic_sim()} within the \code{rand_control_diff} parameter.  When
#'   \code{x} is generated by \code{simple_sim()},
#'   \code{rand_control_diff_val} is ignored.
#' @param print_chg_warn A parameter not used by the user, but is used by
#'   \code{plot()} to ensure warnings are not printed twice.
#' @param ...	further arguments passed to or from other methods.
#'
#' @return \code{print()} returns a two dimensional array of simulation results.
#'
#' @examples
#' #Run a Weibull simulation, using simple_sim().
#' #For meaningful results, trial_reps needs to be much larger than 2.
#' weibull_test <- simple_sim(trial_reps = 2, outcome_type = "weibull",
#'                            subj_per_arm = c(50, 100, 150, 200),
#'                            effect_vals = c(0.6, 1, 1.4),
#'                            control_parms = c(2.82487,3),
#'                            time_vec = NULL, censor_value = NULL,
#'                            alpha = 0.05, get_var = TRUE,
#'                            get_bias = TRUE, get_mse = TRUE,
#'                            seedval=123, quietly=TRUE)
#'
#' #Tabulate the simulation results for power.
#' test_table <- print(x=weibull_test, measure="power",
#'                     tab_type=NULL, subj_per_arm_val=NULL, a0_val=NULL,
#'                     effect_val=NULL, rand_control_diff_val=NULL)
#' print(test_table)
#'
#' #Tabulate the simulation results for estimates.
#' print(x=weibull_test, measure="est")
#'
#' #Tabulate the simulation results for variance.
#' print(x=weibull_test, measure="var")
#'
#' #Tabulate the simulation results for bias.
#' print(x=weibull_test, measure="bias")
#'
#' #Tabulate the simulation results for mse.
#' print(x=weibull_test, measure="mse")
#'
#' \donttest{
#' #Run another weibull simulation, using historic_sim().
#' #Note: historic_sim() can take a while to run.
#' #Generate a sample of historical data for use in example.
#' set.seed(2250)
#' SampleHistData <- genweibulldata(sample_size=60, scale1=2.82487,
#'                                  hazard_ratio=0.6, common_shape=3,
#'                                  censor_value=3)
#' histdata <- subset(SampleHistData, subset=(treatment==0))
#' histdata$id <- histdata$id+10000
#'
#' #For meaningful results, trial_reps needs to be larger than 100.
#' weibull_test2 <- historic_sim(trial_reps = 100, outcome_type = "weibull",
#'                               subj_per_arm = c(50, 100, 150, 200, 250),
#'                               a0_vals = c(0, 0.33, 0.67, 1),
#'                               effect_vals = c(0.6, 1, 1.4),
#'                               rand_control_diff = c(0.8, 1, 1.2),
#'                               hist_control_data = histdata, time_vec = NULL,
#'                               censor_value = 3, alpha = 0.05, get_var = TRUE,
#'                               get_bias = TRUE, get_mse = TRUE, seedval=123,
#'                               quietly=TRUE)
#'
#' #Tabulate the simulation results for power.
#' test_table <- print(x=weibull_test2, measure="power",
#'                     tab_type="WX|YZ", effect_val=0.6,
#'                     rand_control_diff_val=1.0)
#' print(test_table)
#'
#' #Tabulate the simulation results for estimates.
#' print(x=weibull_test2, measure="est", tab_type="WX|YZ",
#'       effect_val=0.6, rand_control_diff_val=1.0)
#'
#' #Tabulate the simulation results for variance.
#' print(x=weibull_test2, measure="var", tab_type="WX|YZ",
#'       effect_val=0.6, rand_control_diff_val=1.0)
#'
#' #Tabulate the simulation results for bias.
#' print(x=weibull_test2, measure="bias", tab_type="WX|YZ",
#'       effect_val=0.6, rand_control_diff_val=1.0)
#'
#' #Tabulate the simulation results for mse.
#' print(x=weibull_test2, measure="mse", tab_type="WX|YZ",
#'       effect_val=0.6, rand_control_diff_val=1.0)
#' }
#'
#' \donttest{
#' #Run a Bernoulli simulation, using historic_sim().
#' #Generate a sample of historical Bernoulli data for use in example.
#' set.seed(2250)
#' samplehistdata <- genbernoullidata(sample_size=60, prob1=0.6, odds_ratio=0.6)
#' histdata <- subset(samplehistdata, subset=(treatment==0))
#' histdata$id <- histdata$id+10000
#'
#' #For meaningful results, trial_reps needs to be larger than 100.
#' bernoulli_test <- historic_sim(trial_reps = 100, outcome_type = "bernoulli",
#'                               subj_per_arm = c(150),
#'                               a0_vals = c(1.0),
#'                               effect_vals = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),
#'                               rand_control_diff = c(1.8),
#'                               hist_control_data = histdata, time_vec = NULL,
#'                               censor_value = 3, alpha = 0.05, get_var = TRUE,
#'                               get_bias = TRUE, get_mse = TRUE, seedval=123,
#'                               quietly=TRUE)
#' test_table <- print(x=bernoulli_test, measure="power",
#'                     tab_type=NULL, effect_val=NULL,
#'                     subj_per_arm_val=NULL)
#' print(test_table)
#'
#' #If only one or two of the subj_per_arm, a0_vals, effect_vals, or
#' #rand_control_diff parameters have length greater than 1, then
#' #only bayes_ctd_array and measure parameters are needed.
#' #Tabulate the simulation results for estimates.
#' print(x=bernoulli_test, measure="est")
#'
#' #Tabulate the simulation results for variance.
#' print(x=bernoulli_test, measure="var")
#'
#' #Tabulate the simulation results for bias.
#' print(x=bernoulli_test, measure="bias")
#'
#' #Tabulate the simulation results for mse.
#' print(x=bernoulli_test, measure="mse")
#' }
#'
#' @export
print.bayes_ctd_array <- function(x = NULL, measure = "power", tab_type = "WX|YZ", subj_per_arm_val = NULL,
                                        a0_val = NULL, effect_val = NULL, rand_control_diff_val = NULL, print_chg_warn = 1, ...) {

  # (W): subj_per_arm,
  # (X): a0_val,
  # (Y): effect_val,
  # (Z): rand_control_diff
  # tab_type='WX|YZ', Table of Design Characteristics: Sample Size by a0
  # tab_type='WY|XZ', Table of Design Characteristics: Sample Size by Effect
  # tab_type='WZ|XY', Table of Design Characteristics: Sample Size by Control Differences
  # tab_type='XY|WZ', Table of Design Characteristics: a0 by Effect
  # tab_type='XZ|WY', Table of Design Characteristics: a0 by Control Diffferences
  # tab_type='YZ|WX', Table of Design Characteristics: Effect by Control Diffferences
  # tab_type='ZX|WY', Table of Design Characteristics: Control Diffferences by a0
  # tab_type='XW|YZ', Table of Design Characteristics: a0 by Sample Size
  # tab_type='YW|XZ', Table of Design Characteristics: Effect by Sample Size
  # tab_type='YX|WZ', Table of Design Characteristics: Effect by a0
  # tab_type='ZW|XY', Table of Design Characteristics: Control Differences by Sample Size
  # tab_type='ZY|WX', Table of Design Characteristics: Control Differences by Effect

  #Add simple code to help out the user with tab_type on matters where order is
  #not important and can be easily fixed.
  tab_type = toupper(tab_type)
  if (!((length(tab_type)== 0) && (typeof(tab_type) == "character"))){
    if (x$objtype == 'historic') {
      if (tab_type=='WX|ZY') tab_type <- 'WX|YZ'
      if (tab_type=='WY|ZX') tab_type <- 'WY|XZ'
      if (tab_type=='WZ|YX') tab_type <- 'WZ|XY'
      if (tab_type=='XY|ZW') tab_type <- 'XY|WZ'
      if (tab_type=='XZ|YW') tab_type <- 'XZ|WY'
      if (tab_type=='YZ|XW') tab_type <- 'YZ|WX'
      if (tab_type=='ZX|YW') tab_type <- 'ZX|WY'
      if (tab_type=='XW|ZY') tab_type <- 'XW|YZ'
      if (tab_type=='YW|ZX') tab_type <- 'YW|XZ'
      if (tab_type=='YX|ZW') tab_type <- 'YX|WZ'
      if (tab_type=='ZW|YX') tab_type <- 'ZW|XY'
      if (tab_type=='ZY|XW') tab_type <- 'ZY|WX'
    }
  }
  #Do the same for measure
  measure = tolower(measure)
  #The above code was added very late in initial production, so time was not
  # available to go through all code and remove code no longer necessary.


  #------------- Go through all the high level checks for proper input. -------------#

  if (x$objtype == 'simple') {
    tab_type <- "WY|XZ"
    a0_val <- 0
    rand_control_diff_val <- 1
    subj_per_arm_val <- x$subj_per_arm[1]
    effect_val <- x$effect_vals[1]
    if (print_chg_warn == 1) {
      print("Since simple_sim was used, tab_type was set to WY|XZ")
      print("Values for tab_type, subj_per_arm_val, a0_val, effect_val, and rand_control_diff_val were ignored")
      print("This works towards putting all results in a single table, effect by sample size")
    }
  }

  if (x$objtype == 'realsimple') {
    if (length(x$subj_per_arm) > 1) {
      if (print_chg_warn == 1) {
        print("Since only subj_per_arm vector has more than 1 element, tab_type was set to WX|YZ")
        print("This works towards putting all results in a single table")
      }
      tab_type <- "WX|YZ"
      a0_val <- x$a0_vals[1]
      effect_val <- x$effect_vals[1]
      rand_control_diff_val <- x$rand_control_diff[1]
    }
    if (length(x$a0_vals) > 1) {
      if (print_chg_warn == 1) {
        print("Since only a0_vals vector has more than 1 element, tab_type was set to XY|WZ")
        print("This works towards putting all results in a single table")
      }
      tab_type <- "XY|WZ"
      subj_per_arm_val <- x$subj_per_arm[1]
      effect_val <- x$effect_vals[1]
      rand_control_diff_val <- x$rand_control_diff[1]
    }
    if (length(x$effect_vals) > 1) {
      if (print_chg_warn == 1) {
        print("Since only effect_vals vector has more than 1 element, tab_type was set to YZ|WX")
        print("This works towards putting all results in a single table")
      }
      tab_type <- "YZ|WX"
      subj_per_arm_val <- x$subj_per_arm[1]
      a0_val <- x$a0_vals[1]
      rand_control_diff_val <- x$rand_control_diff[1]
    }
    if (length(x$rand_control_diff) > 1) {
      if (print_chg_warn == 1) {
        print("Since only rand_control_diff vector has more than 1 element, tab_type was set to ZX|WY")
        print("This works towards putting all results in a single table")
      }
      tab_type <- "ZX|WY"
      subj_per_arm_val <- x$subj_per_arm[1]
      a0_val <- x$a0_vals[1]
      effect_val <- x$effect_vals[1]
    }
  }

  # Check to see if two dimensions of bayes_ctd_array have only 1 element each. Is so set tab_type to an appropriate
  # value.
  if (x$objtype == 'historic') {
    if (length(x$effect_vals) == 1 & length(x$rand_control_diff) == 1) {
      chkflag1 <- 0
      if(!((length(tab_type)== 0) && (typeof(tab_type) == "character"))){
        if ((tab_type != "WX|YZ") & (tab_type != "XW|YZ")) {
          tab_type <- "WX|YZ"
          chkflag1 <- 1
        }
      }
      if( (length(tab_type)== 0) && (typeof(tab_type) == "character") ){
        tab_type <- "WX|YZ"
        chkflag1 <- 1
      }
      effect_val <- x$effect_vals[1]
      rand_control_diff_val <- x$rand_control_diff[1]
      if (print_chg_warn == 1 & chkflag1 == 1) {
        print("Since effect_vals and rand_control_diff vectors only have 1 element each and tab_type not equal to 'WX|YZ' or 'XW|YZ', tab_type was set to WX|YZ")
        print("This works towards putting all results in a single table")
      }
    }
    if (length(x$a0_vals) == 1 & length(x$rand_control_diff) == 1) {
      chkflag1 <- 0
      if(!((length(tab_type)== 0) && (typeof(tab_type) == "character"))){
        if ((tab_type != "WY|XZ") & (tab_type != "YW|XZ")) {
          tab_type <- "WY|XZ"
          chkflag1 <- 1
        }
      }
      if((length(tab_type)== 0) && (typeof(tab_type) == "character")){
        tab_type <- "WY|XZ"
        chkflag1 <- 1
      }
      a0_val <- x$a0_vals[1]
      rand_control_diff_val <- x$rand_control_diff[1]
      if (print_chg_warn == 1 & chkflag1 == 1) {
        print("Since a0 and rand_control_diff vectors only have 1 element each and tab_type not equal to 'WY|XZ' or 'YW|XZ', tab_type was set to WY|XZ")
        print("This works towards putting all results in a single table")
      }
    }
    if (length(x$subj_per_arm) == 1 & length(x$rand_control_diff) == 1) {
      chkflag1 <- 0
      if(!((length(tab_type)== 0) && (typeof(tab_type) == "character"))){
        if ((tab_type != "XY|WZ") & (tab_type != "YX|WZ")) {
          tab_type <- "XY|WZ"
          chkflag1 <- 1
        }
      }
      if((length(tab_type)== 0) && (typeof(tab_type) == "character")){
        tab_type <- "XY|WZ"
        chkflag1 <- 1
      }
      subj_per_arm_val <- x$subj_per_arm
      rand_control_diff_val <- x$rand_control_diff
      if (print_chg_warn == 1 & chkflag1 == 1) {
        print("Since sample size and rand_control_diff vectors only have 1 element each and tab_type not equal to 'XY|WZ' or 'YX|WZ', tab_type was set to XY|WZ")
        print("This works towards putting all results in a single table")
      }
    }
    if (length(x$subj_per_arm) == 1 & length(x$effect_vals) == 1) {
      chkflag1 <- 0
      if(!((length(tab_type)== 0) && (typeof(tab_type) == "character"))){
        if ((tab_type != "XZ|WY") & (tab_type != "ZX|WY")) {
          tab_type <- "XZ|WY"
          chkflag1 <- 1
        }
      }
      if((length(tab_type)== 0) && (typeof(tab_type) == "character")){
        tab_type <- "XZ|WY"
        chkflag1 <- 1
      }
      subj_per_arm_val <- x$subj_per_arm
      effect_val <- x$effect_vals
      if (print_chg_warn == 1 & chkflag1 == 1) {
        print("Since sample size and effect_vals vectors only have 1 element each and tab_type not equal to 'XZ|WY' or 'ZX|WY', tab_type was set to XZ|WY")
        print("This works towards putting all results in a single table")
      }
    }
    if (length(x$subj_per_arm) == 1 & length(x$a0_vals) == 1) {
      chkflag1 <- 0
      if(!((length(tab_type)== 0) && (typeof(tab_type) == "character"))){
        if ((tab_type != "YZ|WX") & (tab_type != "ZY|WX")) {
          tab_type <- "YZ|WX"
          chkflag1 <- 1
        }
      }
      if((length(tab_type)== 0) && (typeof(tab_type) == "character")){
        tab_type <- "YZ|WX"
        chkflag1 <- 1
      }
      subj_per_arm_val <- x$subj_per_arm
      a0_val <- x$a0_vals
      if (print_chg_warn == 1 & chkflag1 == 1) {
        print("Since sample size and a0 vectors only have 1 element each and tab_type not equal to 'YZ|WX' or 'ZY|WX', tab_type was set to YZ|WX")
        print("This works towards putting all results in a single table")
      }
    }
    if (length(x$a0_vals) == 1 & length(x$effect_vals) == 1) {
      chkflag1 <- 0
      if(!((length(tab_type)== 0) && (typeof(tab_type) == "character"))){
        if ((tab_type != "WZ|XY") & (tab_type != "ZW|XY")) {
          tab_type <- "WZ|XY"
          chkflag1 <- 1
        }
      }
      if((length(tab_type)== 0) && (typeof(tab_type) == "character")){
        tab_type <- "WZ|XY"
        chkflag1 <- 1
      }
      a0_val <- x$a0_vals
      effect_val <- x$effect_vals
      if (print_chg_warn == 1 & chkflag1 == 1) {
        print("Since a0 and effect_vals vectors only have 1 element each and tab_type not equal to 'WZ|XY' or 'ZW|XY', tab_type was set to WZ|XY")
        print("This works towards putting all results in a single table")
      }
    }
  }

  print_error_checks(x, measure, tab_type, subj_per_arm_val, a0_val, effect_val, rand_control_diff_val)

  #Select the data array depending on user defined value of measure.
  if (tolower(measure) == "power")
    dataarray <- x$data$power
  if (tolower(measure) == "est")
    dataarray <- x$data$est
  if (tolower(measure) == "var")
    dataarray <- x$data$var
  if (tolower(measure) == "bias")
    dataarray <- x$data$bias
  if (tolower(measure) == "mse")
    dataarray <- x$data$mse
  if (toupper(tab_type) == "WX|YZ") {
    # Table of Design Characteristics: Sample Size by a0
    effect_val_levels <- x$effect_vals
    rand_control_diff_levels <- x$rand_control_diff
    # check to make sure effect_val is in effect_vals and
    # rand_control_diff_val is in rand_control_diff
    effect_val_loc <- match(effect_val, effect_val_levels)
    rand_control_diff_loc <- match(rand_control_diff_val, rand_control_diff_levels)
    return_table <- dataarray[,, effect_val_loc, rand_control_diff_loc]
  }
  if (toupper(tab_type) == "WY|XZ") {
    # Table of Design Characteristics: Sample Size by Effect
    a0_val_levels <- x$a0_vals
    rand_control_diff_levels <- x$rand_control_diff
    # check to make sure a0_val is in a0_vals and
    # rand_control_diff_val is in rand_control_diff
    a0_val_loc <- match(a0_val, a0_val_levels)
    rand_cntrl_diff_loc <- match(rand_control_diff_val, rand_control_diff_levels)
    return_table <- dataarray[, a0_val_loc,, rand_cntrl_diff_loc]
  }
  if (toupper(tab_type) == "WZ|XY") {
    # Table of Design Characteristics: Sample Size by Control Differences
    a0_val_levels <- x$a0_vals
    effect_val_levels <- x$effect_vals
    # check to make sure effect_val is in effect_vals and
    # a0_val is in a0_vals
    a0_val_loc <- match(a0_val, a0_val_levels)
    effect_val_loc <- match(effect_val, effect_val_levels)
    return_table <- dataarray[, a0_val_loc, effect_val_loc, ]
  }
  if (toupper(tab_type) == "XY|WZ") {
    # Table of Design Characteristics: a0 by Effect
    subj_per_arm_levels <- x$subj_per_arm
    rand_control_diff_levels <- x$rand_control_diff
    # check to make sure subj_per_arm_val is in subj_per_arm and
    # rand_control_diff_val is in rand_control_diff
    subj_per_arm_loc <- match(subj_per_arm_val, subj_per_arm_levels)
    rand_control_diff_loc <- match(rand_control_diff_val, rand_control_diff_levels)
    return_table <- dataarray[subj_per_arm_loc,,, rand_control_diff_loc]
  }
  if (toupper(tab_type) == "XZ|WY") {
    # Table of Design Characteristics: a0 by Control Diffferences
    subj_per_arm_levels <- x$subj_per_arm
    effect_val_levels <- x$effect_vals
    # check to make sure effect_val is in effect_vals and
    # subj_per_arm_val is in subj_per_arm
    subj_per_arm_loc <- match(subj_per_arm_val, subj_per_arm_levels)
    effect_val_loc <- match(effect_val, effect_val_levels)
    return_table <- dataarray[subj_per_arm_loc,, effect_val_loc, ]
  }
  if (toupper(tab_type) == "YZ|WX") {
    # Table of Design Characteristics: Effect by Control Diffferences
    subj_per_arm_levels <- x$subj_per_arm
    a0_val_levels <- x$a0_vals
    # check to make sure subj_per_arm_val is in subj_per_arm and
    # a0_val is in a0_vals
    subj_per_arm_loc <- match(subj_per_arm_val, subj_per_arm_levels)
    a0_val_loc <- match(a0_val, a0_val_levels)
    return_table <- dataarray[subj_per_arm_loc, a0_val_loc,, ]
  }
  if (toupper(tab_type) == "XW|YZ") {
    # Table of Design Characteristics: a0 by sample size
    effect_val_levels <- x$effect_vals
    rand_control_diff_levels <- x$rand_control_diff
    # check to make sure rand_control_diff_val is in rand_control_diff and
    # effect_val is in effect_vals
    effect_val_loc <- match(effect_val, effect_val_levels)
    rand_control_diff_loc <- match(rand_control_diff_val, rand_control_diff_levels)
    return_table <- t(dataarray[,,effect_val_loc,rand_control_diff_loc])
  }
  if (toupper(tab_type) == "YW|XZ") {
    # Table of Design Characteristics: Effect by sample size
    a0_val_levels <- x$a0_vals
    rand_control_diff_levels <- x$rand_control_diff
    # check to make sure rand_control_diff_val is in rand_control_diff and
    # a0_val is in a0_vals
    a0_val_loc <- match(a0_val, a0_val_levels)
    rand_control_parm_loc <- match(rand_control_diff_val, rand_control_diff_levels)
    return_table <- t(dataarray[,a0_val_loc,,rand_control_parm_loc])
  }
  if (toupper(tab_type) == "YX|WZ") {
    # Table of Design Characteristics: Effect by a0
    subj_per_arm_levels <- x$subj_per_arm
    rand_control_diff_levels <- x$rand_control_diff
    # check to make sure subj_per_arm_val is in subj_per_arm and
    # rand_control_diff_val is in rand_control_diff
    subj_per_arm_loc <- match(subj_per_arm_val, subj_per_arm_levels)
    rand_control_parm_loc <- match(rand_control_diff_val, rand_control_diff_levels)
    return_table <- dataarray[subj_per_arm_loc,,,rand_control_parm_loc]
  }
  if (toupper(tab_type) == "ZW|XY") {
    # Table of Design Characteristics: Control Differences by sample size
    a0_val_levels <- x$a0_vals
    effect_val_levels <- x$effect_vals
    # check to make sure a0_val is in a0_vals and
    # effect_val is in effect_vals
    a0_val_loc <- match(a0_val, a0_val_levels)
    effect_val_loc <- match(effect_val, effect_val_levels)
    return_table <- t(dataarray[,a0_val_loc,effect_val_loc, ])
  }
  if (toupper(tab_type) == "ZX|WY") {
    # Table of Design Characteristics: Control Differences by a0
    subj_per_arm_levels <- x$subj_per_arm
    effect_val_levels <- x$effect_vals
    # check to make sure subj_per_arm_val is in subj_per_arm and
    # effect_val is in effect_vals
    subj_per_arm_loc <- match(subj_per_arm_val, subj_per_arm_levels)
    effect_val_loc <- match(effect_val, effect_val_levels)
    return_table <- t(dataarray[subj_per_arm_loc,,effect_val_loc, ])
  }
  if (toupper(tab_type) == "ZY|WX") {
    # Table of Design Characteristics: Control Differences by Effect
    subj_per_arm_levels <- x$subj_per_arm
    a0_val_levels <- x$a0_vals
    # check to make sure subj_per_arm_val is in subj_per_arm and
    # a0_val is in a0_vals
    subj_per_arm_loc <- match(subj_per_arm_val, subj_per_arm_levels)
    a0_val_loc <- match(a0_val, a0_val_levels)
    return_table <- t(dataarray[subj_per_arm_loc,a0_val_loc,, ])
  }

  if (toupper(tab_type) == "YX|WZ") {
    # Table of Design Characteristics: Effect by a0
    subj_per_arm_levels <- x$subj_per_arm
    rand_control_diff_levels <- x$rand_control_diff
    # check to make sure subj_per_arm_val is in subj_per_arm and
    # rand_control_diff_val is in rand_control_diff
    subj_per_arm_loc <- match(subj_per_arm_val, subj_per_arm_levels)
    rand_control_diff_val_loc <- match(rand_control_diff_val, rand_control_diff_levels)
    return_table <- t(dataarray[subj_per_arm_loc,,,rand_control_diff_val_loc ])
  }

  return_table
}



#' Plot Data from Two Arm Bayesian Clinical Trial Simulation.
#'
#' \code{plot.bayes_ctd_array()} takes an S3 object of class \code{bayes_ctd_array}, and
#' creates a line plot from a one or two dimensional slice of the data generated by a
#' clinical trial simulation using \code{historic_sim()} or \code{simple_sim()}.  The
#' plotted results can be smoothed or unsmoothed.
#'
#' If the object of class \code{bayes_ctd_array} is created by \code{historic_sim()},
#' the function \code{plot()} allows the user to create line plots of user-specified
#' 1- or 2- dimensional slices of the simulation results based on slicing code
#' described below. If the object of class \code{bayes_ctd_array} is created by
#' \code{simple_sim()}, a basic plot of characteristic by sample size and effect is created.
#'
#' If the object of class \code{bayes_ctd_array} is created by \code{simple_sim()}, then
#' all four trial characteristics (\code{subj_per_arm_val}, \code{a0_vals},
#' \code{effect_val}, and \code{rand_control_diff_val}) can be ignored as can the
#' parameter defining what type of plot to create through the parameter \code{tab_type}.
#' A call to \code{plot()} will require the user to specify a measure (power, est,
#' var, bias, or mse).
#'
#' If the object of class \code{bayes_ctd_array} is created by \code{historic_sim()},
#' when calling \code{plot()} the user must specify a measure to plot
#' (power, est, var, bias, or mse) and may be required to specify a plot type through
#' the \code{tab_type} parameter.  A plot type, \code{tab_type}, will be required if
#' 3 of the 4 trial characteristics are equal to a vector of 2 or more values.  This
#' plot type specification uses the letters W, X, Y, and Z.  The letter W represents
#' the subject per arm dimension.  The letter X represents the a0 dimension.  The
#' letter Y represents the effect dimension.  The letter Z represents the control
#' difference dimension.  To plot a slice of the 4-dimensional array, these letters
#' are put into an AB|CD pattern just like in \code{print()}.  The two letters
#' to the right of the vertical bar define which variables are held constant.  The two
#' letters to the left of the vertical bar define which variables are going to show up
#' in the plot.  The first letter defines the x-axis variable and the second letter
#' defines the stratification variable.  The result is a plot of power, estimate,
#' variance, bias, or mse by the trial characteristic represented by the first letter.
#' On this plot, one line will be created for each value of the trial characteristic
#' represented by the second letter.  For example if tab_type equals \code{WX|YZ},
#' then effect and control differences will be held constant, while sample size will be
#' represented along the horizontal axis and a0 values will be represented by separate
#' lines.  The actual values that are plotted on the y-axis depend on what measure is
#' requested in the parameter \code{measure}.
#'
#' \itemize{
#'   \item \code{tab_type='WX|YZ'}, Sample Size by a0
#'   \item \code{tab_type='WY|XZ'}, Sample Size by Effect
#'   \item \code{tab_type='WZ|XY'}, Sample Size by Control Differences
#'   \item \code{tab_type='XY|WZ'}, a0 by Effect
#'   \item \code{tab_type='XZ|WY'}, a0 by Control Differences
#'   \item \code{tab_type='YZ|WX'}, Effect by Control Differences
#'   \item \code{tab_type='ZX|WY'}, Control Differences by a0
#'   \item \code{tab_type='XW|YZ'}, a0 by Sample Size
#'   \item \code{tab_type='YW|XZ'}, Effect by Sample Size
#'   \item \code{tab_type='YX|WZ'}, Effect by a0
#'   \item \code{tab_type='ZW|XY'}, Control Differences by Sample Size
#'   \item \code{tab_type='ZY|WX'}, Control Differences by Effect
#' }
#' It is very important to populate the values of \code{subj_per_arm_val},
#' \code{a0_val}, \code{effect_val}, and \code{rand_control_diff_val} correctly given
#' the value of tab_type, when the object of class \code{bayes_ctd_array} is created by
#' \code{historic_sim()} and at least 3 of the four parameters have more than one
#' value.  On, the other hand, if 2 or more of the four parameters have only one value,
#' then \code{subj_per_arm_val}, \code{a0_vals}, \code{effect_val},
#' \code{rand_control_diff_val}, as well as \code{tab_type} can be ignored.  If the last
#' two letters are \code{YZ}, then \code{effect_val} and \code{rand_control_diff_val}
#' must be populated.  If the last two letters are \code{XZ}, then \code{a0_val} and
#' \code{rand_control_diff_val} must be populated.  If the last two letters are \code{XY},
#' then \code{a0_val} and \code{effect_val} must be populated.  If the last two letters
#' are \code{WZ}, then \code{sample_val} and \code{rand_control_diff_val} must be
#' populated.  If the last two letters are \code{WY}, then \code{sample_size_val} and
#' \code{effect_val} must be populated.  If the last two letters are \code{WX}, then
#' \code{sample_size_val} and \code{a0_val} must be populated.
#'
#' If the object of class \code{bayes_ctd_array} is created by \code{simple_sim()}, the
#' parameters \code{tab_type}, \code{subj_per_arm_val}, \code{a0_val}, \code{effect_val},
#' and \code{rand_control_diff_val} are ignored.
#'
#' @param x Name of object of class \code{bayes_ctd_array} containing
#'   data from clinical trial simulation.
#' @param measure Must be equal to \code{power}, \code{est}, \code{var}, \code{bias},
#'   or \code{mse}.  Default is \code{power}.  Case does not matter.
#' @param tab_type A character string that must equal \code{WX|YZ}, \code{WY|XZ},
#'   \code{WZ|XY}, \code{XY|WZ}, \code{XZ|WY}, \code{YZ|WX}, \code{ZX|WY}, \code{XW|YZ},
#'   \code{YW|XZ}, \code{YX|WZ}, \code{ZW|XY}, \code{ZX|WY}, \code{ZY|WX} when
#'   \code{x} is generated by \code{historic_sim()}.  Default is
#'   \code{WX|YZ}.  When \code{x} is generated by \code{simple_sim()},
#'   \code{tab_type} is ignored.
#' @param smooth A true/false parameter indicating whether smoothed results
#'   should be plotted. Note, smoothing of simulation results requires the length of
#'   \code{subj_per_arm_val} or \code{a0_val} or \code{effect_val} or
#'   \code{rand_control_diff_val}, whichever populates the x-axis on the graph to
#'   contain enough elements to justify the smoothing.  No checking occurs to
#'   determine if enough elements are present to justify smoothing.  Default is
#'   \code{FALSE}.
#' @param plot_out A true/false parameter indicating whether the plot should be
#'   produced.  This parameter is useful if the user only wants a table of smoothed
#'   values.  Default is \code{TRUE}.
#' @param subj_per_arm_val Must be non-missing, if \code{x} is generated
#'   by \code{historic_sim()} and sample size is being held constant.
#'   If \code{x} is generated by \code{historic_sim()} and sample size
#'   is being held constant, \code{subj_per_arm_val} must equal a value submitted
#'   to \code{historic_sim()} within the \code{subj_per_arm} parameter.  When
#'   \code{x} is generated by \code{simple_sim()}, \code{subj_per_arm_val}
#'   is ignored.
#' @param a0_val Must be non-missing, if \code{x} is generated
#'   by \code{historic_sim()} and a0, the power prior parameter, is being held
#'   constant.  If \code{x} is generated by \code{historic_sim()} and
#'   a0 is being held constant, \code{a0_val} must equal a value submitted
#'   to \code{historic_sim()} within the \code{a0_val} parameter.  When
#'   \code{x} is generated by \code{simple_sim()}, \code{a0_val} is
#'   ignored.
#' @param effect_val  Must be non-missing, if \code{x} is generated
#'   by \code{historic_sim()} and effect is being held constant.  If
#'   \code{x} is generated by \code{historic_sim()} and effect is being
#'   held constant, \code{effect_val} must equal a value submitted to
#'   \code{historic_sim()} within the \code{effect_vals} parameter.  When
#'   \code{x} is generated by \code{simple_sim()}, \code{effect_val} is
#'   ignored.
#' @param rand_control_diff_val Must be non-missing, if \code{x} is
#'   generated by \code{historic_sim()} and differences between randomized
#'   and historical controls are being held constant.  If \code{x}
#'   is generated by \code{historic_sim()} and control differences are being
#'   held constant, \code{rand_control_diff_val} must equal a value submitted to
#'   \code{historic_sim()} within the \code{rand_control_diff} parameter.  When
#'   \code{x} is generated by \code{simple_sim()},
#'   \code{rand_control_diff_val} is ignored.
#' @param span The \code{span} parameter value for a call \code{loess()}.  Default is 0.75.  If
#'   \code{span} is a single number, then that value will be used to smooth the data in all
#'   columns of the table being plotted.  If \code{span} is a vector, then it must have length
#'   equal to the number of columns being plotted.
#' @param degree The \code{degree} parameter value for a call \code{loess()}.  Default is 2.
#'   The value of \code{degree} will be used for all columns being plotted.
#' @param family The \code{family} parameter value for a call \code{loess()}.  Default is
#'   "\code{gaussian}".  The value of \code{family} will be used for all columns being plotted.
#' @param title Title for the plot.
#' @param ylim Lower and upper limits for y-axis of plot.
#' @param ...	further arguments passed to or from other methods.
#'
#' @return \code{plot()} returns a plot for a two dimensional array of simulation
#'   results.  If \code{smooth} is \code{TRUE}, then the plot is based on a smoothed
#'   version of the simulation results. If \code{smooth} is \code{FALSE}, then the plot
#'   is based on the raw data from the simulation results.  What actually is plotted
#'   depends on the value of \code{measure}.  If \code{plot_out} is \code{FALSE}, the
#'   plot is not created.  This option is useful when the user wants a table of smoothed
#'   simulation results but does not want the plot. Smoothing of simulation results
#'   requires the length of \code{subj_per_arm_val} or \code{a0_val} or \code{effect_val}
#'   or \code{rand_control_diff_val}, whichever populates the x-axis on the graph to
#'   contain enough elements to justify the smoothing.  No checking occurs to
#'   determine if enough elements are present to justify smoothing.
#'
#' @examples
#' #Run a Weibull simulation, using simple_sim().
#' #For meaningful results, trial_reps needs to be much larger than 2.
#' weibull_test <- simple_sim(trial_reps = 2, outcome_type = "weibull",
#'                            subj_per_arm = c(50, 100, 150, 200),
#'                            effect_vals = c(0.6, 1),
#'                            control_parms = c(2.82487,3), time_vec = NULL,
#'                            censor_value = NULL, alpha = 0.05,
#'                            get_var = TRUE, get_bias = TRUE, get_mse = TRUE,
#'                            seedval=123, quietly=TRUE)
#'
#' #Create a plot of the power simulation results.
#' plot(x=weibull_test, measure="power", tab_type=NULL,
#'      smooth=FALSE, plot_out=TRUE, subj_per_arm_val=NULL, a0_val=NULL,
#'      effect_val=NULL, rand_control_diff_val=NULL)
#' #Create a plot of the hazard ratio simulation results.
#' plot(x=weibull_test, measure="est", tab_type=NULL,
#'      smooth=FALSE, plot_out=TRUE, subj_per_arm_val=NULL, a0_val=NULL,
#'      effect_val=NULL, rand_control_diff_val=NULL)
#' #Create a plot of the hazard ratio variance simulation results.
#' plot(x=weibull_test, measure="var", tab_type=NULL,
#'      smooth=FALSE, plot_out=TRUE, subj_per_arm_val=NULL, a0_val=NULL,
#'      effect_val=NULL, rand_control_diff_val=NULL)
#' #Create a plot of the hazard ratio bias simulation results.
#' plot(x=weibull_test, measure="bias", tab_type=NULL,
#'      smooth=FALSE, plot_out=TRUE, subj_per_arm_val=NULL, a0_val=NULL,
#'      effect_val=NULL, rand_control_diff_val=NULL)
#' #Create a plot of the hazard ratio mse simulation results.
#' plot(x=weibull_test, measure="mse", tab_type=NULL,
#'      smooth=FALSE, plot_out=TRUE, subj_per_arm_val=NULL, a0_val=NULL,
#'      effect_val=NULL, rand_control_diff_val=NULL)
#'
#' \donttest{
#' #Run a second Weibull simulation, using simple_sim() and smooth the plot.
#' #For meaningful results, trial_reps needs to be larger than 100.
#' weibull_test2 <- simple_sim(trial_reps = 100, outcome_type = "weibull",
#'                             subj_per_arm = c(50, 75, 100, 125, 150, 175, 200, 225, 250),
#'                             effect_vals = c(0.6, 1, 1.4),
#'                             control_parms = c(2.82487,3), time_vec = NULL,
#'                             censor_value = NULL, alpha = 0.05, get_var = TRUE,
#'                             get_bias = TRUE, get_mse = TRUE, seedval=123,
#'                             quietly=TRUE)
#'
#' #Create a plot of the power simulation results.
#' plot(x=weibull_test2, measure="power", tab_type=NULL,
#'      smooth=TRUE, plot_out=TRUE, subj_per_arm_val=NULL, a0_val=NULL,
#'      effect_val=NULL, rand_control_diff_val=NULL, span=c(1,1,1))
#' }
#'
#' \donttest{
#' #Run a third weibull simulation, using historic_sim().
#' #Note: historic_sim() can take a while to run.
#' #Generate a sample of historical data for use in example.
#' set.seed(2250)
#' SampleHistData <- genweibulldata(sample_size=60, scale1=2.82487,
#'                                  hazard_ratio=0.6, common_shape=3,
#'                                  censor_value=3)
#' histdata <- subset(SampleHistData, subset=(treatment==0))
#' histdata$id <- histdata$id+10000
#'
#' #For meaningful results, trial_reps needs to be larger than 100.
#' weibull_test3 <- historic_sim(trial_reps = 100, outcome_type = "weibull",
#'                               subj_per_arm = c(50, 100, 150, 200, 250),
#'                               a0_vals = c(0, 0.33, 0.67, 1),
#'                               effect_vals = c(0.6, 1, 1.4),
#'                               rand_control_diff = c(0.8, 1, 1.2),
#'                               hist_control_data = histdata, time_vec = NULL,
#'                               censor_value = 3, alpha = 0.05, get_var = TRUE,
#'                               get_bias = TRUE, get_mse = TRUE, seedval=123,
#'                               quietly=TRUE)
#'
#' #Create a plot of the power simulation results.
#' plot(x=weibull_test3, measure="power", tab_type="WX|YZ",
#'      smooth=FALSE, plot_out=TRUE, effect_val=0.6,
#'      rand_control_diff_val=1.0)
#' }
#'
#' \donttest{
#' #Run a Gaussian simulation, using historic_sim()
#' #Generate a sample of historical Gaussian data for use in example.
#' set.seed(2250)
#' samplehistdata <- gengaussiandata(sample_size=60, mu1=25, mean_diff=0, common_sd=3)
#' histdata <- subset(samplehistdata, subset=(treatment==0))
#' histdata$id <- histdata$id+10000
#'
#' #For meaningful results, trial_reps needs to be larger than 100.
#' gaussian_test <- historic_sim(trial_reps = 100, outcome_type = "gaussian",
#'                              subj_per_arm = c(150),
#'                              a0_vals = c(1.0),
#'                              effect_vals = c(0.15),
#'                              rand_control_diff = c(-4.0,-3.5,-3.0,-2.5,-2.0,
#'                                                    -1.5,-1.0,-0.5,0,0.5,1.0),
#'                              hist_control_data = histdata, time_vec = NULL,
#'                              censor_value = 3, alpha = 0.05, get_var = TRUE,
#'                              get_bias = TRUE, get_mse = TRUE, seedval=123,
#'                              quietly=TRUE)
#' test_table <- print(x=gaussian_test, measure="power",
#'                          tab_type=NULL, effect_val=NULL,
#'                          subj_per_arm_val=NULL)
#' print(test_table)
#' #Create a plot of the power simulation results.
#' plot(x=gaussian_test, measure="power", tab_type=NULL,
#'      smooth=TRUE, plot_out=TRUE, effect_val=NULL,
#'      rand_control_diff_val=NULL)
#' }
#'
#' \donttest{
#' #Generate a sample of historical pwe data for use in example.
#' set.seed(2250)
#' nvalHC <- 60
#' time.vec <- c(0.3,0.9,1.5,2.1,2.4)
#' lambdaHC.vec <- c(0.19,0.35,0.56,0.47,0.38,0.34)
#' censor.value <- 3
#'
#' SampleHistData <- genpwedata(nvalHC, lambdaHC.vec, 1.0, time.vec, censor.value)
#' histdata <- subset(SampleHistData, subset=(treatment==0))
#' histdata$indicator <- 2 #If set to 2, then historical controls will be collapsed with
#' #randomized controls, when time_vec is re-considered and
#' #potentially restructured.  If set to 1, then historical
#' #controls will be treated as a separated cohort when
#' #time_vec is being assessed for restructuring.
#' histdata$id <- histdata$id+10000
#'
#' #Run a pwe simulation, using historic_sim().
#' #For meaningful results, trial_reps needs to be larger than 100.
#' pwe_test <- historic_sim(trial_reps = 100, outcome_type = "pwe",
#'                         subj_per_arm = c(25,50,75,100,125,150,175,200,225,250),
#'                         a0_vals = c(1.0),
#'                         effect_vals = c(0.6),
#'                         rand_control_diff = c(1.8),
#'                         hist_control_data = histdata, time_vec = time.vec,
#'                         censor_value = 3, alpha = 0.05, get_var = TRUE,
#'                         get_bias = TRUE, get_mse = TRUE, seedval=123,
#'                         quietly=TRUE)
#'
#' #Create a plot of the power simulation results.
#' plot(x=pwe_test, measure="power", tab_type=NULL,
#'      smooth=TRUE, plot_out=TRUE, effect_val=NULL,
#'      rand_control_diff_val=NULL)
#' }
#'
#' @export
plot.bayes_ctd_array <- function(x = NULL, measure = "power", tab_type = "WX|YZ", smooth = FALSE,
                                       plot_out = TRUE, subj_per_arm_val = NULL, a0_val = NULL, effect_val = NULL, rand_control_diff_val = NULL, span = 0.75,
                                       degree = 2, family = "gaussian", title=NULL, ylim=NULL, ...) {
  # (W): subj_per_arm,
  # (X): a0_val,
  # (Y): effect_val,
  # (Z): rand_control_diff
  # tab_type='WX|YZ', Line Plot of Design Characteristics: Sample Size by a0
  # tab_type='WY|XZ', Line Plot of Design Characteristics: Sample Size by Effect
  # tab_type='WZ|XY', Line Plot of Design Characteristics: Sample Size by Control Differences
  # tab_type='XY|WZ', Line Plot of Design Characteristics: a0 by Effect
  # tab_type='XZ|WY', Line Plot of Design Characteristics: a0 by Control Diffferences
  # tab_type='YZ|WX', Line Plot of Design Characteristics: Effect by Control Diffferences
  # tab_type='ZX|WY', Line Plot of Design Characteristics: Control Diffferences by a0
  # tab_type='XW|YZ', Line Plot of Design Characteristics: a0 by Sample Size
  # tab_type='YW|XZ', Line Plot of Design Characteristics: Effect by Sample Size
  # tab_type='YX|WZ', Line Plot of Design Characteristics: Effect by a0
  # tab_type='ZW|XY', Line Plot of Design Characteristics: Control Differences by Sample Size
  # tab_type='ZY|WX', Line Plot of Design Characteristics: Control Differences by Effect

  #Add simple code to help out the user with tab_type on matters where order is
  #not important and can be easily fixed.
  tab_type = toupper(tab_type)
  if (!((length(tab_type)== 0) && (typeof(tab_type) == "character"))){
    if (x$objtype == 'historic') {
      if (tab_type=='WX|ZY') tab_type <- 'WX|YZ'
      if (tab_type=='WY|ZX') tab_type <- 'WY|XZ'
      if (tab_type=='WZ|YX') tab_type <- 'WZ|XY'
      if (tab_type=='XY|ZW') tab_type <- 'XY|WZ'
      if (tab_type=='XZ|YW') tab_type <- 'XZ|WY'
      if (tab_type=='YZ|XW') tab_type <- 'YZ|WX'
      if (tab_type=='ZX|YW') tab_type <- 'ZX|WY'
      if (tab_type=='XW|ZY') tab_type <- 'XW|YZ'
      if (tab_type=='YW|ZX') tab_type <- 'YW|XZ'
      if (tab_type=='YX|ZW') tab_type <- 'YX|WZ'
      if (tab_type=='ZW|YX') tab_type <- 'ZW|XY'
      if (tab_type=='ZY|XW') tab_type <- 'ZY|WX'
    }
  }
  #Do the same for measure
  measure = tolower(measure)
  #The above code was added very late in initial production, so time was not
  # available to go through all code and remove code no longer necessary.

  #------------- Go through all the high level checks for proper input. -------------#
  if (x$objtype == 'simple') {
    tab_type <- "WY|XZ"
    a0_val <- 0
    rand_control_diff_val <- 1
    subj_per_arm_val <- x$subj_per_arm[1]
    effect_val <- x$effect_vals[1]
    print("Since simple_sim was used, tab_type was set to WY|XZ")
    print("Values for tab_type, subj_per_arm_val, a0_val, effect_val, and rand_control_diff_val were ignored")
    print("This works towards putting all results on a single plot, effect by sample size")
  }

  if (x$objtype == 'realsimple') {
    if (length(x$subj_per_arm) > 1) {
      print("Since only subj_per_arm vector has more than 1 element, tab_type was set to WX|YZ")
      print("This works towards putting all results in a single table")
      tab_type <- "WX|YZ"
      a0_val <- x$a0_vals[1]
      effect_val <- x$effect_vals[1]
      rand_control_diff_val <- x$rand_control_diff[1]
    }
    if (length(x$a0_vals) > 1) {
      print("Since only a0_vals vector has more than 1 element, tab_type was set to XY|WZ")
      print("This works towards putting all results in a single table")
      tab_type <- "XY|WZ"
      subj_per_arm_val <- x$subj_per_arm[1]
      effect_val <- x$effect_vals[1]
      rand_control_diff_val <- x$rand_control_diff[1]
    }
    if (length(x$effect_vals) > 1) {
      print("Since only effect_vals vector has more than 1 element, tab_type was set to YZ|WX")
      print("This works towards putting all results in a single table")
      tab_type <- "YZ|WX"
      subj_per_arm_val <- x$subj_per_arm[1]
      a0_val <- x$a0_vals[1]
      rand_control_diff_val <- x$rand_control_diff[1]
    }
    if (length(x$rand_control_diff) > 1) {
      print("Since only rand_control_diff vector has more than 1 element, tab_type was set to ZX|WY")
      print("This works towards putting all results in a single table")
      tab_type <- "ZX|WY"
      subj_per_arm_val <- x$subj_per_arm[1]
      a0_val <- x$a0_vals[1]
      effect_val <- x$effect_vals[1]
    }
  }

  if (x$objtype == 'historic') {
    # Check to see if two dimensions of x have only 1 element each. Is so set tab_type to an appropriate
    # value.
    if (length(x$effect_vals) == 1 & length(x$rand_control_diff) == 1) {
      chkflag1 <- 0
      if(!((length(tab_type)== 0) && (typeof(tab_type) == "character"))){
        if ((tab_type != "WX|YZ") & (tab_type != "XW|YZ")) {
          tab_type <- "WX|YZ"
          chkflag1 <- 1
        }
      }
      if( (length(tab_type)== 0) && (typeof(tab_type) == "character")){
        tab_type <- "WX|YZ"
        chkflag1 <- 1
      }
      effect_val <- x$effect_vals[1]
      rand_control_diff_val <- x$rand_control_diff[1]
      if (chkflag1 == 1) {
        print("Since effect_vals and rand_control_diff vectors only have 1 element each and tab_type not equal to 'WX|YZ' or 'XW|YZ', tab_type was set to WX|YZ")
        print("This works towards putting all results on a single plot")
      }
    }
    if (length(x$a0_vals) == 1 & length(x$rand_control_diff) == 1) {
      chkflag1 <- 0
      if(!((length(tab_type)== 0) && (typeof(tab_type) == "character"))){
        if ((tab_type != "WY|XZ") & (tab_type != "YW|XZ")) {
          tab_type <- "WY|XZ"
          chkflag1 <- 1
        }
      }
      if( (length(tab_type)== 0) && (typeof(tab_type) == "character") ){
        tab_type <- "WY|XZ"
        chkflag1 <- 1
      }
      a0_val <- x$a0_vals[1]
      rand_control_diff_val <- x$rand_control_diff[1]
      if (chkflag1 == 1) {
        print("Since a0 and rand_control_diff vectors only have 1 element each and tab_type not equal to 'WY|XZ' or 'YW|XZ', tab_type was set to WY|XZ")
        print("This works towards putting all results on a single plot")
      }
    }
    if (length(x$subj_per_arm) == 1 & length(x$rand_control_diff) == 1) {
      chkflag1 <- 0
      if(!((length(tab_type)== 0) && (typeof(tab_type) == "character"))){
        if ((tab_type != "XY|WZ") & (tab_type != "YX|WZ")) {
          tab_type <- "XY|WZ"
          chkflag1 <- 1
        }
      }
      if( (length(tab_type)== 0) && (typeof(tab_type) == "character") ){
        tab_type <- "XY|WZ"
        chkflag1 <- 1
      }
      subj_per_arm_val <- x$subj_per_arm
      rand_control_diff_val <- x$rand_control_diff
      if (chkflag1 == 1) {
        print("Since sample size and rand_control_diff vectors only have 1 element each and tab_type not equal to 'XY|WZ' or 'YX|WZ', tab_type was set to XY|WZ")
        print("This works towards putting all results on a single plot")
      }
    }
    if (length(x$subj_per_arm) == 1 & length(x$effect_vals) == 1) {
      chkflag1 <- 0
      if(!((length(tab_type)== 0) && (typeof(tab_type) == "character"))){
        if ((tab_type != "XZ|WY") & (tab_type != "ZX|WY")) {
          tab_type <- "XZ|WY"
          chkflag1 <- 1
        }
      }
      if( (length(tab_type)== 0) && (typeof(tab_type) == "character") ){
        tab_type <- "XZ|WY"
        chkflag1 <- 1
      }
      subj_per_arm_val <- x$subj_per_arm
      effect_val <- x$effect_vals
      if (chkflag1 == 1) {
        print("Since sample size and effect_vals vectors only have 1 element each and tab_type not equal to 'XZ|WY' or 'ZX|WY', tab_type was set to XZ|WY")
        print("This works towards putting all results on a single plot")
      }
    }
    if (length(x$subj_per_arm) == 1 & length(x$a0_vals) == 1) {
      chkflag1 <- 0
      if(!((length(tab_type)== 0) && (typeof(tab_type) == "character"))){
        if ((tab_type != "YZ|WX") & (tab_type != "ZY|WX")) {
          tab_type <- "YZ|WX"
          chkflag1 <- 1
        }
      }
      if( (length(tab_type)== 0) && (typeof(tab_type) == "character") ){
        tab_type <- "YZ|WX"
        chkflag1 <- 1
      }
      subj_per_arm_val <- x$subj_per_arm
      a0_val <- x$a0_vals
      if (chkflag1 == 1) {
        print("Since sample size and a0 vectors only have 1 element each and tab_type not equal to 'YZ|WX' or 'ZY|WX', tab_type was set to YZ|WX")
        print("This works towards putting all results on a single plot")
      }
    }
    if (length(x$a0_vals) == 1 & length(x$effect_vals) == 1) {
      chkflag1 <- 0
      if(!((length(tab_type)== 0) && (typeof(tab_type) == "character"))){
        if ((tab_type != "WZ|XY") & (tab_type != "ZW|XY")) {
          tab_type <- "WZ|XY"
          chkflag1 <- 1
        }
      }
      if( (length(tab_type)== 0) && (typeof(tab_type) == "character") ){
        tab_type <- "WZ|XY"
        chkflag1 <- 1
      }
      a0_val <- x$a0_vals
      effect_val <- x$effect_vals
      if (chkflag1 == 1) {
        print("Since a0 and effect_vals vectors only have 1 element each and tab_type not equal to 'WZ|XY' or 'ZW|XY', tab_type was set to WZ|XY")
        print("This works towards putting all results on a single plot")
      }
    }
  }

  plot_error_checks(x, measure, tab_type, smooth, plot_out, subj_per_arm_val, a0_val, effect_val,
                          rand_control_diff_val, span)

  #If title is NULL, create a title from the two parameters that are fixed.
  if (is.null(title)){
    # (W): subj_per_arm_val,
    # (X): a0_val,
    # (Y): effect_val,
    # (Z): rand_control_diff_val
    #WX|YZ, WY|XZ, WZ|XY, XY|WZ, XZ|WY, YZ|WX, ZX|WY, XW|YZ, YW|XZ, YX|WZ, ZW|XY, ZX|WY, or ZY|WX

    if (tab_type == 'WX|YZ'){
      title = paste0('effect_val =', effect_val,' and rand_control_diff_val =', rand_control_diff_val)
    }
    if (tab_type == 'WY|XZ'){
      title = paste0('a0_val =', a0_val,' and rand_control_diff_val =', rand_control_diff_val)
    }
    if (tab_type == 'WZ|XY'){
      title = paste0('a0_val =', a0_val,' and effect_val =', effect_val)
    }
    if (tab_type == 'XY|WZ'){
      title = paste0('subj_per_arm_val =', subj_per_arm_val,' and rand_control_diff_val =', rand_control_diff_val)
    }
    if (tab_type == 'XZ|WY'){
      title = paste0('subj_per_arm_val =', subj_per_arm_val,' and effect_val =', effect_val)
    }
    if (tab_type == 'YZ|WX'){
      title = paste0('subj_per_arm_val =', subj_per_arm_val,' and a0_val =', a0_val)
    }
    if (tab_type == 'ZX|WY'){
      title = paste0('subj_per_arm_val =', subj_per_arm_val,' and effect_val =', effect_val)
    }
    if (tab_type == 'XW|YZ'){
      title = paste0('effect_val =', effect_val,' and rand_control_diff_val =', rand_control_diff_val)
    }
    if (tab_type == 'YW|XZ'){
      title = paste0('a0_val =', a0_val,' and rand_control_diff_val =', rand_control_diff_val)
    }
    if (tab_type == 'YX|WZ'){
      title = paste0('subj_per_arm_val =', subj_per_arm_val,' and rand_control_diff_val =', rand_control_diff_val)
    }
    if (tab_type == 'ZW|XY'){
      title = paste0('a0_val =', a0_val,' and effect_val =', effect_val)
    }
    if (tab_type == 'ZX|WY'){
      title = paste0('subj_per_arm_val =', subj_per_arm_val,' and effect_val =', effect_val)
    }
    if (tab_type == 'ZY|WX'){
      title = paste0('subj_per_arm_val =', subj_per_arm_val,' and a0_val =', a0_val)
    }
  }

  if (tab_type == "WX|YZ") {
    # Line Plot of Design Characteristics: Sample Size by a0
    # First generate the table using a call to print.
    table_mat <- print(x = x, measure = measure, tab_type = tab_type, effect_val = effect_val,
                             rand_control_diff_val = rand_control_diff_val, print_chg_warn = 0)
    table_mat <- data.frame(table_mat)
    if (dim(table_mat)[2]==1){
      tempcolname <- paste('X', as.character(x$a0_vals), sep='')
      colnames(table_mat) <- tempcolname
    }
    # Use a call to stats::loess to smooth simulation results is requested.
    if (smooth == TRUE) {
      colcnt <- dim(table_mat)[2]
      if (length(span) == 1)  #If span is single number, then use it to smooth all columns.
        span <- rep(span, colcnt)
      y <- as.numeric(rownames(table_mat))
      for (i in 1:colcnt) {
        loessfit <- stats::loess(table_mat[, i] ~ y, span = span[i], degree = degree, family = family)
        table_mat[, i] <- stats::predict(loessfit, data.frame(y = y))
      }
    }
    # If no plot needed, then return the table.  A user should only use this
    # option when they want to use the plot method to get a table of
    # smoothed results.
    if (plot_out == FALSE) {
      return(table_mat)
    }
    if (plot_out == TRUE) {
      #Shape the data for plotting.
      table_mat$row_names <- as.numeric(rownames(table_mat))
      table_dat <- reshape2::melt(table_mat, id = "row_names")
      colnames(table_dat) <- c("sample_size", "a0_val", measure)
      table_dat$sample_size <- as.numeric(table_dat$sample_size)

      #Construct proper values for the line variable.
      avec <- NULL
      if (table_dat$a0_val[1] == "table_mat") {
        table_dat$a0_val <- "0"
      } else if (table_dat$a0_val[1] != "table_mat") {
        for (i in 1:dim(table_dat)[1]) {
          avec[i] <- as.numeric(unlist(strsplit(as.character(table_dat[i, 2]), split = "X", fixed = TRUE)))[2]
        }
        table_dat$a0_val <- avec
        table_dat$a0_val <- as.character(table_dat$a0_val)
      }

      #Create the appropriate graph depending on the value of measure.
      if (measure == "power") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'power', group = 'a0_val', color = 'a0_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Sample Size', y='Power') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "est") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'est', group = 'a0_val', color = 'a0_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Sample Size', y='Estimate') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "var") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'var', group = 'a0_val', color = 'a0_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Sample Size', y='Variance') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "bias") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'bias', group = 'a0_val', color = 'a0_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Sample Size', y='Bias')  +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "mse") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'mse', group = 'a0_val', color = 'a0_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Sample Size', y='MSE')  +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
    }
  }
  if (tab_type == "WY|XZ") {
    # Line Plot of Design Characteristics: Sample Size by Effect
    # First generate the table using a call to print.
    table_mat <- print(x = x, measure = measure, tab_type = "WY|XZ", a0_val = a0_val,
                             rand_control_diff_val = rand_control_diff_val, print_chg_warn = 0)
    table_mat <- data.frame(table_mat)
    if (dim(table_mat)[2]==1){
      tempcolname <- paste('X', as.character(x$effect_vals), sep='')
      colnames(table_mat) <- tempcolname
    }
    # Use a call to stats::loess to smooth simulation results is requested.
    if (smooth == TRUE) {
      colcnt <- dim(table_mat)[2]
      if (length(span) == 1) #If span is single number, then use it to smooth all columns.
        span <- rep(span, colcnt)
      y <- as.numeric(rownames(table_mat))
      for (i in 1:colcnt) {
        loessfit <- stats::loess(table_mat[, i] ~ y, span = span[i], degree = degree, family = family)
        table_mat[, i] <- stats::predict(loessfit, data.frame(y = y))
      }
    }
    # If no plot needed, then return the table.  A user should only use this
    # option when they want to use the plot method to get a table of
    # smoothed results.
    if (plot_out == FALSE) {
      return(table_mat)
    }
    if (plot_out == TRUE) {
      #Shape the data for plotting.
      table_mat$row_names <- as.numeric(rownames(table_mat))
      table_dat <- reshape2::melt(table_mat, id = "row_names")
      colnames(table_dat) <- c("sample_size", "effect_val", measure)
      table_dat$sample_size <- as.numeric(table_dat$sample_size)

      #Construct proper values for the line variable.
      avec <- NULL
      for (i in 1:dim(table_dat)[1]) {
        avec[i] <- as.numeric(unlist(strsplit(as.character(table_dat[i, 2]), split = "X", fixed = TRUE)))[2]
      }
      table_dat$effect_val <- avec
      table_dat$effect_val <- as.character(table_dat$effect_val)

      #Create the appropriate graph depending on the value of measure.
      if (measure == "power") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'power', group = 'effect_val', color = 'effect_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='Sample Size', y='Power') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "est") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'est', group = 'effect_val', color = 'effect_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='Sample Size', y='Estimate') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "var") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'var', group = 'effect_val', color = 'effect_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='Sample Size', y='Variance') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "bias") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'bias', group = 'effect_val', color = 'effect_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='Sample Size', y='Bias') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "mse") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'mse', group = 'effect_val', color = 'effect_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='Sample Size', y='MSE') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
    }
  }
  if (tab_type == "WZ|XY") {
    # Line Plot of Design Characteristics: Sample Size by Control Differences
    # First generate the table using a call to print.
    table_mat <- print(x = x, measure = measure, tab_type = "WZ|XY", a0_val = a0_val,
                             effect_val = effect_val, print_chg_warn = 0)
    table_mat <- data.frame(table_mat)
    if (dim(table_mat)[2]==1){
      tempcolname <- paste('X', as.character(x$rand_control_diff), sep='')
      colnames(table_mat) <- tempcolname
    }
    # Use a call to stats::loess to smooth simulation results is requested.
    if (smooth == TRUE) {
      colcnt <- dim(table_mat)[2]
      if (length(span) == 1) #If span is single number, then use it to smooth all columns.
        span <- rep(span, colcnt)
      y <- as.numeric(rownames(table_mat))
      for (i in 1:colcnt) {
        loessfit <- stats::loess(table_mat[, i] ~ y, span = span[i], degree = degree, family = family)
        table_mat[, i] <- stats::predict(loessfit, data.frame(y = y))
      }
    }
    # If no plot needed, then return the table.  A user should only use this
    # option when they want to use the plot method to get a table of
    # smoothed results.
    if (plot_out == FALSE) {
      return(table_mat)
    }
    if (plot_out == TRUE) {
      #Shape the data for plotting.
      table_mat$row_names <- as.numeric(rownames(table_mat))
      table_dat <- reshape2::melt(table_mat, id = "row_names")
      colnames(table_dat) <- c("sample_size", "control_differences", measure)
      table_dat$sample_size <- as.numeric(table_dat$sample_size)

      #Construct proper values for the line variable.
      avec <- NULL
      if (table_dat$control_differences[1] == "table_mat") {
        table_dat$control_differences <- "1"
      } else if (table_dat$control_differences[1] != "table_mat") {
        for (i in 1:dim(table_dat)[1]) {
          avec[i] <- as.numeric(unlist(strsplit(as.character(table_dat[i, 2]), split = "X", fixed = TRUE)))[2]
        }
        table_dat$control_differences <- avec
        table_dat$control_differences <- as.character(table_dat$control_differences)
      }

      #Create the appropriate graph depending on the value of measure.
      if (measure == "power") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'power', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='Sample Size', y='Power') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "est") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'est', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='Sample Size', y='Estimate') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "var") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'var', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='Sample Size', y='Variance')+
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "bias") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'bias', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='Sample Size', y='Bias') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "mse") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'sample_size', y = 'mse', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='Sample Size', y='MSE') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }

    }
  }

  if (tab_type == "XY|WZ") {
    # Line Plot of Design Characteristics: a0 by Effect
    # First generate the table using a call to print.
    table_mat <- print(x = x, measure = measure, tab_type = "XY|WZ", subj_per_arm_val = subj_per_arm_val,
                             rand_control_diff_val = rand_control_diff_val, print_chg_warn = 0)
    table_mat <- data.frame(table_mat)
    if (dim(table_mat)[2]==1){
      tempcolname <- paste('X', as.character(x$effect_vals), sep='')
      colnames(table_mat) <- tempcolname
    }
    # Use a call to stats::loess to smooth simulation results is requested.
    if (smooth == TRUE) {
      colcnt <- dim(table_mat)[2]
      if (length(span) == 1) #If span is single number, then use it to smooth all columns.
        span <- rep(span, colcnt)
      y <- as.numeric(rownames(table_mat))
      for (i in 1:colcnt) {
        loessfit <- stats::loess(table_mat[, i] ~ y, span = span[i], degree = degree, family = family)
        table_mat[, i] <- stats::predict(loessfit, data.frame(y = y))
      }
    }
    # If no plot needed, then return the table.  A user should only use this
    # option when they want to use the plot method to get a table of
    # smoothed results.
    if (plot_out == FALSE) {
      return(table_mat)
    }
    if (plot_out == TRUE) {
      #Shape the data for plotting.
      table_mat$row_names <- as.numeric(rownames(table_mat))
      table_dat <- reshape2::melt(table_mat, id = "row_names")
      colnames(table_dat) <- c("a0_val", "effect_val", measure)
      table_dat$a0_val <- as.numeric(table_dat$a0_val)

      #Construct proper values for the line variable.
      avec <- NULL
      for (i in 1:dim(table_dat)[1]) {
        avec[i] <- as.numeric(unlist(strsplit(as.character(table_dat[i, 2]), split = "X", fixed = TRUE)))[2]
      }
      table_dat$effect_val <- avec
      table_dat$effect_val <- as.character(table_dat$effect_val)

      #Create the appropriate graph depending on the value of measure.
      if (measure == "power") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'power', group = 'effect_val', color = 'effect_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='a0 Value', y='Power') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "est") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'est', group = 'effect_val', color = 'effect_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='a0 Value', y='Estimate') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "var") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'var', group = 'effect_val', color = 'effect_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='a0 Value', y='Variance') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "bias") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'bias', group = 'effect_val', color = 'effect_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='a0 Value', y='Bias') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "mse") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'mse', group = 'effect_val', color = 'effect_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='a0 Value', y='MSE') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
    }
  }

  if (tab_type == "XZ|WY") {
    # Line Plot of Design Characteristics: a0 by Control Diffferences
    # First generate the table using a call to print.
    table_mat <- print(x = x, measure = measure, tab_type = "XZ|WY", subj_per_arm_val = subj_per_arm_val,
                             effect_val = effect_val, print_chg_warn = 0)
    table_mat <- data.frame(table_mat)
    if (dim(table_mat)[2]==1){
      tempcolname <- paste('X', as.character(x$rand_control_diff), sep='')
      colnames(table_mat) <- tempcolname
    }
    # Use a call to stats::loess to smooth simulation results is requested.
    if (smooth == TRUE) {
      colcnt <- dim(table_mat)[2]
      if (length(span) == 1) #If span is single number, then use it to smooth all columns.
        span <- rep(span, colcnt)
      y <- as.numeric(rownames(table_mat))
      for (i in 1:colcnt) {
        loessfit <- stats::loess(table_mat[, i] ~ y, span = span[i], degree = degree, family = family)
        table_mat[, i] <- stats::predict(loessfit, data.frame(y = y))
      }
    }
    # If no plot needed, then return the table.  A user should only use this
    # option when they want to use the plot method to get a table of
    # smoothed results.
    if (plot_out == FALSE) {
      return(table_mat)
    }
    if (plot_out == TRUE) {
      #Shape the data for plotting.
      table_mat$row_names <- as.numeric(rownames(table_mat))
      table_dat <- reshape2::melt(table_mat, id = "row_names")
      colnames(table_dat) <- c("a0_val", "control_differences", measure)
      table_dat$a0_val <- as.numeric(table_dat$a0_val)

      #Construct proper values for the line variable.
      avec <- NULL
      for (i in 1:dim(table_dat)[1]) {
        avec[i] <- as.numeric(unlist(strsplit(as.character(table_dat[i, 2]), split = "X", fixed = TRUE)))[2]
      }
      table_dat$control_differences <- avec
      table_dat$control_differences <- as.character(table_dat$control_differences)

      #Create the appropriate graph depending on the value of measure.
      if (measure == "power") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'power', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='a0 Value', y='Power') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "est") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'est', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='a0 Value', y='Estimate') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "var") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'var', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='a0 Value', y='Variance') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "bias") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'bias', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='a0 Value', y='Bias') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "mse") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'mse', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='a0 Value', y='MSE') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }

    }
  }

  if (tab_type == "YZ|WX") {
    # Line Plot of Design Characteristics: Effect by Control Diffferences
    # First generate the table using a call to print.
    table_mat <- print(x = x, measure = measure, tab_type = "YZ|WX", subj_per_arm_val = subj_per_arm_val,
                             a0_val = a0_val, print_chg_warn = 0)
    table_mat <- data.frame(table_mat)
    if (dim(table_mat)[2]==1){
      tempcolname <- paste('X', as.character(x$rand_control_diff), sep='')
      colnames(table_mat) <- tempcolname
    }
    # Use a call to stats::loess to smooth simulation results is requested.
    if (smooth == TRUE) {
      colcnt <- dim(table_mat)[2]
      if (length(span) == 1) #If span is single number, then use it to smooth all columns.
        span <- rep(span, colcnt)
      y <- as.numeric(rownames(table_mat))
      for (i in 1:colcnt) {
        loessfit <- stats::loess(table_mat[, i] ~ y, span = span[i], degree = degree, family = family)
        table_mat[, i] <- stats::predict(loessfit, data.frame(y = y))
      }
    }
    # If no plot needed, then return the table.  A user should only use this
    # option when they want to use the plot method to get a table of
    # smoothed results.
    if (plot_out == FALSE) {
      return(table_mat)
    }
    if (plot_out == TRUE) {
      #Shape the data for plotting.
      table_mat$row_names <- as.numeric(rownames(table_mat))
      table_dat <- reshape2::melt(table_mat, id = "row_names")
      colnames(table_dat) <- c("effect_val", "control_differences", measure)
      table_dat$effect_val <- as.numeric(table_dat$effect_val)

      #Construct proper values for the line variable.
      avec <- NULL
      if (table_dat$control_differences[1] == "table_mat") {
        table_dat$control_differences <- "1"
      } else if (table_dat$control_differences[1] != "table_mat") {
        for (i in 1:dim(table_dat)[1]) {
          avec[i] <- as.numeric(unlist(strsplit(as.character(table_dat[i, 2]), split = "X", fixed = TRUE)))[2]
        }
        table_dat$control_differences <- avec
        table_dat$control_differences <- as.character(table_dat$control_differences)
      }

      #Create the appropriate graph depending on the value of measure.
      if (measure == "power") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'power', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='Effect Value', y='Power') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "est") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'est', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='Effect Value', y='Estimate') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "var") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'var', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='Effect Value', y='Variance') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "bias") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'bias', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='Effect Value', y='Bias') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "mse") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'mse', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='Effect Value', y='MSE') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
    }
  }

  if (tab_type == "ZX|WY") {
    # Line Plot of Design Characteristics: Control Diffferences by a0
    # First generate the table using a call to print.
    table_mat <- print(x = x, measure = measure, tab_type = "ZX|WY", subj_per_arm_val = subj_per_arm_val,
                             effect_val = effect_val, print_chg_warn = 0)
    #The following if block makes the 1D case work, if it occurs.
    #The return of a 1D vector from print(tab_type="ZX|WY")
    #  which is a 1byK row vector with class=matrix causes problems
    #when I make it into a dataframe.  The following code re-classifies
    #the object as a vector if necessary.
    if (is.matrix(table_mat) & dim(table_mat)[1]==1){
      colnamevec <- colnames(table_mat)
      table_mat <- as.vector(table_mat)
      names(table_mat) <- colnamevec
    }
    table_mat <- data.frame(table_mat)
    if (dim(table_mat)[2]==1){
      tempcolname <- paste('X', as.character(x$a0_vals), sep='')
      colnames(table_mat) <- tempcolname
    }
    # Use a call to stats::loess to smooth simulation results is requested.
    if (smooth == TRUE) {
      colcnt <- dim(table_mat)[2]
      if (length(span) == 1) #If span is single number, then use it to smooth all columns.
        span <- rep(span, colcnt)
      y <- as.numeric(rownames(table_mat))
      for (i in 1:colcnt) {
        loessfit <- stats::loess(table_mat[, i] ~ y, span = span[i], degree = degree, family = family)
        table_mat[, i] <- stats::predict(loessfit, data.frame(y = y))
      }
    }
    # If no plot needed, then return the table.  A user should only use this
    # option when they want to use the plot method to get a table of
    # smoothed results.
    if (plot_out == FALSE) {
      return(table_mat)
    }
    if (plot_out == TRUE) {
      #Shape the data for plotting.
      table_mat$row_names <- as.numeric(rownames(table_mat))
      table_dat <- reshape2::melt(table_mat, id = "row_names")
      colnames(table_dat) <- c("control_differences", "a0_val", measure)
      table_dat$control_differences <- as.numeric(table_dat$control_differences)

      #Construct proper values for the line variable.
      avec <- NULL
      if (table_dat$a0_val[1] == "table_mat") {
        table_dat$a0_val <- "1"
      } else if (table_dat$a0_val[1] != "table_mat") {
        for (i in 1:dim(table_dat)[1]) {
          avec[i] <- as.numeric(unlist(strsplit(as.character(table_dat[i, 2]), split = "X", fixed = TRUE)))[2]
        }
        table_dat$a0_val <- avec
        table_dat$a0_val <- as.character(table_dat$a0_val)
      }

      #Create the appropriate graph depending on the value of measure.
      if (measure == "power") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'power', group = 'a0_val', color = 'a0_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Control Differences', y='Power') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "est") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'est', group = 'a0_val', color = 'a0_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Control Differences', y='Estimate') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "var") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'var', group = 'a0_val', color = 'a0_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Control Differences', y='Variance') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "bias") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'bias', group = 'a0_val', color = 'a0_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Control Differences', y='Bias') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "mse") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'mse', group = 'a0_val', color = 'a0_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Control Differences', y='MSE') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
    }
  }

  if (tab_type == "XW|YZ") {
    # Line Plot of Design Characteristics: a0 by sample size
    # First generate the table using a call to print.
    table_mat <- print(x = x, measure = measure, tab_type = "XW|YZ", effect_val = effect_val,
                             rand_control_diff_val = rand_control_diff_val, print_chg_warn = 0)
    #The following if block makes the 1D case work, if it occurs.
    #The return of a 1D vector from print(tab_type="ZX|WY")
    #  which is a 1byK row vector with class=matrix causes problems
    #when I make it into a dataframe.  The following code re-classifies
    #the object as a vector if necessary.
    if (is.matrix(table_mat) & dim(table_mat)[1]==1){
      colnamevec <- colnames(table_mat)
      table_mat <- as.vector(table_mat)
      names(table_mat) <- colnamevec
    }
    table_mat <- data.frame(table_mat)
    if (dim(table_mat)[2]==1){
      tempcolname <- paste('X', as.character(x$subj_per_arm), sep='')
      colnames(table_mat) <- tempcolname
    }
    # Use a call to stats::loess to smooth simulation results is requested.
    if (smooth == TRUE) {
      colcnt <- dim(table_mat)[2]
      if (length(span) == 1) #If span is single number, then use it to smooth all columns.
        span <- rep(span, colcnt)
      y <- as.numeric(rownames(table_mat))
      for (i in 1:colcnt) {
        loessfit <- stats::loess(table_mat[, i] ~ y, span = span[i], degree = degree, family = family)
        table_mat[, i] <- stats::predict(loessfit, data.frame(y = y))
      }
    }
    # If no plot needed, then return the table.  A user should only use this
    # option when they want to use the plot method to get a table of
    # smoothed results.
    if (plot_out == FALSE) {
      return(table_mat)
    }
    if (plot_out == TRUE) {
      #Shape the data for plotting.
      table_mat$row_names <- as.numeric(rownames(table_mat))
      table_dat <- reshape2::melt(table_mat, id = "row_names")
      colnames(table_dat) <- c("a0_val", "subj_per_arm", measure)
      table_dat$a0_val <- as.numeric(table_dat$a0_val)

      #Construct proper values for the line variable.
      avec <- NULL
      if (table_dat$subj_per_arm[1] == "table_mat") {
        table_dat$subj_per_arm <- "1"
      } else if (table_dat$subj_per_arm[1] != "table_mat") {
        for (i in 1:dim(table_dat)[1]) {
          avec[i] <- as.numeric(unlist(strsplit(as.character(table_dat[i, 2]), split = "X", fixed = TRUE)))[2]
        }
        table_dat$subj_per_arm <- avec
        table_dat$subj_per_arm <- as.character(table_dat$subj_per_arm)
      }

      #Create the appropriate graph depending on the value of measure.
      if (measure == "power") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'power', group = 'subj_per_arm', color = 'subj_per_arm')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Subjects per Arm', x='a0_val', y='Power') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "est") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'est', group = 'subj_per_arm', color = 'subj_per_arm')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Subjects per Arm', x='a0_val', y='Estimate') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "var") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'var', group = 'subj_per_arm', color = 'subj_per_arm')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Subjects per Arm', x='a0_val', y='Variance') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "bias") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'bias', group = 'subj_per_arm', color = 'subj_per_arm')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Subjects per Arm', x='a0_val', y='Bias') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "mse") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'a0_val', y = 'mse', group = 'subj_per_arm', color = 'subj_per_arm')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Subjects per Arm', x='a0_val', y='MSE') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
    }
  }

  if (tab_type == "YW|XZ") {
    # Line Plot of Design Characteristics: Effect by sample size
    # First generate the table using a call to print.
    table_mat <- print(x = x, measure = measure, tab_type = "YW|XZ",
                             a0_val = a0_val, rand_control_diff_val=rand_control_diff_val, print_chg_warn = 0)
    #The following if block makes the 1D case work, if it occurs.
    #The return of a 1D vector from print(tab_type="ZX|WY")
    #  which is a 1byK row vector with class=matrix causes problems
    #when I make it into a dataframe.  The following code re-classifies
    #the object as a vector if necessary.
    if (is.matrix(table_mat) & dim(table_mat)[1]==1){
      colnamevec <- colnames(table_mat)
      table_mat <- as.vector(table_mat)
      names(table_mat) <- colnamevec
    }
    table_mat <- data.frame(table_mat)
    if (dim(table_mat)[2]==1){
      tempcolname <- paste('X', as.character(x$subj_per_arm), sep='')
      colnames(table_mat) <- tempcolname
    }
    # Use a call to stats::loess to smooth simulation results is requested.
    if (smooth == TRUE) {
      colcnt <- dim(table_mat)[2]
      if (length(span) == 1) #If span is single number, then use it to smooth all columns.
        span <- rep(span, colcnt)
      y <- as.numeric(rownames(table_mat))
      for (i in 1:colcnt) {
        loessfit <- stats::loess(table_mat[, i] ~ y, span = span[i], degree = degree, family = family)
        table_mat[, i] <- stats::predict(loessfit, data.frame(y = y))
      }
    }
    # If no plot needed, then return the table.  A user should only use this
    # option when they want to use the plot method to get a table of
    # smoothed results.
    if (plot_out == FALSE) {
      return(table_mat)
    }
    if (plot_out == TRUE) {
      #Shape the data for plotting.
      table_mat$row_names <- as.numeric(rownames(table_mat))
      table_dat <- reshape2::melt(table_mat, id = "row_names")
      colnames(table_dat) <- c("effect_val", "subj_per_arm", measure)
      table_dat$effect_val <- as.numeric(table_dat$effect_val)

      #Construct proper values for the line variable.
      avec <- NULL
      if (table_dat$subj_per_arm[1] == "table_mat") {
        table_dat$subj_per_arm <- "1"
      } else if (table_dat$subj_per_arm[1] != "table_mat") {
        for (i in 1:dim(table_dat)[1]) {
          avec[i] <- as.numeric(unlist(strsplit(as.character(table_dat[i, 2]), split = "X", fixed = TRUE)))[2]
        }
        table_dat$subj_per_arm <- avec
        table_dat$subj_per_arm <- as.character(table_dat$subj_per_arm)
      }

      #Create the appropriate graph depending on the value of measure.
      if (measure == "power") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'power', group = 'subj_per_arm', color = 'subj_per_arm')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Subjects per Arm', x='Effect Value', y='Power') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "est") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'est', group = 'subj_per_arm', color = 'subj_per_arm')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Subjects per Arm', x='Effect Value', y='Estimate') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "var") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'var', group = 'subj_per_arm', color = 'subj_per_arm')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Subjects per Arm', x='Effect Value', y='Variance') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "bias") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'bias', group = 'subj_per_arm', color = 'subj_per_arm')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Subjects per Arm', x='Effect Value', y='Bias') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "mse") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'mse', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='Effect Value', y='MSE') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
    }
  }

  if (tab_type == "ZW|XY") {
    # Line Plot of Design Characteristics: Control differences by sample size
    # First generate the table using a call to print.
    table_mat <- print(x = x, measure = measure, tab_type = "ZW|XY",
                             a0_val = a0_val, effect_val=effect_val, print_chg_warn = 0)
    #The following if block makes the 1D case work, if it occurs.
    #The return of a 1D vector from print(tab_type="ZX|WY")
    #  which is a 1byK row vector with class=matrix causes problems
    #when I make it into a dataframe.  The following code re-classifies
    #the object as a vector if necessary.
    if (is.matrix(table_mat) & dim(table_mat)[1]==1){
      colnamevec <- colnames(table_mat)
      table_mat <- as.vector(table_mat)
      names(table_mat) <- colnamevec
    }
    table_mat <- data.frame(table_mat)
    if (dim(table_mat)[2]==1){
      tempcolname <- paste('X', as.character(x$subj_per_arm), sep='')
      colnames(table_mat) <- tempcolname
    }
    # Use a call to stats::loess to smooth simulation results is requested.
    if (smooth == TRUE) {
      colcnt <- dim(table_mat)[2]
      if (length(span) == 1) #If span is single number, then use it to smooth all columns.
        span <- rep(span, colcnt)
      y <- as.numeric(rownames(table_mat))
      for (i in 1:colcnt) {
        loessfit <- stats::loess(table_mat[, i] ~ y, span = span[i], degree = degree, family = family)
        table_mat[, i] <- stats::predict(loessfit, data.frame(y = y))
      }
    }
    # If no plot needed, then return the table.  A user should only use this
    # option when they want to use the plot method to get a table of
    # smoothed results.
    if (plot_out == FALSE) {
      return(table_mat)
    }
    if (plot_out == TRUE) {
      #Shape the data for plotting.
      table_mat$row_names <- as.numeric(rownames(table_mat))
      table_dat <- reshape2::melt(table_mat, id = "row_names")
      colnames(table_dat) <- c("control_differences", "subj_per_arm", measure)
      table_dat$control_differences <- as.numeric(table_dat$control_differences)

      #Construct proper values for the line variable.
      avec <- NULL
      if (table_dat$subj_per_arm[1] == "table_mat") {
        table_dat$subj_per_arm <- "1"
      } else if (table_dat$subj_per_arm[1] != "table_mat") {
        for (i in 1:dim(table_dat)[1]) {
          avec[i] <- as.numeric(unlist(strsplit(as.character(table_dat[i, 2]), split = "X", fixed = TRUE)))[2]
        }
        table_dat$subj_per_arm <- avec
        table_dat$subj_per_arm <- as.character(table_dat$subj_per_arm)
      }

      #Create the appropriate graph depending on the value of measure.
      if (measure == "power") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'power', group = 'subj_per_arm', color = 'subj_per_arm')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Subjects per Arm', x='Control Differences', y='Power') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "est") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'est', group = 'subj_per_arm', color = 'subj_per_arm')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Subjects per Arm', x='Control Differences', y='Estimate') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "var") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'var', group = 'subj_per_arm', color = 'subj_per_arm')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Subjects per Arm', x='Control Differences', y='Variance') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "bias") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'bias', group = 'subj_per_arm', color = 'subj_per_arm')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Subjects per Arm', x='Control Differences', y='Bias') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "mse") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'mse', group = 'control_differences', color = 'control_differences')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='Control\nDifferences', x='Control Differences', y='MSE') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
    }
  }

  if (tab_type == "ZY|WX") {
    # Line Plot of Design Characteristics: a0 by Effect
    # First generate the table using a call to print.
    table_mat <- print(x = x, measure = measure, tab_type = "ZY|WX", subj_per_arm_val = subj_per_arm_val,
                             a0_val = a0_val, print_chg_warn = 0)
    #The following if block makes the 1D case work, if it occurs.
    #The return of a 1D vector from print(tab_type="ZX|WY")
    #  which is a 1byK row vector with class=matrix causes problems
    #when I make it into a dataframe.  The following code re-classifies
    #the object as a vector if necessary.
    if (is.matrix(table_mat) & dim(table_mat)[1]==1){
      colnamevec <- colnames(table_mat)
      table_mat <- as.vector(table_mat)
      names(table_mat) <- colnamevec
    }
    table_mat <- data.frame(table_mat)
    if (dim(table_mat)[2]==1){
      tempcolname <- paste('X', as.character(x$effect_vals), sep='')
      colnames(table_mat) <- tempcolname
    }
    # Use a call to stats::loess to smooth simulation results is requested.
    if (smooth == TRUE) {
      colcnt <- dim(table_mat)[2]
      if (length(span) == 1) #If span is single number, then use it to smooth all columns.
        span <- rep(span, colcnt)
      y <- as.numeric(rownames(table_mat))
      for (i in 1:colcnt) {
        loessfit <- stats::loess(table_mat[, i] ~ y, span = span[i], degree = degree, family = family)
        table_mat[, i] <- stats::predict(loessfit, data.frame(y = y))
      }
    }
    # If no plot needed, then return the table.  A user should only use this
    # option when they want to use the plot method to get a table of
    # smoothed results.
    if (plot_out == FALSE) {
      return(table_mat)
    }
    if (plot_out == TRUE) {
      #Shape the data for plotting.
      table_mat$row_names <- as.numeric(rownames(table_mat))
      table_dat <- reshape2::melt(table_mat, id = "row_names")
      colnames(table_dat) <- c("control_differences", "effect_val", measure)
      table_dat$control_differences <- as.numeric(table_dat$control_differences)

      #Construct proper values for the line variable.
      avec <- NULL
      for (i in 1:dim(table_dat)[1]) {
        avec[i] <- as.numeric(unlist(strsplit(as.character(table_dat[i, 2]), split = "X", fixed = TRUE)))[2]
      }
      table_dat$effect_val <- avec
      table_dat$effect_val <- as.character(table_dat$effect_val)

      #Create the appropriate graph depending on the value of measure.
      if (measure == "power") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'power', group = 'effect_val', color = 'effect_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='Control Differences', y='Power') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "est") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'est', group = 'effect_val', color = 'effect_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='Control Differences', y='Estimate') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "var") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'var', group = 'effect_val', color = 'effect_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='Control Differences', y='Variance') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "bias") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'bias', group = 'effect_val', color = 'effect_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='Control Differences', y='Bias') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "mse") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'control_differences', y = 'mse', group = 'effect_val', color = 'effect_val')) + ggplot2::geom_line() +
          ggplot2::geom_point() + ggplot2::labs(color='Effect Value', x='Control Differences', y='MSE') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
    }
  }

  if (tab_type == "YX|WZ") {
    # Line Plot of Design Characteristics: Effect by a0
    # First generate the table using a call to print.
    table_mat <- print(x = x, measure = measure, tab_type = "YX|WZ", subj_per_arm_val = subj_per_arm_val,
                             rand_control_diff_val = rand_control_diff_val, print_chg_warn = 0)
    #The following if block makes the 1D case work, if it occurs.
    #The return of a 1D vector from print(tab_type="ZX|WY")
    #  which is a 1byK row vector with class=matrix causes problems
    #when I make it into a dataframe.  The following code re-classifies
    #the object as a vector if necessary.
    if (is.matrix(table_mat) & dim(table_mat)[1]==1){
      colnamevec <- colnames(table_mat)
      table_mat <- as.vector(table_mat)
      names(table_mat) <- colnamevec
    }
    table_mat <- data.frame(table_mat)
    if (dim(table_mat)[2]==1){
      tempcolname <- paste('X', as.character(x$a0_vals), sep='')
      colnames(table_mat) <- tempcolname
    }
    # Use a call to stats::loess to smooth simulation results is requested.
    if (smooth == TRUE) {
      colcnt <- dim(table_mat)[2]
      if (length(span) == 1) #If span is single number, then use it to smooth all columns.
        span <- rep(span, colcnt)
      y <- as.numeric(rownames(table_mat))
      for (i in 1:colcnt) {
        loessfit <- stats::loess(table_mat[, i] ~ y, span = span[i], degree = degree, family = family)
        table_mat[, i] <- stats::predict(loessfit, data.frame(y = y))
      }
    }
    # If no plot needed, then return the table.  A user should only use this
    # option when they want to use the plot method to get a table of
    # smoothed results.
    if (plot_out == FALSE) {
      return(table_mat)
    }
    if (plot_out == TRUE) {
      #Shape the data for plotting.
      table_mat$row_names <- as.numeric(rownames(table_mat))
      table_dat <- reshape2::melt(table_mat, id = "row_names")
      colnames(table_dat) <- c("effect_val", "a0_val", measure)
      table_dat$effect_val <- as.numeric(table_dat$effect_val)

      #Construct proper values for the line variable.
      avec <- NULL
      if (table_dat$a0_val[1] == "table_mat") {
        table_dat$a0_val <- "1"
      } else if (table_dat$a0_val[1] != "table_mat") {
        for (i in 1:dim(table_dat)[1]) {
          avec[i] <- as.numeric(unlist(strsplit(as.character(table_dat[i, 2]), split = "X", fixed = TRUE)))[2]
        }
        table_dat$a0_val <- avec
        table_dat$a0_val <- as.character(table_dat$a0_val)
      }

      #Create the appropriate graph depending on the value of measure.
      if (measure == "power") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'power', group = 'a0_val', color = 'a0_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Effect Value', y='Power') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "est") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'est', group = 'a0_val', color = 'a0_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Effect Value', y='Estimate') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "var") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'var', group = 'a0_val', color = 'a0_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Effect Value', y='Variance') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "bias") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'bias', group = 'a0_val', color = 'a0_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Effect Value', y='Bias') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
      if (measure == "mse") {
        p <- ggplot2::ggplot(data = table_dat, ggplot2::aes_string(x = 'effect_val', y = 'mse', group = 'a0_val', color = 'a0_val')) +
          ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::labs(color='a0 Value', x='Effect Value', y='MSE') +
          ggplot2::ggtitle(title)
        if (!is.null(ylim)){
          p <- p + ggplot2::ylim(ylim)
        }
        return(p)
      }
    }
  }
}
