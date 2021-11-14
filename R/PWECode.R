#' Generating function for Piece-wise Exponential Data.
#'
#' \code{genpwedata()} function used mainly internally by
#' \code{pwetrialsimulator()} function to generate data for a two-arm
#' clinical trial, experimental and control groups.  Can be used to generate
#' random trial data.
#'
#' @param sample_size  Number of subjects per arm.
#' @param lambda_vec  Set of lambdas passed to \code{eha::rpch()} through the
#'   levels parameter.  Used only in control arm.
#' @param hazard_ratio Desired Hazard Ratio between experimental and control groups.
#' @param time_vec Set of cutpoints passed to \code{eha::rpch()} through the
#'   cuts parameter.
#' @param censor_value Value at which time-to-event data are right censored.
#'
#' @return \code{genpwedata()} returns a data frame with columns: 'id', 'treatment',
#'   'event_time', 'status', and 'indicator'.
#'
#' @examples
#' nvalHC <- 60
#' time.vec <- c(0.3,0.9,1.5,2.1,2.4)
#' lambdaHC.vec <- c(0.19,0.35,0.56,0.47,0.38,0.34)
#' censor.value <- 3
#'
#' SampleHistData <- genpwedata(nvalHC, lambdaHC.vec, 1.0, time.vec, censor.value)
#' SampleHistData
#' @export
genpwedata <- function(sample_size, lambda_vec, hazard_ratio, time_vec, censor_value) {

    # --------------------------------------------------------------- #
	# The function genpwedata simulates a balanced clinical trial with
	# 'sample_size' subjects per arm using a piecewise exponential
	# distribution.  'lambda_vec' is the vector of exponential hazard
	# rates for each interval.  'time_vec' identifies the time
	# boundaries of each interval. 'hazard_ratio' defines the
	# relationship between group hazards in each interval.  At present
	# a constant hazard ratio is used. 'censor_value' is the value when
    # right censoring occurs.  As of 4/11/2018, genpwedata only
	# generates data with right censoring.  Random right censoring
	# is not incorporated.
	#
	# In the code below time1, lambda_vec, test1, status1, etc. are
	# data for the control goup.
	# In the code below time2, lambda_vec2, test2, status2, etc. are
	# data for the experimental group.
    # --------------------------------------------------------------- #

	# The user specified interval hazard vector, lambda_vec, and the
	# user specified vector of time cut points, are used to
	# generate random event times from a piecewise exponential.
	# Before I generate the event times for the experimental group
	# I need to generate the hazards for the experimental group.
	# For both control and experimental groups, each interval can
	# have its own hazard; however, in every interval the hazard_ratio
    # is equal to the user specified hazard ratio.

    lambdavec2 <- hazard_ratio * lambda_vec

    # Create event times for both groups
    time1 <- eha::rpch(n = sample_size, cuts = time_vec, levels = lambda_vec)
    time2 <- eha::rpch(n = sample_size, cuts = time_vec, levels = lambdavec2)

    # Create variables needed for simulation.
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
		# Need to add explaination of indicator1 & indicator2.
        indicator1 <- rep(2, sample_size)
        indicator2 <- rep(3, sample_size)
    }
	#Create status variable if censor_value is not specified (in such a case
	# status = 1 for all observeations.)
    if (is.null(censor_value) == TRUE) {
        status1 <- rep(1, sample_size)
        status2 <- rep(1, sample_size)
		# Need to add explaination of indicator1 & indicator2.
        indicator1 <- rep(2, sample_size)
        indicator2 <- rep(3, sample_size)
    }

	#Take all data created above and put into a data frame that contains
	# the required variables.
    subjid <- seq(from = 1, to = 2 * sample_size)
    trt <- c(rep(0, sample_size), rep(1, sample_size))
    event_time <- c(time1, time2)
    status <- c(status1, status2)
    indicator <- c(indicator1, indicator2)

    gendata <- data.frame(subjid, trt, event_time, status, indicator)
    colnames(gendata) <- c("id", "treatment", "event_time", "status", "indicator")

    return(gendata)
}


#' Log-likelihood function for two-arm trial with historical data using Piece-wise
#' Exponential (pwe) distribution.
#'
#' \code{pwe_loglike()} function used only used internally by
#' \code{pwetrialsimulator()} function to estimate pwe model parameters
#' when clinical trial involves experimental and control groups as well as historical
#' control data.  The pwe log-likelihood is calculated by using \code{eha::dpch()}.
#' pwe_loglike() should not be called directly by user.
#'
#' @param params  A vector of p + 1 parameters.  The last element is
#'   the log hazard ratio.  The first p elements are the estimated log hazards
#'   within each section created by \code{time_vec}.  If \code{time_vec} has
#'   m elements, then p = m + 1.
#' @param time_vec Vector of cut points used to divide time-scale into sections
#'   within which hazard is constant.
#' @param randdata  Dataset of randomly generated trial data.
#' @param hist_data Dataset of historical data.
#' @param a0 Power prior parameter: 0 implies historical data is ignored and 1 implies
#'   all information in historical data is used.
#'
#' @return \code{pwe_loglike()} returns a value of the loglikelihood function
#'   given a set of pwe parameters, randomly generated trial data, and observed
#'   historical data.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
pwe_loglike <- function(params, time_vec, randdata, hist_data, a0) {

    # --------------------------------------------------------------- #
    #  This function calculates the piecewise exponential
    #  log-likelihood given a vector of parameter values, a vector
    #  of time cutpoints to define intervals, a dataset of randomized
    #  trial data (two arms, no covariates beyond treatment), and a
    #  dataset of historical control data.
    #  The hazards within each interval is assumed to be the same in
    #  both control groups.
    #  For each interval the hazard in the experimental group is equal
    #  to the hazard ratio times the control hazard for that interval.
    #  The parameter vector is of length p+1, implying p time
    #  intervals and one parameter for the hazard ratio.  The
    #  parameters are on the log scale, so the first p parameters
    #  are log hazards for the time intervals and the last parameter
    #  is the log hazard ratio.
    # --------------------------------------------------------------- #

    # Extract the vector of hazards for control group from the parameter vector.
    # Note the parameters are on the log scale.
    params_c_i <- exp(params[-length(params)])
    # Extract the vector of hazards for the experimental group from the parameter vector.
    params_e_i <- exp(params[-length(params)] + params[length(params)])

    # Calculate the log-likelihood values for all randomized observations.
    ll_R <- (randdata$status == 1 & randdata$treatment == 0) * eha::dpch(randdata$event_time, cuts = time_vec,
                                                                         levels = params_c_i, log = TRUE) +
            (randdata$status == 0 & randdata$treatment == 0) * eha::ppch(randdata$event_time, cuts = time_vec,
                                                                         levels = params_c_i, lower.tail = FALSE,
                                                                         log.p = TRUE) +
            (randdata$status == 1 & randdata$treatment == 1) * eha::dpch(randdata$event_time, cuts = time_vec,
                                                                         levels = params_e_i, log = TRUE) +
            (randdata$status == 0 & randdata$treatment == 1) * eha::ppch(randdata$event_time, cuts = time_vec,
                                                                         levels = params_e_i, lower.tail = FALSE,
                                                                         log.p = TRUE)

    # Calculate the loglikelihood values for all historical control observations.
    ll_H <- (hist_data$status == 1) * eha::dpch(hist_data$event_time, cuts = time_vec, levels = params_c_i,
                                                log = TRUE) +
            (hist_data$status == 0) * eha::ppch(hist_data$event_time, cuts = time_vec, levels = params_c_i,
                                                lower.tail = FALSE, log.p = TRUE)

    # Calculate the overall log likelihood by adding the randomized log-likelihood to the historical control
    # log-likelihood by a0, where a0 is the power prior parameter.  This a0 value is defined by the
    # user and not estimated via object function optimization.
    ll <- sum(ll_R) + a0 * sum(ll_H)

    # Return the sum of all individual elements to the negative log-likelihood
    return(-ll)
}


#' Simulate a single randomized trial using a piece-wise exponential outcome and information from
#' historical controls.
#'
#' \code{pwetrialsimulator()} function only used internally by
#' \code{pwe_sim()} function to run a single trial simulation involving historical
#' control data and a piece-wise exponential outcome.
#'
#' The simulation of a trial with a piece-wise exponential outcome involving historical control data returns
#' an estimate of the hazard ratio as well as an estimate of the log hazard ratio variance.
#' Finally the simulation returns an indication of whether or not the simulated trial led to
#' a rejection of the null hypothesis (1) or not (0).
#'
#' \code{pwetrialsimulator()} should not be called directly by user.
#'
#' @param sample_size_val Number of subjects per arm.
#' @param hist_data Dataset of historical data.  Historical datasets must have 5 columns:
#'   id, treatment, event_time, status, and indicator.  The value of treatment should be 0.  The
#'   values of event_time must be positive.  The values of status must be 0 (right
#'   censored event) or 1 (observed event).  The value of indicator must be 2.
#' @param lambda_vec_val  Set of lambdas for the randomized control arm that will be
#'   passed to \code{eha::rpch()} through the levels parameter.
#' @param time_vec_val  Set of cutpoints for the simulated trial data that will be
#'   passed to \code{eha::rpch()} through the cuts parameter.
#' @param hazard_ratio_val Desired Hazard Ratio between randomized experimental and control groups.
#' @param censor_value Value at which time-to-event data are right censored.
#' @param a0_val A power prior parameter ranging from 0 to 1, where 0
#'   implies no information from historical data should be used, 1 implies all of
#'   the information from historical data should be used.  A value between 0 and 1
#'   implies that a proportion of the information from historical data will be used.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#'
#' @return \code{pwetrialsimulator()} returns a vector of simulation results. The
#'   first element is an estimated hazard ratio, the second element is the estimated
#'   variance of the log hazard ratio, and the third element is a 0/1 variable indicator
#'   whether or not the trial rejected the null hypothesis (1) or failed to reject
#'   the null hypothesis (0).
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
pwetrialsimulator <- function(sample_size_val, hist_data, lambda_vec_val, time_vec_val, hazard_ratio_val, censor_value,
    a0_val, alpha) {

  # --------------------------------------------------------------- #
  #  This function simulates a two-arm Bayesian trial where
  #  historical data is utilized in the parameter estimation.
  # --------------------------------------------------------------- #

  # First, Generate pwe trial data given the user defined trial characteristics.
    sampleranddata <- genpwedata(sample_size = sample_size_val, lambda_vec = lambda_vec_val, hazard_ratio = hazard_ratio_val,
        time_vec = time_vec_val, censor_value = censor_value)
    if (sum(sampleranddata$event_time == censor_value) == dim(sampleranddata)[1]) {
        stop("Simulated trial data must have at least one observation that is not right censored.")
    }

    # I have a problem to solve, what if I have an interval that has no events?
    # To solve this problem, I check to see if each time interval has enough events among
    # controls and experimental group.
    # Using the variable 'indicator', which is equal to 2 for controls and 3 for
    # experimental group, I create a table of indicator value by interval that identifies the
    # number of events that have occurred in each interval by 'indicator' combination.
    # If each cell of this table does not contain at least 2 events, then I reduce the number
    # of intervals by one, reconstruct the table and reassess.  The process stops if
    # all table cells contain at least 2 events or the number of intervals equals two.
    # Note that historical controls and randomized controls both need to have 'indicator'
    # set to 2, so both control groups will be combined into the same row of the assessment
    # table.  If the historical controls have a different value for 'indicator' then
    # the assessment table will have three rows and the number of intervals will be
    # heavily influenced by the historical data.
    # Regardless of the number of rows in the assessment table, each cell in the all_trts
    # by interval_test table needs to have at least 2 events.

    # First, get a list of event times for all observed events.
    all_events <- c(hist_data$event_time[hist_data$status == 1], sampleranddata$event_time[sampleranddata$status ==
                                                                                             1])
    # Next, get a list of the indicator values for all observed events.
    all_trts <- c(hist_data$indicator[hist_data$status == 1], sampleranddata$indicator[sampleranddata$status == 1])

    # Set the flag used to indicator when an acceptable set of time intervals are identified.
    flag1 <- 0

    # First check the user defined intervals for sufficient number of events in each interval
    # Note, findInterval uses time_vec_val to divide the real line into p + 1 intervals,
    # where p is equal to the number of elements in time_vec_val.  The return value of
    # findInterval() is 0,1,2,3,etc., where 0 is the interval 0 to the first element of
    # time_vec_val, 1 is the interval between the first and second elements of time_vec_val,
    # etc.
    interval_test <- findInterval(all_events, time_vec_val)
    test_table <- table(all_trts, interval_test)
    for (m in 1:dim(test_table)[1]) {
      for (n in 1:dim(test_table)[2]) {
        if (test_table[m, n] < 2 & flag1 == 0) flag1 <- 1
      }
    }

    # If user defined set of time cutoffs are sufficient use them.
    if (flag1 == 0) new_breaks <- time_vec_val

    #If user defined intervals do not create an assessment table with at least 5 events, then
    # start with a set of cutpoints equal in number to the length of the user defined set, but
    # spaced so that the number of events in each interval are about the same.
    if (flag1 == 1){
      new_breaks <- NULL
      splits <- length(time_vec_val) + 1 #Ensures the number of initial new breaks with equal length(time_vec_val)
      repeat {
        #For present number of cutpoints, identify cutpoints with equal proportions of
        # events within time intervals created by cutpoints.
        temp_breaks <- stats::quantile(all_events, prob = seq(from = 1/splits, to = 1 - 1/splits, by = 1/splits))
        #For all events, find the interval number (interval 1, interval 2, etc.) that
        # contains the event.
        interval_test <- findInterval(all_events, temp_breaks)
        # Create a table of group source (historical controls, randomized controls, or
        # or experimental group) by interval counts.
        test_table <- table(all_trts, interval_test)
        # If each cell in test_table contains at least 5 events, then cut point
        # selection is used to estimate hazard ratio.
        for (m in 1:dim(test_table)[1]) {
          for (n in 1:dim(test_table)[2]) {
            if (test_table[m, n] < 2 & flag1 == 1)
              flag1 <- 0
          }
        }

        # If we only have two intervals, stop.  Note: splits equals the number of intervals not number of cutpoints.
        if (flag1 == 1 | splits == 2) {
          break
        } else {
          splits <- splits - 1
          flag1 <- 1
        }
      }
      new_breaks <- as.vector(temp_breaks)
    }
    new_breaks2 <- c(0, new_breaks, censor_value)
    # The values in new_breaks create length(new_breaks)+1 time intervals. Each must have a
    # modeled hazard parameter.  Also the last parameter passed to the model fitting routine
    # will be the log hazard ratio.

    # Generate initial values for your call to optim()
    chk_model <- eha::pchreg(survival::Surv(event_time, status) ~ treatment, data = sampleranddata, cuts = new_breaks2)
    rc_hazards <- c(chk_model$hazards)
    rand_hr <- exp(chk_model$coefficients)

    init.parms <- as.vector(c(log(rc_hazards), log(rand_hr)))

    # Generate the Bayesian CLT based parameter estimates needed for inference on hazard ratio.
    fitmod <- stats::optim(par = init.parms, fn = pwe_loglike, time_vec = new_breaks, randdata = sampleranddata, hist_data = hist_data,
                    a0 = a0_val, method = "Nelder-Mead", hessian = TRUE)

    #Extract model parameters and statistics
    hr_loc <- length(init.parms)
    covar_mat <- solve(fitmod$hessian)

    hazard_ratio <- exp(fitmod$par[hr_loc])

    lower_pwe_hazard_ratio <- exp(fitmod$par[hr_loc] - stats::qnorm(1 - alpha/2) * sqrt(covar_mat[hr_loc, hr_loc]))
    upper_pwe_hazard_ratio <- exp(fitmod$par[hr_loc] + stats::qnorm(1 - alpha/2) * sqrt(covar_mat[hr_loc, hr_loc]))

    # Make a decision about the simulated trial, reject or fail to reject null hypothesis.
    reject <- ifelse(((lower_pwe_hazard_ratio > 1) | (upper_pwe_hazard_ratio < 1)), 1, 0)
    output <- c(hazard_ratio, covar_mat[hr_loc, hr_loc], reject)

    # Return the mean ratio, the estimated variance of the log mean ratio, and the trial decision.
    names(output) <- c("pwe_hazard_ratio", "log_hr_var", "reject")
    return(output)

}


#' Simulate a single randomized trial using a piece-wise exponential outcome but not including any
#' information from historical controls.
#'
#' \code{pwetrialsimulatornohist()} function only used internally by
#' \code{simple_pwe_sim()} function to estimate piece-wise exponential model parameters
#' when clinical trial involves experimental and control groups but no historical control
#' data.
#'
#' The simulation of a trial with a piece-wise exponential outcome involving no historical control
#' data returns an estimate of the hazard ratio as well as an estimate of the log hazard ratio variance.
#' Finally the simulation returns an indication of whether or not the simulated trial led to
#' a rejection of the null hypothesis (1) or not (0).
#'
#' \code{pwetrialsimulatornohist()} should not be called directly by user.
#'
#' @param sample_size_val Number of subjects per arm.
#' @param lambda_vec_val  Set of lambdas for the randomized control arm that will be
#'   passed to \code{eha::rpch()} through the levels parameter.
#' @param time_vec_val  Set of cutpoints for the simulated trial data that will be
#'   passed to \code{eha::rpch()} through the cuts parameter.
#' @param hazard_ratio_val Desired Hazard Ratio between randomized experimental and control groups.
#' @param censor_value Value at which time-to-event data are right censored.
#' @param alpha A number ranging between 0 and 1 that defines the acceptable Type 1
#'   error rate. Default is 0.05.
#'
#' @return \code{pwetrialsimulatornohist()} returns a vector of simulation results. The
#'   first element is an estimated hazard ratio, the second element is the estimated
#'   variance of the log hazard ratio, and the third element is a 0/1 variable indicator
#'   whether or not the trial rejected the null hypothesis (1) or failed to reject
#'   the null hypothesis (0).
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
pwetrialsimulatornohist <- function(sample_size_val, lambda_vec_val, time_vec_val, hazard_ratio_val, censor_value, alpha) {

    # --------------------------------------------------------------- #
    #  This function simulates a two-arm Bayesian trial where no
    #  historical data is utilized in the parameter estimation.
    #  Only randomized data is used in parameter estimation.
    # --------------------------------------------------------------- #

    # Generate pwe trial data
    sampleranddata <- genpwedata(sample_size = sample_size_val, lambda_vec = lambda_vec_val, hazard_ratio = hazard_ratio_val,
                                 time_vec = time_vec_val, censor_value = censor_value)
    if (sum(sampleranddata$event_time == censor_value) == dim(sampleranddata)[1]) {
      stop("Simulated trial data must have at least one observation that is not right censored.")
    }

    # I have a problem to solve, what if I have an interval that has no events?
    # To solve this problem, I check to see if each time interval has enough events among
    # controls and experimental group.
    # Using the variable 'indicator', which is equal to 2 for controls and 3 for
    # experimental group, I create a table of indicator value by interval that identifies the
    # number of events that have occurred in each interval by 'indicator' combination.
    # If each cell of this table does not contain at least 2 events, then I reduce the number
    # of intervals by one, reconstruct the table and reassess.  The process stops if
    # all table cells contain at least 2 events or the number of intervals equals two.
    # Note that historical controls and randomized controls both need to have 'indicator'
    # set to 2, so both control groups will be combined into the same row of the assessment
    # table.  If the historical controls have a different value for 'indicator' then
    # the assessment table will have three rows and the number of intervals will be
    # heavily influenced by the historical data.
    # Regardless of the number of rows in the assessment table, each cell in the all_trts
    # by interval_test table needs to have at least 2 events.

    # First, get a list of event times for all observed events.
    all_events <- c(sampleranddata$event_time[sampleranddata$status == 1])
    # Next, get a list of the indicator values for all observed events.
    all_trts <- c(sampleranddata$indicator[sampleranddata$status == 1])

    # Set the flag used to indicator when an acceptable set of time intervals are identified.
    flag1 <- 0

    #First check the user defined intervals for sufficient number of events in each interval
    # Note, findInterval uses time_vec_val to divide the real line into p + 1 intervals,
    # where p is equal to the number of elements in time_vec_val.  The return value of
    # findInterval() is 0,1,2,3,etc., where 0 is the interval 0 to the first element of
    # time_vec_val, 1 is the interval between the first and second elements of time_vec_val,
    # etc.
    interval_test <- findInterval(all_events, time_vec_val)
    test_table <- table(all_trts, interval_test)
    for (m in 1:dim(test_table)[1]) {
      for (n in 1:dim(test_table)[2]) {
        if (test_table[m, n] < 2 & flag1 == 0) flag1 <- 1
      }
    }

    # If user defined set of time cutoffs are sufficient use them.
    if (flag1 == 0) new_breaks <- time_vec_val

    # If user defined intervals do not create an assessment table with at least 5 events, then
    # start with a set of cutpoints equal in number to the length of the user defined set, but
    # spaced so that the number of events in each interval are about the same.
    if (flag1 == 1){
      new_breaks <- NULL
      splits <- length(time_vec_val) + 1 #Ensures the number of initial new breaks with equal length(time_vec_val)
      repeat {
        # For present number of cutpoints, identify cutpoints with equal proportions of
        # events within time intervals created by cutpoints.
        temp_breaks <- stats::quantile(all_events, prob = seq(from = 1/splits, to = 1 - 1/splits, by = 1/splits))
        # For all events, find the interval number (interval 1, interval 2, etc.) that
        # contains the event.
        interval_test <- findInterval(all_events, temp_breaks)
        # Create a table of group source (historical controls, randomized controls, or
        # or experimental group) by interval counts.
        test_table <- table(all_trts, interval_test)
        # If each cell in test_table contains at least 5 events, then cut point
        # selection is used to estimate hazard ratio.
        for (m in 1:dim(test_table)[1]) {
          for (n in 1:dim(test_table)[2]) {
            if (test_table[m, n] < 2 & flag1 == 1)
              flag1 <- 0
          }
        }

        # If we only have two intervals, stop.  Note: splits equals the number of intervals not number of cutpoints.
        if (flag1 == 1 | splits == 2) {
          break
        } else {
          splits <- splits - 1
          flag1 <- 1
        }
      }
      new_breaks <- as.vector(temp_breaks)
    }
    new_breaks2 <- c(0, new_breaks, censor_value)
    # The values in new_breaks create length(new_breaks)+1 time intervals. Each must have a
    # modeled hazard parameter.  Also the last parameter passed to the model fitting routine
    # will be the log hazard ratio.

    # Generate initial values for your call to optim()
    chk_model <- eha::pchreg(survival::Surv(event_time, status) ~ treatment, data = sampleranddata, cuts = new_breaks2)

    # Extract model parameters and statistics
    rc_hazards <- c(chk_model$hazards)
    rand_hr <- exp(chk_model$coefficients)

    hazard_ratio <- exp(chk_model$coefficients)

    lower_pwe_hazard_ratio <- exp(chk_model$coefficients - stats::qnorm(1 - alpha/2) * sqrt(chk_model$var))
    upper_pwe_hazard_ratio <- exp(chk_model$coefficients + stats::qnorm(1 - alpha/2) * sqrt(chk_model$var))

    # Make a decision about the simulated trial, reject or fail to reject null hypothesis.
    reject <- ifelse(((lower_pwe_hazard_ratio > 1) | (upper_pwe_hazard_ratio < 1)), 1, 0)
    output <- c(hazard_ratio, chk_model$var, reject)

    # Return the hazard ratio, the estimated variance of the log hazard ratio, and the trial decision.
    names(output) <- c("pwe_hazard_ratio", "log_hr_var", "reject")
    return(output)

}


#' Repeated Two Arm Bayesian Clinical Trial Simulation with Historical Data and
#' Piece-wise Exponential Outcome.
#'
#' \code{pwe_sim()} function only used internally by \code{historic_sim()}
#' function to run a set of trial simulations involving historical
#' control data and a piece-wise exponential outcome.  User defined simulation parameters are
#' used to generate a set of trial scenarios.  Each scenario is simulated multiple
#' times and then means are taken to calculate estimates of power, hazard ratios,
#' and other user requested study summary statistics like variance of hazard
#' ratio, bias (on hazard ratio scale), and mse (on hazard ratio scale).
#' The number of repeated simulations is defined by the user.
#'
#' \code{pwe_sim()} should not be called directly by user.
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
#' @param effect_vals A vector of hazard ratios (randomized experimental over control),
#'   all of which must be positive.
#' @param rand_control_diff For piece-wise exponential outcomes this is a vector of hazard ratios
#'   (randomized controls over historical controls) that represent differences
#'   between randomized and historical controls.
#' @param hist_control_data A dataset of historical data.  Default is \code{NULL}.
#'   Historical datasets must have 4 columns: id, treatment, event_time, status, and
#'   indicator.  The value of treatment should be 0.  The values of event_time must
#'   be positive.  The values of status must be 0 (right censored event) or
#'   1 (observed event).  The value of indicator must be 2.
#' @param time_vec_val A vector of time values which are used to create time periods
#'   within which the exponential hazard is constant.  Only used for piecewise
#'   exponential models.  Default is \code{NULL}.
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
#' @return \code{pwe_sim()} returns an S3 object of class bayes_ctd_array.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
pwe_sim <- function(trial_reps=100, subj_per_arm, a0_vals, effect_vals,
                    rand_control_diff, hist_control_data, time_vec_val,
                    censor_value, alpha=0.05, get_var=FALSE, get_bias=FALSE,
                    get_mse=FALSE, quietly=TRUE) {

    # --------------------------------------------------------------- #
    # For a set of user specified scenarios (defined by combinations
    # of user specified parameters), simulate "trial_reps" trials
    # and estimate power, hazard ratio estimate, and if requested by user:
    # variance of hazard ratio, bias, and mse.  Using a piece-wise
    # exponential oucome and incorporating data from historical controls.
    # --------------------------------------------------------------- #
    time_vec_val2 <- c(0, time_vec_val, censor_value)

    # Need to take the historical data and generate distributional parameter estimates
    hist_data <- hist_control_data
    chk_model <- eha::pchreg(survival::Surv(event_time, status) ~ 1, data = hist_data, cuts = time_vec_val2)
    hc_hazards <- c(chk_model$hazards)

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
        # Need to adjust the randomized control lambda parameters given the historical control lambda parameters and the
        # hazard ratios given in rand_control_diff
        adj_lambda <- rand_control_diff[diffs] * hc_hazards

        for (effvals in 1:length(effect_vals)) {
            for (a0vals in 1:length(a0_vals)) {
                for (sizes in 1:length(subj_per_arm)) {
                  if (!quietly){
                    cat("\r", c(subj_per_arm[sizes], a0_vals[a0vals], effect_vals[effvals], rand_control_diff[diffs]))
                  }
                  # For each combination of rand_control_diff, effect_vals, a0_val, and subj_per_arm, simulate the trial
                  #trial_reps times and then calculate the mean reject rate to estimate power.  For bias, work on the
                  #hazard ratio scale and take the mean of all differences between estimated hazard ratios and the
                  #true hazard ratio.  For mse, calculate the mean of squared differences between the
                  #estimated hazard ratios and the true hazard ratio value.
                  collect <- matrix(rep(0, 3 * trial_reps), ncol = 3)
                  for (k in 1:trial_reps) {
                    # sample_size_val will be equal to both arms
                    collect[k, ] <- pwetrialsimulator(sample_size_val = subj_per_arm[sizes], hist_data = hist_data,
                      lambda_vec_val = adj_lambda, time_vec_val = time_vec_val, hazard_ratio_val = effect_vals[effvals],
                      censor_value = censor_value, a0_val = a0_vals[a0vals], alpha = alpha)
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

    #Lines 637 through 940 simply apply names to the dimensions of array created by the
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
#' Piece-wise Exponential Outcome.
#'
#' \code{simple_pwe_sim()} function only used internally by \code{simple_sim()}
#' function to run a set of trial simulations involving no historical
#' control data and a piece-wise exponential outcome.  User defined simulation parameters are
#' used to generate a set of trial scenarios.  Each scenario is simulated multiple
#' times and then means are taken to calculate estimates of power, hazard ratios,
#' and other user requested study summary statistics like variance of hazard
#' ratio, bias (on hazard ratio scale), and mse (on hazard ratio scale).
#' The number of repeated simulations is defined by the user.
#'
#' \code{simple_pwe_sim()} should not be called directly by user.
#'
#' @param trial_reps Number of trials to replicate within each combination of
#'   a0_val, subj_per_arm, effect_vals, and rand_control_diff.  As the number
#'   of trials increases, the precision of the estimate will increase. Default is
#'   100.
#' @param subj_per_arm A vector of sample sizes, all of which must be positive
#'   integers.
#' @param effect_vals A vector of hazard ratios (randomized experimental over control),
#'   all of which must be positive.
#' @param time_vec_val Set of cutpoints for the simulated trial data that will be
#'   passed to \code{eha::rpch()} through the cuts parameter.
#' @param rc_hazards Set of lambdas for the randomized control arm that will be
#'   passed to \code{eha::rpch()} through the levels parameter.
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
#' @return \code{simple_pwe_sim()} returns an S3 object of class bayes_ctd_array.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
simple_pwe_sim <- function(trial_reps=100, subj_per_arm, effect_vals, time_vec_val,
                           rc_hazards, censor_value, alpha=0.05, get_var=FALSE,
                           get_bias=FALSE, get_mse=FALSE, quietly=TRUE) {

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
                  #hazard ratio scale and take the mean of all differences between estimated hazard ratios and the
                  #true hazard ratio.  For mse, calculate the mean of squared differences between the estimated
                  #hazard ratios and the true hazard ratio value.  Note that rand_control_diff is set to 1 and
                  #a0_val is set to 0.
                  collect <- matrix(rep(0, 3 * trial_reps), ncol = 3)
                  for (k in 1:trial_reps) {
                    # sample_size_val will be equal in both arms
                    collect[k, ] <- pwetrialsimulatornohist(sample_size_val = subj_per_arm[sizes], lambda_vec_val = rc_hazards,
                      time_vec_val = time_vec_val, hazard_ratio_val = effect_vals[effvals], censor_value = censor_value,
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

    #Lines 1073 through 1376 simply apply names to the dimensions of array created by the
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

