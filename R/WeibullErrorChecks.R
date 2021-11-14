
#' Checks for Errors when Outcome is Weibull.
#'
#' \code{weibull_error_checks()} function used only used internally by
#' \code{historic_sim()} function to check for proper input.  Not to be called
#' directly by user.
#'
#' @param effect_vals See help page for \code{historic_sim()}.
#' @param hist_control_data See help page for \code{historic_sim()}.
#' @param rand_control_diff See help page for \code{historic_sim()}.
#' @param censor_value See help page for \code{historic_sim()}.
#' @param alpha See help page for \code{historic_sim()}.
#'
#' @return \code{weibull_error_checks()} returns messages when
#' \code{historic_sim()} function inputs are incorrectly specified.
#' Not to be called directly by user.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
weibull_error_checks <- function(effect_vals, hist_control_data, rand_control_diff, censor_value, alpha) {

    hist_chk <- class(hist_control_data)
    rcp_chk <- is.null(rand_control_diff)

    # Need to check that Effect_vals is a list of positive numbers
    for (eff_val in effect_vals) {
        if (!is.numeric(eff_val))
            stop("historic_sim() requires effect_vals to be numeric.")
        if (eff_val <= 0)
            stop("For Weibull, historic_sim() requires effect_vals to represent Hazard Ratios which must be positive numbers")
    }
    # Need to check that hist_control_data is a data frame and has the correct columns.
    if (hist_chk == "data.frame") {
        colnamevals <- colnames(hist_control_data)
        colnamechk <- (colnamevals == c("id", "treatment", "event_time", "status"))
        if (sum(colnamechk) != 4) {
            stop("historic_sim() requires hist_control_data to have columns: id, treatment, event_time, and status.")
        }
        if (sum(colnamechk) == 4) {
            if (!is.numeric(hist_control_data$event_time)) {
                stop("historic_sim() requires hist_control_data$event_time to be positive numeric data.")
            }
            if (is.numeric(hist_control_data$event_time)) {
                if (min(hist_control_data$event_time) <= 0) {
                  stop("historic_sim() requires hist_control_data$event_time to be positive numeric data.")
                }
            }
            if (!is.numeric(hist_control_data$treatment)) {
                stop("historic_sim() requires hist_control_data$treatment to be numeric 0/1 data.")
            }
            if (is.numeric(hist_control_data$treatment)) {
                trt_levels <- names(table(hist_control_data$treatment))
                if (length(trt_levels) > 2 | (trt_levels[1] != "0" & trt_levels[2] != "1")) {
                  stop("historic_sim() requires hist_control_data$treatment to be numeric 0/1 data.")
                }
            }
            if (!is.numeric(hist_control_data$status)) {
                stop("historic_sim() requires hist_control_data$status to be numeric 0/1 data.")
            }
            if (is.numeric(hist_control_data$status)) {
                trt_levels <- names(table(hist_control_data$status))
                if (length(trt_levels) > 2 | (trt_levels[1] != "0" & trt_levels[2] != "1")) {
                  stop("historic_sim() requires hist_control_data$status to be numeric 0/1 data.")
                }
            }
        }
    }
    # If not NULL, need to check that rand_control_diff is positive.
    if (rcp_chk == FALSE) {
        for (rand_cp in rand_control_diff) {
            if (!is.numeric(rand_cp))
                stop("historic_sim() requires rand_control_diff to be numeric.")
            if (rand_cp <= 0)
                stop("historic_sim() requires rand_control_diff to be Hazard Ratios, so they must be positive.")
        }
    }
    # Need to check that censor_value is NULL or a non-negative number
    if (!is.null(censor_value) == TRUE) {
        if (!is.numeric(censor_value))
            stop("historic_sim() requires censor_value to be numeric.")
        if (censor_value <= 0)
            stop("historic_sim() requires censor_value to be positive.")
    }
    # Need to check that alpha ranges between 0 and 1
    if (!is.null(alpha) == TRUE) {
        if (!is.numeric(alpha))
            stop("historic_sim() requires alpha to be numeric.")
        if (alpha <= 0 | alpha >= 1)
            stop("historic_sim() requires alpha to be between 0 and 1 but not equal to 0 or 1.")
    }
}


#' Checks for Errors when Outcome is Weibull.
#'
#' \code{weibull_error_checks_simple()} function used only used internally by
#' \code{simple_sim()} function to check for proper input.  Not to be called
#' directly by user.
#'
#' @param effect_vals See help page for \code{simple_sim()}.
#' @param control_parms See help page for \code{simple_sim()}.
#' @param censor_value See help page for \code{simple_sim()}.
#' @param alpha See help page for \code{simple_sim()}.

#'
#' @return \code{weibull_error_checks_simple()} returns messages when
#' \code{simple_sim()} function inputs are incorrectly specified.
#' Not to be called directly by user.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
weibull_error_checks_simple <- function(effect_vals, control_parms, censor_value, alpha) {

    chk_parm_lgth <- length(control_parms)

    # Need to check that effect_vals is a list of positive numbers
    for (eff_val in effect_vals) {
        if (!is.numeric(eff_val))
            stop("simple_sim() requires effect_vals to be numeric.")
        if (eff_val <= 0)
            stop("For Weibull, simple_sim() requires effect_vals to represent Hazard Ratios which must be positive numbers")
    }
    if (chk_parm_lgth != 2) {
        stop("simple_sim() requires two elements in control_parms, first=Weibull scale parameter for rweibull(), second= Weibull shape parameter for rweibull()")
    }
    if (chk_parm_lgth == 2) {
        if (!is.numeric(control_parms[1]))
            stop("simple_sim() requires Weibull scale parameter for controls to be numeric.")
        if (!is.numeric(control_parms[2]))
            stop("simple_sim() requires Weibull shape parameter for controls to be numeric.")
        if (control_parms[1] <= 0)
            stop("simple_sim() requires Weibull scale parameter for controls to be positive.")
        if (control_parms[2] <= 0)
            stop("simple_sim() requires Weibull shape parameter for controls to be positive.")
    }
    # Need to check that censor_value is NULL or a non-negative number
    if (!is.null(censor_value) == TRUE) {
        if (!is.numeric(censor_value))
            stop("simple_sim() requires censor_value to be numeric.")
        if (censor_value <= 0)
            stop("simple_sim() requires censor_value to be positive.")
    }
    # Need to check that alpha ranges between 0 and 1
    if (!is.null(alpha) == TRUE) {
        if (!is.numeric(alpha))
            stop("simple_sim() requires alpha to be numeric.")
        if (alpha <= 0 | alpha >= 1)
            stop("simple_sim() requires alpha to be between 0 and 1 but not equal to 0 or 1.")
    }
}
