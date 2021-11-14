
#' Checks for Errors when Outcome is Piece-wise Exponential.
#'
#' \code{pwe_error_checks()} function used only used internally by
#' \code{historic_sim()} function to check for proper input.  Not to be called
#' directly by user.
#'
#' @param effect_vals See help page for \code{historic_sim()}.
#' @param hist_control_data See help page for \code{historic_sim()}.
#' @param rand_control_diff See help page for \code{historic_sim()}.
#' @param time_vec See help page for \code{historic_sim()}.
#' @param censor_value See help page for \code{historic_sim()}.
#' @param alpha See help page for \code{historic_sim()}.

#'
#' @return \code{pwe_error_checks()} returns messages when
#' \code{historic_sim()} function inputs are incorrectly specified.
#' Not to be called directly by user.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
pwe_error_checks <- function(effect_vals, hist_control_data, rand_control_diff, time_vec, censor_value, alpha) {

    hist_chk <- class(hist_control_data)
    rcp_chk <- is.null(rand_control_diff)
    tv_chk <- is.null(time_vec)

    # Need to check that Effect_vals is a list of positive numbers
    for (eff_val in effect_vals) {
        if (!is.numeric(eff_val))
            stop("historic_sim() requires effect_vals to be numeric.")
        if (eff_val <= 0)
            stop("For pwe, historic_sim() requires effect_vals to represent Hazard Ratios which must be positive numbers")
    }
    # Need to check that hist_control_data is a data frame and has the correct columns.
    if (hist_chk == "data.frame") {
        colnamevals <- colnames(hist_control_data)
        colnamechk <- (colnamevals == c("id", "treatment", "event_time", "status", "indicator"))
        if (sum(colnamechk) != 5) {
            stop("historic_sim() requires hist_control_data to have columns: id, treatment, event_time, status, and indicator.")
        }
        if (sum(colnamechk) == 5) {
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
            if (!is.numeric(hist_control_data$indicator)) {
                stop("historic_sim() requires hist_control_data$indicator to be a vector of 2s.")
            }
            if (2*dim(hist_control_data)[1] != sum(hist_control_data$indicator)) {
                stop("historic_sim() requires hist_control_data$indicator to be a vector of 2s.")
            }
        }
    }
    # If not NULL, need to check that rand_control_diff is positive.
    if (rcp_chk == FALSE) {
        for (Rand_cp in rand_control_diff) {
            if (!is.numeric(Rand_cp))
                stop("historic_sim() requires rand_control_diff to be numeric.")
            if (Rand_cp <= 0)
                stop("historic_sim() requires rand_control_diff to be Hazard Ratios, so they must be positive.")
        }
    }
    if (tv_chk == TRUE) {
        stop("historic_sim() requires at least one value in time_vec for a pwe model.")
    }
    if (tv_chk == FALSE) {
        for (tvcp in time_vec) {
            if (!is.numeric(tvcp))
                stop("historic_sim() requires time_vec to be numeric.")
            if (tvcp <= 0)
                stop("historic_sim() requires time_vec to be positive.")
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


#' Checks for Errors when Outcome is Piece-wise Exponential.
#'
#' \code{pwe_error_checks_simple()} function used only used internally by
#' \code{simple_sim()} function to check for proper input.  Not to be called
#' directly by user.
#'
#' @param effect_vals See help page for \code{simple_sim()}.
#' @param time_vec See help page for \code{simple_sim()}.
#' @param control_parms See help page for \code{simple_sim()}.
#' @param censor_value See help page for \code{simple_sim()}.
#' @param alpha See help page for \code{simple_sim()}.

#'
#' @return \code{pwe_error_checks_simple()} returns messages when
#' \code{simple_sim()} function inputs are incorrectly specified.
#' Not to be called directly by user.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
pwe_error_checks_simple <- function(effect_vals, time_vec, control_parms, censor_value, alpha) {

    tv_chk <- is.null(time_vec)
    tv_len <- length(time_vec)
    rv_chk <- is.null(control_parms)
    rv_len <- length(control_parms)

    # Need to check that Effect_vals is a list of positive numbers
    for (eff_val in effect_vals) {
        if (!is.numeric(eff_val))
            stop("simple_sim() requires effect_vals to be numeric.")
        if (eff_val <= 0)
            stop("For pwe, simple_sim() requires effect_sizes to represent Hazard Ratios which must be positive numbers")
    }
    if (tv_chk == TRUE) {
        stop("simple_sim() requires at least one value in time_vec for a pwe model.")
    }
    if (tv_chk == FALSE) {
        for (tvcp in time_vec) {
            if (!is.numeric(tvcp))
                stop("simple_sim() requires time_vec to be numeric.")
            if (tvcp <= 0)
                stop("simple_sim() requires time_vec to be positive.")
        }
    }
    if (rv_chk == TRUE) {
        stop("simple_sim() requires at least two value in control_parms for a pwe model.")
    }
    if (rv_chk == FALSE) {
        if (length(control_parms) < 2) {
            stop("simple_sim() requires at least two value in control_parms for a pwe model.")
        }
        if (rv_len != (tv_len + 1)) {
            stop("simple_sim() requires control_parms to have one more element than time_vec")
        }
        for (rvcp in control_parms) {
            if (!is.numeric(rvcp))
                stop("simple_sim() requires control_parms to be numeric.")
            if (rvcp <= 0)
                stop("simple_sim() requires control_parms to be positive.")
        }
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
