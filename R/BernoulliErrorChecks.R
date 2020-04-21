
#' Checks for Errors when Outcome is Bernoulli.
#'
#' \code{bernoulli_error_checks()} function used only used internally by
#' \code{historic_sim()} function to check for proper input.  Not to be called
#' directly by user.
#'
#' @param effect_vals See help page for \code{historic_sim()}.
#' @param hist_control_data See help page for \code{historic_sim()}.
#' @param rand_control_diff See help page for \code{historic_sim()}.
#' @param alpha See help page for \code{historic_sim()}.

#'
#' @return \code{bernoulli_error_checks()} returns messages when
#' \code{historic_sim()} function inputs are incorrectly specified.
#' Not to be called directly by user.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
bernoulli_error_checks <- function(effect_vals, hist_control_data, rand_control_diff, alpha) {

    hist_chk <- class(hist_control_data)
    rcp_chk <- is.null(rand_control_diff)

    # Need to check that EffectSize is a list of positive numbers
    for (eff_val in effect_vals) {
        if (!is.numeric(eff_val))
            stop("historic_sim() requires effect_vals to be numeric.")
        if (eff_val <= 0)
            stop("For Bernoulli, historic_sim() requires effect_vals to represent Odds Ratios which must be positive numbers")
    }
    # Need to check that hist_control_data is a data frame and has the correct columns.
    if (hist_chk == "data.frame") {
        colnamevals <- colnames(hist_control_data)
        colnamechk <- (colnamevals == c("id", "treatment", "y"))
        if (sum(colnamechk) != 3) {
            stop("historic_sim() requires hist_control_data to have columns: ID, Treatment, and Y.")
        }
        if (sum(colnamechk) == 3) {
            if (!is.numeric(hist_control_data$y)) {
                stop("historic_sim() requires hist_control_data$Y to be numeric.")
            }
            if (!is.numeric(hist_control_data$treatment)) {
                stop("historic_sim() requires hist_control_data$treatment to be numeric 0/1 data.")
            }
            if (is.numeric(hist_control_data$treatment)) {
                trt_levels <- names(table(hist_control_data$treatment))
                if (length(trt_levels) > 2 | (trt_levels[1] != "0" & trt_levels[2] != "1")) {
                  stop("historic_sim() requires hist_control_data$Treatment to be numeric 0/1 data.")
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
                stop("historic_sim() requires rand_control_diff to be Odds Ratios, so they must be positive.")
        }
    }
    # Need to check that alpha ranges between 0 and 1
    if (!is.null(alpha) == TRUE) {
        if (!is.numeric(alpha))
            stop("historic_sim() requires alpha to be numeric.")
        if (alpha <= 0 | alpha >= 1)
            stop("historic_sim() requires alpha to be between 0 and 1 but not equal to 0 or 1.")
    }
}


#' Checks for Errors when Outcome is Bernoulli.
#'
#' \code{bernoulli_error_checks_simple()} function used only used internally by
#' \code{simple_sim()} function to check for proper input.  Not to be called
#' directly by user.
#'
#' @param effect_vals See help page for \code{simple_sim()}.
#' @param control_parms See help page for \code{simple_sim()}.
#' @param alpha See help page for \code{simple_sim()}.

#'
#' @return \code{bernoulli_error_checks_simple()} returns messages when
#' \code{simple_sim()} function inputs are incorrectly specified.
#' Not to be called directly by user.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
bernoulli_error_checks_simple <- function(effect_vals, control_parms, alpha) {

    chk_parm_lgth <- length(control_parms)

    # Need to check that EffectSize is a list of positive numbers
    for (eff_val in effect_vals) {
        if (!is.numeric(eff_val))
            stop("simple_sim() requires effect_vals to be numeric.")
        if (eff_val <= 0)
            stop("For Bernoulli, simple_sim() requires effect_vals to represent Odds Ratios which must be positive numbers")
    }

    if (chk_parm_lgth != 1) {
        stop("simple_sim() requires a single elements in control_parms which is the probability parameter for rbinom() among controls")
    }
    if (chk_parm_lgth == 1) {
        if (!is.numeric(control_parms[1]))
            stop("simple_sim() requires probability parameter for controls to be numeric.")
        if (control_parms[1] <= 0)
            stop("simple_sim() requires probability parameter for controls to be between 0 and 1.")
        if (control_parms[1] >= 1)
            stop("simple_sim() requires probability parameter for controls to be between 0 and 1.")
    }

    # Need to check that alpha ranges between 0 and 1
    if (!is.null(alpha) == TRUE) {
        if (!is.numeric(alpha))
            stop("simple_sim() requires alpha to be numeric.")
        if (alpha <= 0 | alpha >= 1)
            stop("simple_sim() requires alpha to be between 0 and 1 but not equal to 0 or 1.")
    }
}
