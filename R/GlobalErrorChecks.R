
#' Checks \code{historic_sim()} Function Call for Global Type Errors.
#'
#' \code{global_error_checks()} function used only used internally by \code{historic_sim()}
#' function to check for proper input.  Not to be called directly by user.
#'
#' @param outcome_type See help page for \code{historic_sim()}.
#' @param subj_per_arm See help page for \code{historic_sim()}.
#' @param a0_vals See help page for \code{historic_sim()}.
#' @param hist_control_data See help page for \code{historic_sim()}.
#' @param rand_control_diff See help page for \code{historic_sim()}.
#' @param get_var See help page for \code{historic_sim()}.
#' @param get_bias See help page for \code{historic_sim()}.
#' @param get_mse See help page for \code{historic_sim()}.
#'
#' @return \code{global_error_checks()} returns messages when \code{historic_sim()} function
#' inputs are incorrectly specified.  Not to be called directly by user.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
global_error_checks <- function(outcome_type, a0_vals, subj_per_arm, hist_control_data, rand_control_diff, get_var, get_bias,
    get_mse) {

    if (!(tolower(outcome_type) %in% c("weibull", "lognormal", "pwe", "gaussian", "bernoulli", "poisson"))) {
        stop("historic_sim only supports the following outcome distributions:\nWeibull, Lognormal,Piecewsie Exponential (PWE), Gaussian, Logistic, and Poisson.\noutcome_type values must be equal to weibull, lognormal, pwe, gaussian, logistic, or poisson")
    }

    # a0_val must be numeric and between 0 and 1 inclusively.
    for (pow_parm in a0_vals) {
        if (!is.numeric(pow_parm))
            stop("historic_sim requires a0_val to be numeric.")
        if (is.numeric(pow_parm) & (pow_parm < 0 | pow_parm > 1))
            stop("historic_sim requires a0_val to be between 0 and 1 inclusively.")
    }

    # subj_per_arm must be a list of non-negative values.
    temp_s <- as.integer(subj_per_arm)
    for (sam_siz in subj_per_arm) {
        if (!is.numeric(sam_siz))
            stop("historic_sim requires subj_per_arm to be numeric.")
        if (sam_siz < 1)
            stop("historic_sim requires at least 1 randomized subject per arm. Negative numbers or zero not allowed.")
    }
    sam_siz_chk <- sum(temp_s == subj_per_arm)
    if (sam_siz_chk != length(subj_per_arm))
        stop("historic_sim requires subj_per_arm to be an integer.")

    # hist_control_data must be a data frame or a list of numeric values.
    hist_chk <- class(hist_control_data)
    if (hist_chk != "data.frame") {
        stop("historic_sim requires hist_control_data to be a data frame.")
    }

    # Need to check that rand_control_diff is either NULL or a list of reasonable numbers
    rcp_chk <- is.null(rand_control_diff)
    if (rcp_chk == FALSE) {
        for (rand_cp in rand_control_diff) {
            if (!is.numeric(rand_cp))
                stop("historic_sim requires rand_control_diff to be numeric.")
        }
    }

    # Need to check that get_var is TRUE/FALSE
    if (get_var != TRUE & get_var != FALSE) {
        stop("historic_sim requires get_var to be either TRUE or FALSE")
    }

    # Need to check that get_bias is TRUE/FALSE
    if (get_bias != TRUE & get_bias != FALSE) {
        stop("historic_sim requires get_bias to be either TRUE or FALSE")
    }

    # Need to check that get_mse is TRUE/FALSE
    if (get_mse != TRUE & get_mse != FALSE) {
        stop("historic_sim requires get_mse to be either TRUE or FALSE")
    }

}


#' Checks \code{simple_sim()} Function Call for Global Type Errors.
#'
#' \code{global_error_checks_simple()} function used only used internally by
#' \code{simple_sim()} function to check for proper input.  Not to be called
#' directly by user.
#'
#' @param outcome_type See help page for \code{simple_sim()}.
#' @param subj_per_arm See help page for \code{simple_sim()}.
#' @param get_var See help page for \code{simple_sim()}.
#' @param get_bias See help page for \code{simple_sim()}.
#' @param get_mse See help page for \code{simple_sim()}.
#'
#' @return \code{global_error_checks_simple()} returns messages when
#' \code{simple_sim()} function inputs are incorrectly specified.  Not to be called
#' directly by user.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
global_error_checks_simple <- function(outcome_type, subj_per_arm, get_var, get_bias, get_mse) {

    if (!(tolower(outcome_type) %in% c("weibull", "lognormal", "pwe", "gaussian", "bernoulli", "poisson"))) {
        stop("simple_sim only supports the following outcome distributions:\nWeibull, Lognormal,Piecewsie Exponential (PWE), Gaussian, Logistic, and Poisson.\noutcome_type values must be equal to weibull, lognormal, pwe, gaussian, logistic, or poisson")
    }

    # subj_per_arm must be a list of non-negative values.
    temp_s <- as.integer(subj_per_arm)
    for (sam_siz in subj_per_arm) {
        if (!is.numeric(sam_siz))
            stop("simple_sim requires subj_per_arm to be numeric.")
        if (sam_siz < 1)
            stop("simple_sim requires at least 1 randomized subject per arm. Negative numbers or zero not allowed.")
    }
    sam_siz_chk <- sum(temp_s == subj_per_arm)
    if (sam_siz_chk != length(subj_per_arm))
        stop("simple_sim requires subj_per_arm to be an integer.")

    # Need to check that get_var is TRUE/FALSE
    if (get_var != TRUE & get_var != FALSE) {
        stop("simple_sim requires get_var to be either TRUE or FALSE")
    }

    # Need to check that get_bias is TRUE/FALSE
    if (get_bias != TRUE & get_bias != FALSE) {
        stop("simple_sim requires get_bias to be either TRUE or FALSE")
    }

    # Need to check that get_mse is TRUE/FALSE
    if (get_mse != TRUE & get_mse != FALSE) {
        stop("simple_sim requires get_mse to be either TRUE or FALSE")
    }

}

#' Checks \code{print_table()} Function Call for Errors.
#'
#' \code{print_table_error_checks()} function used only internally by
#' \code{print_table()} function to check for proper input.  Not to be called
#' directly by user.
#'
#' @param bayes_ctd_robj See help page for \code{print_table()}.
#' @param measure, See help page for \code{print_table()}.
#' @param tab_type, See help page for \code{print_table()}.
#' @param subj_per_arm_val See help page for \code{print_table()}.
#' @param a0_val See help page for \code{print_table()}.
#' @param effect_val See help page for \code{print_table()}.
#' @param rand_control_diff_val See help page for \code{print_table()}.
#'
#' @return \code{print_table_error_checks()} returns messages when
#' \code{print_table} function inputs are incorrectly specified.  Not to be called
#' directly by user.
#'
#' @examples
#' #None.
#' @keywords internal
#' @noRd
print_table_error_checks <- function(bayes_ctd_robj, measure, tab_type, subj_per_arm_val, a0_val, effect_val, rand_control_diff_val) {

    if (!(tolower(measure) %in% c("power", "est", "var", "bias", "mse"))) {
        stop("print_table() requires measure to be equal to power, est, var, bias, mse")
    }

    if ((tolower(measure) == "var") & is.null(bayes_ctd_robj$data$var)){
        stop("Table of variance cannot be printed or plotted, because variance estimates were not requested in the trial simulation.")
    }
    if ((tolower(measure) == "bias") & is.null(bayes_ctd_robj$data$bias)){
        stop("Table of bias cannot be printed or plotted, because bias estimates were not requested in the trial simulation.")
    }
    if ((tolower(measure) == "mse") & is.null(bayes_ctd_robj$data$mse)){
       stop("Table of mse cannot be printed or plotted, because mse estimates were not requested in the trial simulation.")
    }


    if (!(toupper(tab_type) %in% c("WX|YZ", "WY|XZ", "WZ|XY", "XY|WZ", "XZ|WY", "YZ|WX", "ZX|WY", "XW|YZ", "YW|XZ", "YX|WZ", "ZW|XY", "ZX|WY", "ZY|WX"))) {
        stop("print_table() requires tab_type to be equal to WX|YZ, WY|XZ, WZ|XY, XY|WZ, XZ|WY, YZ|WX, ZX|WY, XW|YZ, YW|XZ, YX|WZ, ZW|XY, ZX|WY, or ZY|WX")
    }

    if (toupper(tab_type) == "WX|YZ") {
        if (is.null(effect_val) | is.null(rand_control_diff_val)) {
            stop("print_table() requires effect_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
        }
    }
    if (toupper(tab_type) == "WY|XZ") {
        if (is.null(a0_val) | is.null(rand_control_diff_val)) {
            stop("print_table() requires a0_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
        }
    }
    if (toupper(tab_type) == "WZ|XY") {
        if (is.null(a0_val) | is.null(effect_val)) {
            stop("print_table() requires a0_val and effect_val to be non-missing and equal to a value used in your simulation.")
        }
    }
    if (toupper(tab_type) == "XY|WZ") {
        if (is.null(subj_per_arm_val) | is.null(rand_control_diff_val)) {
            stop("print_table() requires subj_per_arm_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
        }
    }
    if (toupper(tab_type) == "XZ|WY") {
        if (is.null(subj_per_arm_val) | is.null(effect_val)) {
            stop("print_table() requires subj_per_arm_val and effect_val to be non-missing and equal to a value used in your simulation.")
        }
    }
    if (toupper(tab_type) == "YZ|WX") {
        if (is.null(subj_per_arm_val) | is.null(a0_val)) {
            stop("print_table() requires subj_per_arm_val and a0_val to be non-missing and equal to a value used in your simulation.")
        }
    }
    if (toupper(tab_type) == "ZX|WY") {
        if (is.null(subj_per_arm_val) | is.null(effect_val)) {
            stop("print_table() requires subj_per_arm_val and effect_val to be non-missing and equal to a value used in your simulation.")
        }
    }
    if (toupper(tab_type) == "XW|YZ") {
      if (is.null(effect_val) | is.null(rand_control_diff_val)) {
        stop("print_table() requires effect_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
      }
    }
    if (toupper(tab_type) == "YW|XZ") {
      if (is.null(a0_val) | is.null(rand_control_diff_val)) {
        stop("print_table() requires a0_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
      }
    }
    if (toupper(tab_type) == "YX|WZ") {
      if (is.null(subj_per_arm_val) | is.null(rand_control_diff_val)) {
        stop("print_table() requires subj_per_arm_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
      }
    }
    if (toupper(tab_type) == "ZW|XY") {
      if (is.null(a0_val) | is.null(effect_val)) {
        stop("print_table() requires a0_val and effect_val to be non-missing and equal to a value used in your simulation.")
      }
    }
    if (toupper(tab_type) == "ZX|WY") {
      if (is.null(subj_per_arm_val) | is.null(effect_val)) {
        stop("print_table() requires subj_per_arm_val and effect_val to be non-missing and equal to a value used in your simulation.")
      }
    }
    if (toupper(tab_type) == "ZY|WX") {
      if (is.null(subj_per_arm_val) | is.null(a0_val)) {
        stop("print_table() requires subj_per_arm_val and a0_val to be non-missing and equal to a value used in your simulation.")
      }
    }

    subj_per_arm_vec <- bayes_ctd_robj$subj_per_arm
    a0_val_vec <- bayes_ctd_robj$a0_vals
    effect_val_vec <- bayes_ctd_robj$effect_vals
    rand_control_diff_vec <- bayes_ctd_robj$rand_control_diff

    if (!is.null(subj_per_arm_val)) {
        if (is.na(match(subj_per_arm_val, subj_per_arm_vec))) {
            text_val <- "print_table() requires subj_per_arm_val to be equal to one of the following values: "
            for (i in 1:length(subj_per_arm_vec)) {
                text_val <- paste(text_val, subj_per_arm_vec[i], sep = "")
                if (i != length(subj_per_arm_vec)) {
                  text_val <- paste(text_val, ", ", sep = "")
                }
            }
            stop(text_val)
        }
    }

    if (!is.null(a0_val)) {
        if (is.na(match(a0_val, a0_val_vec))) {
            text_val <- "print_table() requires a0_val to be equal to one of the following values: "
            for (i in 1:length(a0_val_vec)) {
                text_val <- paste(text_val, a0_val_vec[i], sep = "")
                if (i != length(a0_val_vec)) {
                  text_val <- paste(text_val, ", ", sep = "")
                }
            }
            stop(text_val)
        }
    }

    if (!is.null(effect_val)) {
        if (is.na(match(effect_val, effect_val_vec))) {
            text_val <- "print_table() requires effect_val to be equal to one of the following values: "
            for (i in 1:length(effect_val_vec)) {
                text_val <- paste(text_val, effect_val_vec[i], sep = "")
                if (i != length(effect_val_vec)) {
                  text_val <- paste(text_val, ", ", sep = "")
                }
            }
            stop(text_val)
        }
    }

    if (!is.null(rand_control_diff_val)) {
        if (is.na(match(rand_control_diff_val, rand_control_diff_vec))) {
            text_val <- "print_table() requires rand_control_diff_val to be equal to one of the following values: "
            for (i in 1:length(rand_control_diff_vec)) {
                text_val <- paste(text_val, rand_control_diff_vec[i], sep = "")
                if (i != length(rand_control_diff_vec)) {
                  text_val <- paste(text_val, ", ", sep = "")
                }
            }
            stop(text_val)
        }
    }

}


#' Checks \code{print.bayes_ctd_array()} Function Call for Errors.
#'
#' \code{print_error_checks()} function used only internally by
#' \code{print.bayes_ctd_array()} function to check for proper input.  Not to be called
#' directly by user.
#'
#' @param bayes_ctd_robj See help page for \code{print_table()}.
#' @param measure, See help page for \code{print_table()}.
#' @param tab_type, See help page for \code{print_table()}.
#' @param subj_per_arm_val See help page for \code{print_table()}.
#' @param a0_val See help page for \code{print_table()}.
#' @param effect_val See help page for \code{print_table()}.
#' @param rand_control_diff_val See help page for \code{print_table()}.
#'
#' @return \code{print_error_checks()} returns messages when
#' \code{printbayes_ctd_array()} function inputs are incorrectly specified.  Not to be called
#' directly by user.
#'
#' @examples
#' #None.
#' @keywords internal
#' @noRd
print_error_checks <- function(bayes_ctd_robj, measure, tab_type, subj_per_arm_val, a0_val, effect_val, rand_control_diff_val) {

  if (!(tolower(measure) %in% c("power", "est", "var", "bias", "mse"))) {
    stop("print.bayes_ctd_array() requires measure to be equal to power, est, var, bias, mse")
  }

  if ((tolower(measure) == "var") & is.null(bayes_ctd_robj$data$var)){
    stop("Table of variance cannot be printed or plotted, because variance estimates were not requested in the trial simulation.")
  }
  if ((tolower(measure) == "bias") & is.null(bayes_ctd_robj$data$bias)){
    stop("Table of bias cannot be printed or plotted, because bias estimates were not requested in the trial simulation.")
  }
  if ((tolower(measure) == "mse") & is.null(bayes_ctd_robj$data$mse)){
    stop("Table of mse cannot be printed or plotted, because mse estimates were not requested in the trial simulation.")
  }


  if (!(toupper(tab_type) %in% c("WX|YZ", "WY|XZ", "WZ|XY", "XY|WZ", "XZ|WY", "YZ|WX", "ZX|WY", "XW|YZ", "YW|XZ", "YX|WZ", "ZW|XY", "ZX|WY", "ZY|WX"))) {
    stop("print.bayes_ctd_array() requires tab_type to be equal to WX|YZ, WY|XZ, WZ|XY, XY|WZ, XZ|WY, YZ|WX, ZX|WY, XW|YZ, YW|XZ, YX|WZ, ZW|XY, ZX|WY, or ZY|WX")
  }

  if (toupper(tab_type) == "WX|YZ") {
    if (is.null(effect_val) | is.null(rand_control_diff_val)) {
      stop("print.bayes_ctd_array() requires effect_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "WY|XZ") {
    if (is.null(a0_val) | is.null(rand_control_diff_val)) {
      stop("print.bayes_ctd_array() requires a0_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "WZ|XY") {
    if (is.null(a0_val) | is.null(effect_val)) {
      stop("print.bayes_ctd_array() requires a0_val and effect_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "XY|WZ") {
    if (is.null(subj_per_arm_val) | is.null(rand_control_diff_val)) {
      stop("print.bayes_ctd_array() requires subj_per_arm_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "XZ|WY") {
    if (is.null(subj_per_arm_val) | is.null(effect_val)) {
      stop("print.bayes_ctd_array() requires subj_per_arm_val and effect_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "YZ|WX") {
    if (is.null(subj_per_arm_val) | is.null(a0_val)) {
      stop("print.bayes_ctd_array() requires subj_per_arm_val and a0_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "ZX|WY") {
    if (is.null(subj_per_arm_val) | is.null(effect_val)) {
      stop("print.bayes_ctd_array() requires subj_per_arm_val and effect_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "XW|YZ") {
    if (is.null(effect_val) | is.null(rand_control_diff_val)) {
      stop("print.bayes_ctd_array() requires effect_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "YW|XZ") {
    if (is.null(a0_val) | is.null(rand_control_diff_val)) {
      stop("print.bayes_ctd_array() requires a0_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "YX|WZ") {
    if (is.null(subj_per_arm_val) | is.null(rand_control_diff_val)) {
      stop("print.bayes_ctd_array() requires subj_per_arm_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "ZW|XY") {
    if (is.null(a0_val) | is.null(effect_val)) {
      stop("print.bayes_ctd_array() requires a0_val and effect_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "ZX|WY") {
    if (is.null(subj_per_arm_val) | is.null(effect_val)) {
      stop("print.bayes_ctd_array() requires subj_per_arm_val and effect_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "ZY|WX") {
    if (is.null(subj_per_arm_val) | is.null(a0_val)) {
      stop("print.bayes_ctd_array() requires subj_per_arm_val and a0_val to be non-missing and equal to a value used in your simulation.")
    }
  }

  subj_per_arm_vec <- bayes_ctd_robj$subj_per_arm
  a0_val_vec <- bayes_ctd_robj$a0_vals
  effect_val_vec <- bayes_ctd_robj$effect_vals
  rand_control_diff_vec <- bayes_ctd_robj$rand_control_diff

  if (!is.null(subj_per_arm_val)) {
    if (is.na(match(subj_per_arm_val, subj_per_arm_vec))) {
      text_val <- "print.bayes_ctd_array() requires subj_per_arm_val to be equal to one of the following values: "
      for (i in 1:length(subj_per_arm_vec)) {
        text_val <- paste(text_val, subj_per_arm_vec[i], sep = "")
        if (i != length(subj_per_arm_vec)) {
          text_val <- paste(text_val, ", ", sep = "")
        }
      }
      stop(text_val)
    }
  }

  if (!is.null(a0_val)) {
    if (is.na(match(a0_val, a0_val_vec))) {
      text_val <- "print.bayes_ctd_array() requires a0_val to be equal to one of the following values: "
      for (i in 1:length(a0_val_vec)) {
        text_val <- paste(text_val, a0_val_vec[i], sep = "")
        if (i != length(a0_val_vec)) {
          text_val <- paste(text_val, ", ", sep = "")
        }
      }
      stop(text_val)
    }
  }

  if (!is.null(effect_val)) {
    if (is.na(match(effect_val, effect_val_vec))) {
      text_val <- "print.bayes_ctd_array() requires effect_val to be equal to one of the following values: "
      for (i in 1:length(effect_val_vec)) {
        text_val <- paste(text_val, effect_val_vec[i], sep = "")
        if (i != length(effect_val_vec)) {
          text_val <- paste(text_val, ", ", sep = "")
        }
      }
      stop(text_val)
    }
  }

  if (!is.null(rand_control_diff_val)) {
    if (is.na(match(rand_control_diff_val, rand_control_diff_vec))) {
      text_val <- "print.bayes_ctd_array() requires rand_control_diff_val to be equal to one of the following values: "
      for (i in 1:length(rand_control_diff_vec)) {
        text_val <- paste(text_val, rand_control_diff_vec[i], sep = "")
        if (i != length(rand_control_diff_vec)) {
          text_val <- paste(text_val, ", ", sep = "")
        }
      }
      stop(text_val)
    }
  }

}



#' Checks \code{plot_table()} Function Call for Errors.
#'
#' \code{plot_table_error_checks()} function used only internally by
#' \code{plot_table()} function to check for proper input.  Not to be called
#' directly by user.
#'
#' @param bayes_ctd_robj See help page for \code{plot_table()}.
#' @param measure See help page for \code{plot_table()}.
#' @param tab_type See help page for \code{plot_table()}.
#' @param subj_per_arm_val See help page for \code{plot_table()}.
#' @param a0_val See help page for \code{print_table()}.
#' @param effect_val See help page for \code{plot_table()}.
#' @param rand_control_diff_val See help page for \code{plot_table()}.
#' @param smooth See help page for \code{plot_table()}.
#' @param plot_out See help page for \code{plot_table()}.
#' @param span See help page for \code{plot_table()}.
#'
#' @return \code{plot_table_error_checks()} returns messages when
#' \code{plot_table()} function inputs are incorrectly specified.  Not to be called
#' directly by user.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
plot_table_error_checks <- function(bayes_ctd_robj, measure, tab_type, smooth, plot_out, subj_per_arm_val, a0_val, effect_val,
    rand_control_diff_val, span) {

    if (!(tolower(measure) %in% c("power", "est", "var", "bias", "mse"))) {
        stop("plot_table() requires measure to be equal to power, est, var, bias, mse")
    }

    if (!(toupper(tab_type) %in% c("WX|YZ", "WY|XZ", "WZ|XY", "XY|WZ", "XZ|WY", "YZ|WX", "ZX|WY", "XW|YZ", "YW|XZ", "YX|WZ", "ZW|XY", "ZX|WY", "ZY|WX"))) {
      stop("plot_table() requires tab_type to be equal to WX|YZ, WY|XZ, WZ|XY, XY|WZ, XZ|WY, YZ|WX, ZX|WY, XW|YZ, YW|XZ, YX|WZ, ZW|XY, ZX|WY, or ZY|WX")
    }

    if (toupper(tab_type) == "WX|YZ") {
        if (is.null(effect_val) | is.null(rand_control_diff_val)) {
            stop("plot_table() requires effect_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
        }
    }
    if (toupper(tab_type) == "WY|XZ") {
        if (is.null(a0_val) | is.null(rand_control_diff_val)) {
            stop("plot_table() requires a0_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
        }
    }
    if (toupper(tab_type) == "WZ|XY") {
        if (is.null(a0_val) | is.null(effect_val)) {
            stop("plot_table() requires a0_val and effect_val to be non-missing and equal to a value used in your simulation.")
        }
    }
    if (toupper(tab_type) == "XY|WZ") {
        if (is.null(subj_per_arm_val) | is.null(rand_control_diff_val)) {
            stop("plot_table() requires subj_per_arm_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
        }
    }
    if (toupper(tab_type) == "XZ|WY") {
        if (is.null(subj_per_arm_val) | is.null(effect_val)) {
            stop("plot_table() requires subj_per_arm_val and effect_val to be non-missing and equal to a value used in your simulation.")
        }
    }
    if (toupper(tab_type) == "YZ|WX") {
        if (is.null(subj_per_arm_val) | is.null(a0_val)) {
            stop("plot_table() requires subj_per_arm_val and a0_val to be non-missing and equal to a value used in your simulation.")
        }
    }
    if (toupper(tab_type) == "ZX|WY") {
      if (is.null(subj_per_arm_val) | is.null(effect_val)) {
           stop("plot_table() requires subj_per_arm_val and effect_val to be non-missing and equal to a value used in your simulation.")
        }
    }
    if (toupper(tab_type) == "XW|YZ") {
      if (is.null(effect_val) | is.null(rand_control_diff_val)) {
        stop("plot_table() requires effect_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
      }
    }
    if (toupper(tab_type) == "YW|XZ") {
      if (is.null(a0_val) | is.null(rand_control_diff_val)) {
        stop("plot_table() requires a0_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
      }
    }
    if (toupper(tab_type) == "YX|WZ") {
      if (is.null(subj_per_arm_val) | is.null(rand_control_diff_val)) {
        stop("plot_table() requires subj_per_arm_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
      }
    }
    if (toupper(tab_type) == "ZW|XY") {
      if (is.null(a0_val) | is.null(effect_val)) {
        stop("plot_table() requires a0_val and effect_val to be non-missing and equal to a value used in your simulation.")
      }
    }
    if (toupper(tab_type) == "ZX|WY") {
      if (is.null(subj_per_arm_val) | is.null(effect_val)) {
        stop("plot_table() requires subj_per_arm_val and effect_val to be non-missing and equal to a value used in your simulation.")
      }
    }
    if (toupper(tab_type) == "ZY|WX") {
      if (is.null(subj_per_arm_val) | is.null(a0_val)) {
        stop("plot_table() requires subj_per_arm_val and a0_val to be non-missing and equal to a value used in your simulation.")
      }
    }

    # Need to check that smooth is TRUE/FALSE
    if (smooth != TRUE & smooth != FALSE) {
        stop("plot_table() requires smooth to be either TRUE or FALSE")
    }

    # Need to check that plot_out is TRUE/FALSE
    if (plot_out != TRUE & plot_out != FALSE) {
        stop("plot_table() requires plot_out to be either TRUE or FALSE")
    }

    # Need to make sure span elements are all positive numbers.
    for (sp in span) {
        if (!is.numeric(sp))
            stop("plot_table() requires span to be numeric.")
        if (sp <= 0)
            stop("plot_table() requires span to be positive")
    }

    subj_per_arm_vec <- bayes_ctd_robj$subj_per_arm
    a0_val_vec <- bayes_ctd_robj$a0_vals
    effect_val_vec <- bayes_ctd_robj$effect_vals
    rand_control_diff_vec <- bayes_ctd_robj$rand_control_diff

    if (!is.null(subj_per_arm_val)) {
        if (is.na(match(subj_per_arm_val, subj_per_arm_vec))) {
            text_val <- "plot_table() requires subj_per_arm_val to be equal to one of the following values: "
            for (i in 1:length(subj_per_arm_vec)) {
                text_val <- paste(text_val, subj_per_arm_vec[i], sep = "")
                if (i != length(subj_per_arm_vec)) {
                  text_val <- paste(text_val, ", ", sep = "")
                }
            }
            stop(text_val)
        }
    }

    if (!is.null(a0_val)) {
        if (is.na(match(a0_val, a0_val_vec))) {
            text_val <- "plot_table() requires a0_val to be equal to one of the following values: "
            for (i in 1:length(a0_val_vec)) {
                text_val <- paste(text_val, a0_val_vec[i], sep = "")
                if (i != length(a0_val_vec)) {
                  text_val <- paste(text_val, ", ", sep = "")
                }
            }
            stop(text_val)
        }
    }

    if (!is.null(effect_val)) {
        if (is.na(match(effect_val, effect_val_vec))) {
            text_val <- "plot_table() requires effect_val to be equal to one of the following values: "
            for (i in 1:length(effect_val_vec)) {
                text_val <- paste(text_val, effect_val_vec[i], sep = "")
                if (i != length(effect_val_vec)) {
                  text_val <- paste(text_val, ", ", sep = "")
                }
            }
            stop(text_val)
        }
    }

    if (!is.null(rand_control_diff_val)) {
        if (is.na(match(rand_control_diff_val, rand_control_diff_vec))) {
            text_val <- "plot_table() requires rand_control_diff_val to be equal to one of the following values: "
            for (i in 1:length(rand_control_diff_vec)) {
                text_val <- paste(text_val, rand_control_diff_vec[i], sep = "")
                if (i != length(rand_control_diff_vec)) {
                  text_val <- paste(text_val, ", ", sep = "")
                }
            }
            stop(text_val)
        }
    }
}




#' Checks \code{plot.bayes_ctd_array()} Function Call for Errors.
#'
#' \code{plot_error_checks()} function used only internally by
#' \code{plot.bayes_ctd_array()} function to check for proper input.  Not to be called
#' directly by user.
#'
#' @param bayes_ctd_robj See help page for \code{plot_table()}.
#' @param measure See help page for \code{plot_table()}.
#' @param tab_type See help page for \code{plot_table()}.
#' @param subj_per_arm_val See help page for \code{plot_table()}.
#' @param a0_val See help page for \code{print_table()}.
#' @param effect_val See help page for \code{plot_table()}.
#' @param rand_control_diff_val See help page for \code{plot_table()}.
#' @param smooth See help page for \code{plot_table()}.
#' @param plot_out See help page for \code{plot_table()}.
#' @param span See help page for \code{plot_table()}.
#'
#' @return \code{plot_error_checks()} returns messages when
#' \code{plot.bayes_ctd_array()} function inputs are incorrectly specified.  Not to be called
#' directly by user.
#'
#' @examples
#' #None
#' @keywords internal
#' @noRd
plot_error_checks <- function(bayes_ctd_robj, measure, tab_type, smooth, plot_out, subj_per_arm_val, a0_val, effect_val,
                                    rand_control_diff_val, span) {

  if (!(tolower(measure) %in% c("power", "est", "var", "bias", "mse"))) {
    stop("plot.bayes_ctd_array() requires measure to be equal to power, est, var, bias, mse")
  }

  if (!(toupper(tab_type) %in% c("WX|YZ", "WY|XZ", "WZ|XY", "XY|WZ", "XZ|WY", "YZ|WX", "ZX|WY", "XW|YZ", "YW|XZ", "YX|WZ", "ZW|XY", "ZX|WY", "ZY|WX"))) {
    stop("plot.bayes_ctd_array() requires tab_type to be equal to WX|YZ, WY|XZ, WZ|XY, XY|WZ, XZ|WY, YZ|WX, ZX|WY, XW|YZ, YW|XZ, YX|WZ, ZW|XY, ZX|WY, or ZY|WX")
  }

  if (toupper(tab_type) == "WX|YZ") {
    if (is.null(effect_val) | is.null(rand_control_diff_val)) {
      stop("plot.bayes_ctd_array() requires effect_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "WY|XZ") {
    if (is.null(a0_val) | is.null(rand_control_diff_val)) {
      stop("plot.bayes_ctd_array() requires a0_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "WZ|XY") {
    if (is.null(a0_val) | is.null(effect_val)) {
      stop("plot.bayes_ctd_array() requires a0_val and effect_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "XY|WZ") {
    if (is.null(subj_per_arm_val) | is.null(rand_control_diff_val)) {
      stop("plot.bayes_ctd_array() requires subj_per_arm_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "XZ|WY") {
    if (is.null(subj_per_arm_val) | is.null(effect_val)) {
      stop("plot.bayes_ctd_array() requires subj_per_arm_val and effect_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "YZ|WX") {
    if (is.null(subj_per_arm_val) | is.null(a0_val)) {
      stop("plot.bayes_ctd_array() requires subj_per_arm_val and a0_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "ZX|WY") {
    if (is.null(subj_per_arm_val) | is.null(effect_val)) {
      stop("plot.bayes_ctd_array() requires subj_per_arm_val and effect_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "XW|YZ") {
    if (is.null(effect_val) | is.null(rand_control_diff_val)) {
      stop("plot.bayes_ctd_array() requires effect_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "YW|XZ") {
    if (is.null(a0_val) | is.null(rand_control_diff_val)) {
      stop("plot.bayes_ctd_array() requires a0_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "YX|WZ") {
    if (is.null(subj_per_arm_val) | is.null(rand_control_diff_val)) {
      stop("plot.bayes_ctd_array() requires subj_per_arm_val and rand_control_diff_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "ZW|XY") {
    if (is.null(a0_val) | is.null(effect_val)) {
      stop("plot.bayes_ctd_array() requires a0_val and effect_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "ZX|WY") {
    if (is.null(subj_per_arm_val) | is.null(effect_val)) {
      stop("plot.bayes_ctd_array() requires subj_per_arm_val and effect_val to be non-missing and equal to a value used in your simulation.")
    }
  }
  if (toupper(tab_type) == "ZY|WX") {
    if (is.null(subj_per_arm_val) | is.null(a0_val)) {
      stop("plot.bayes_ctd_array() requires subj_per_arm_val and a0_val to be non-missing and equal to a value used in your simulation.")
    }
  }

  # Need to check that smooth is TRUE/FALSE
  if (smooth != TRUE & smooth != FALSE) {
    stop("plot.bayes_ctd_array() requires smooth to be either TRUE or FALSE")
  }

  # Need to check that plot_out is TRUE/FALSE
  if (plot_out != TRUE & plot_out != FALSE) {
    stop("plot.bayes_ctd_array() requires plot_out to be either TRUE or FALSE")
  }

  # Need to make sure span elements are all positive numbers.
  for (sp in span) {
    if (!is.numeric(sp))
      stop("plot.bayes_ctd_array() requires span to be numeric.")
    if (sp <= 0)
      stop("plot.bayes_ctd_array() requires span to be positive")
  }

  subj_per_arm_vec <- bayes_ctd_robj$subj_per_arm
  a0_val_vec <- bayes_ctd_robj$a0_vals
  effect_val_vec <- bayes_ctd_robj$effect_vals
  rand_control_diff_vec <- bayes_ctd_robj$rand_control_diff

  if (!is.null(subj_per_arm_val)) {
    if (is.na(match(subj_per_arm_val, subj_per_arm_vec))) {
      text_val <- "plot.bayes_ctd_array() requires subj_per_arm_val to be equal to one of the following values: "
      for (i in 1:length(subj_per_arm_vec)) {
        text_val <- paste(text_val, subj_per_arm_vec[i], sep = "")
        if (i != length(subj_per_arm_vec)) {
          text_val <- paste(text_val, ", ", sep = "")
        }
      }
      stop(text_val)
    }
  }

  if (!is.null(a0_val)) {
    if (is.na(match(a0_val, a0_val_vec))) {
      text_val <- "plot.bayes_ctd_array() requires a0_val to be equal to one of the following values: "
      for (i in 1:length(a0_val_vec)) {
        text_val <- paste(text_val, a0_val_vec[i], sep = "")
        if (i != length(a0_val_vec)) {
          text_val <- paste(text_val, ", ", sep = "")
        }
      }
      stop(text_val)
    }
  }

  if (!is.null(effect_val)) {
    if (is.na(match(effect_val, effect_val_vec))) {
      text_val <- "plot.bayes_ctd_array() requires effect_val to be equal to one of the following values: "
      for (i in 1:length(effect_val_vec)) {
        text_val <- paste(text_val, effect_val_vec[i], sep = "")
        if (i != length(effect_val_vec)) {
          text_val <- paste(text_val, ", ", sep = "")
        }
      }
      stop(text_val)
    }
  }

  if (!is.null(rand_control_diff_val)) {
    if (is.na(match(rand_control_diff_val, rand_control_diff_vec))) {
      text_val <- "plot.bayes_ctd_array() requires rand_control_diff_val to be equal to one of the following values: "
      for (i in 1:length(rand_control_diff_vec)) {
        text_val <- paste(text_val, rand_control_diff_vec[i], sep = "")
        if (i != length(rand_control_diff_vec)) {
          text_val <- paste(text_val, ", ", sep = "")
        }
      }
      stop(text_val)
    }
  }
}
