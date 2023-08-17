
function (model, series, series_length = 25, conditions = NULL,
    exclude_random = TRUE, exclude_terms = NULL, split = NULL,
    sep = "\\.", time_series, transform = NULL, ci_z = 1.96,
    .comparison = NULL)
{
    if (!missing(time_series)) {
        warning("The time_series argument has been deprecated and will be removed in the future. Please use `series` instead.")
        series_q = dplyr::enquo(time_series)
    }
    else {
        time_series = NULL
        series_q <- dplyr::enquo(series)
    }
    .comparison_q <- dplyr::enquo(.comparison)
    if (!is.null(model$dinfo) && (exclude_random || !is.null(exclude_terms))) {
        stop("Excluding random effects and/or terms is not currently supported with discretised models (fitted with discrete = TRUE). Please, set 'exclude_random' to FALSE and/or 'exclude_terms' to NULL.")
    }
    series_name <- rlang::quo_name(series_q)
    outcome_q <- model$formula[[2]]
    cond_terms <- NULL
    if (!is.null(conditions)) {
        for (cond_i in 1:length(conditions)) {
            cond_term <- as.character(rlang::quo_get_expr(conditions[[cond_i]])[2])
            cond_terms <- c(cond_terms, cond_term)
        }
    }
    if (!(rlang::quo_is_null(.comparison_q))) {
        cond_terms <- c(cond_terms, rlang::as_name(.comparison_q))
    }
    random_effects <- list()
    random_effects_terms <- NULL
    re_term <- NULL
    if (exclude_random == TRUE) {
        for (i in 1:length(model[["smooth"]])) {
            smooth_term <- model[["smooth"]][[i]][["term"]][[1]]
            smooth_class <- attr(model$smooth[[i]], "class")[1]
            if (smooth_class == "fs.interaction" && !(smooth_term %in%
                cond_terms)) {
                random_effects <- c(random_effects, list(model$smooth[[i]]$label))
                random_effects_terms <- c(random_effects_terms,
                  model$smooth[[i]]$fterm)
            }
            if (smooth_class == "random.effect" && !(smooth_term %in%
                cond_terms)) {
                random_effects <- c(random_effects, list(model[["smooth"]][[i]]$label))
                random_effects_terms <- c(random_effects_terms,
                  model$smooth[[i]]$term)
                re_term <- c(re_term, model$smooth[[i]]$term)
            }
        }
    }
    if (exclude_random) {
        if (rlang::is_empty(random_effects)) {
            exclude_random_effects <- as.null()
        }
        else {
            exclude_random_effects <- random_effects
        }
    }
    else {
        exclude_random_effects <- as.null()
    }
    exclude_smooths <- as.null()
    excluded_terms <- as.null()
    for (smooth in 1:length(model[["smooth"]])) {
        smooth_class <- attr(model$smooth[[smooth]], "class")[1]
        smooth_term <- model[["smooth"]][[smooth]][["term"]][[1]]
        if (smooth_class == "random.effect") {
            exclude_smooths <- exclude_smooths
            excluded_terms <- excluded_terms
        }
        else if (smooth_term != series_name && !(smooth_term %in%
            cond_terms)) {
            excluded_terms <- c(excluded_terms, smooth_term)
            smooth_label <- model[["smooth"]][[smooth]][["label"]]
            exclude_smooths <- c(exclude_smooths, smooth_label)
        }
        else if (smooth_class == "tensor.smooth") {
            smooth_term <- model[["smooth"]][[smooth]][["term"]][[2]]
            excluded_terms <- c(excluded_terms, smooth_term)
            smooth_label <- model[["smooth"]][[smooth]][["label"]]
            exclude_smooths <- c(exclude_smooths, smooth_label)
        }
    }
    excluded <- as.null()
    if (!is.null(exclude_terms)) {
        for (term in 1:length(exclude_terms)) {
            for (label in 1:length(model[["smooth"]])) {
                smooth_label <- model[["smooth"]][[label]][["label"]]
                if (smooth_label == exclude_terms[term]) {
                  smooth_term <- model[["smooth"]][[label]][["term"]]
                  if (!(smooth_term %in% cond_terms)) {
                    if (length(smooth_term) > 1) {
                      smooth_term_2 <- model[["smooth"]][[label]][["term"]][[2]]
                      excluded <- c(excluded, smooth_term_2)
                    }
                  }
                }
            }
        }
    }
    exclude_these <- c(exclude_random_effects, exclude_smooths,
        exclude_terms)
    var_list <- list()
    for (var in 1:length(model[["var.summary"]])) {
        var_class <- class(model[["var.summary"]][[var]])
        if (var_class == "numeric") {
            if (names(model[["var.summary"]][var]) %in% c(random_effects_terms,
                excluded_terms, excluded)) {
                var_values <- model[["var.summary"]][[var]][[1]]
            }
            else {
                var_values <- seq(model[["var.summary"]][[var]][[1]],
                  model[["var.summary"]][[var]][[3]], length.out = series_length)
            }
        }
        else if (var_class == "factor") {
            if (names(model[["var.summary"]][var]) %in% c(random_effects_terms,
                excluded_terms, excluded)) {
                var_values <- model[["var.summary"]][[var]][[1]]
            }
            else {
                var_values <- levels(model[["var.summary"]][[var]])
            }
        }
        var_values_list <- list(var_values)
        names(var_values_list)[[1]] <- names(model[["var.summary"]][var])
        var_list <- c(var_list, var_values_list)
    }
    fitted_df <- expand.grid(var_list)
    if ("(AR.start)" %in% colnames(fitted_df)) {
        fitted_df$`(AR.start)` <- NULL
    }
    predicted <- stats::predict(model, fitted_df, se.fit = TRUE,
        exclude = exclude_these)
    predicted_tbl <- cbind(fitted_df, predicted) %>% dplyr::mutate(CI_upper = fit +
        ci_z * se.fit, CI_lower = fit - ci_z * se.fit)
    if (!is.null(transform)) {
        trans_fun <- transform
        predicted_tbl$CI_upper <- trans_fun(predicted_tbl$CI_upper)
        predicted_tbl$CI_lower <- trans_fun(predicted_tbl$CI_lower)
        predicted_tbl$fit <- trans_fun(predicted_tbl$fit)
    }
    predicted_tbl <- predicted_tbl %>% dplyr::rename(`:=`(!!rlang::quo_name(outcome_q),
        fit), SE = se.fit)
    if (!is.null(exclude_random_effects)) {
        for (term in 1:length(random_effects_terms)) {
            if (random_effects_terms[term] %in% re_term) {
                random_effects_terms <- random_effects_terms[-term]
            }
        }
        predicted_tbl <- predicted_tbl %>% dplyr::select(-dplyr::one_of(random_effects_terms)) %>%
            unique()
    }
    if (!is.null(exclude_smooths)) {
        for (term in 1:length(excluded_terms)) {
            if (excluded_terms[term] %in% re_term) {
                excluded_terms <- excluded_terms[-term]
            }
        }
        predicted_tbl <- predicted_tbl %>% dplyr::select(-dplyr::one_of(excluded_terms)) %>%
            unique()
    }
    if (!is.null(excluded)) {
        predicted_tbl <- predicted_tbl %>% dplyr::select(-dplyr::one_of(excluded)) %>%
            unique()
    }
    if (!is.null(split)) {
        for (i in 1:length(split)) {
            predicted_tbl <- tidyr::separate(data = predicted_tbl,
                col = names(split)[i], into = split[[i]], sep = sep)
        }
    }
    if (!is.null(conditions)) {
        predicted_tbl <- predicted_tbl %>% dplyr::filter(!!!conditions)
    }
    if (ncol(predicted_tbl) > 5) {
        predicted_tbl <- tidyr::unite(predicted_tbl, ".idx",
            c(-CI_lower, -CI_upper, -SE, -!!rlang::quo_name(outcome_q),
                -!!rlang::quo_name(series_q)), remove = FALSE) %>%
            dplyr::mutate(.idx = as.numeric(as.factor(.idx)))
    }
    return(predicted_tbl)
}

