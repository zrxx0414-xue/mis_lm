# DRM model selection.
# Linear regression simulation and model selection.

library(drmdel)

drm_model_selection <- function(baseline, prediction_1, prediction_2, basis_name,
                                complexity = c(M1 = 1L, M2 = 2L), alpha = 0.05, maxit = 2000L,
                                kernel_bw = "adaptive", dim_basis = "adaptive") {
    stage <- "basis_construction"
    basis_ready <- FALSE
    convergence <- NA_integer_
    optimizer_message <- NA_character_
    d <- NA_integer_
    err <- NA_character_
    error_call <- NA_character_
    warnings <- character()
    decision <- NA_character_
    reason <- "Invalid DRM fit"
    rejected <- NA
    equality_wald <- equality_p <- wald_1 <- wald_2 <- NA_real_
    basis <- if (is.function(basis_name)) "custom" else basis_name

    tryCatch(withCallingHandlers({
        data <- list(baseline, prediction_1, prediction_2)
        basis_func <- if (identical(basis_name, "adaptive")) {
            fit_adaptive_basis(data = data, kernel_bw = kernel_bw, dim_basis = dim_basis)$fun
        } else {
            basis_name
        }
        basis_ready <- TRUE

        stage <- "drm_setup"
        fitter <- drmdel::drmdel
        environment(fitter) <- list2env(list(optim = function(...) {
            stage <<- "drm_optimization"
            opt <- stats::optim(...)
            convergence <<- opt$convergence
            if (!is.null(opt$message)) optimizer_message <<- opt$message
            if (opt$convergence != 0L)
                warning("DRM optimizer did not converge; optim code = ", opt$convergence, call. = FALSE)
            stage <<- "drm_postprocessing"
            opt
        }), parent = environment(fitter))
        pooled <- unlist(data, use.names = FALSE)
        n <- length(pooled)
        fit <- fitter(x = pooled, n_samples = lengths(data), basis_func = basis_func,
                      control = list(maxit = maxit))
        d <- fit$drm_info$d

        stage <- "covariance"
        estimate <- fit$mele
        covariance <- as.matrix(drmdel::meleCov(fit))
        i1 <- 2:(d + 1L)
        i2 <- (d + 3L):(2L * d + 2L)
        wald <- function(b, S) n * drop(crossprod(b, solve(S, b)))

        stage <- "wald_equality"
        diff_cov <- covariance[i1, i1, drop = FALSE] + covariance[i2, i2, drop = FALSE] -
            covariance[i1, i2, drop = FALSE] - covariance[i2, i1, drop = FALSE]
        equality_wald <- wald(estimate[i1] - estimate[i2], diff_cov)
        equality_p <- pchisq(equality_wald, df = d, lower.tail = FALSE)
        stage <- "wald_M1"
        wald_1 <- wald(estimate[i1], covariance[i1, i1, drop = FALSE])
        stage <- "wald_M2"
        wald_2 <- wald(estimate[i2], covariance[i2, i2, drop = FALSE])

        
        stage <- "model_selection"
        rejected <- equality_p < alpha
        if (!rejected) {
            ties <- names(complexity)[complexity == min(complexity)]
            decision <- if (length(ties) == 1L) ties else "Tie"
            reason <- "Failed to reject equality; complexity tie-break"
        } else if (wald_1 != wald_2) {
            decision <- if (wald_1 < wald_2) "M1" else "M2"
            reason <- paste("Rejected equality;", decision, "has smaller individual Wald")
        } else {
            decision <- "Tie"
            reason <- "Rejected equality; individual Wald statistics are equal"
        }
    },
    warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    },
    error = function(e) {
        calls <- vapply(sys.calls(), function(x) paste(deparse(x[[1]]), collapse = ""), character(1))
        if ("adaptive_bw" %in% calls) stage <<- "adaptive_bandwidth"
        else if ("adaptive_dim" %in% calls) stage <<- "adaptive_dimension"
        else if ("fit_adaptive_basis" %in% calls && "eigen" %in% calls) stage <<- "basis_eigensystem"
        error_call <<- if (is.null(conditionCall(e))) stage
        else substr(paste(deparse(conditionCall(e)), collapse = " "), 1L, 400L)
    }), error = function(e) err <<- conditionMessage(e))

    valid <- is.na(err)
    error_type <- NA_character_
    if (!valid) {
        error_type <- if (!is.na(convergence) && convergence != 0) "nonconvergence" else "other"
        patterns <- c(nonfinite_values = "infinite|missing values|non-finite|NA/NaN/Inf",
                      singular_matrix = "singular|reciprocal condition",
                      invalid_dimension = "dimension of the basis")
        hit <- names(patterns)[vapply(patterns, grepl, logical(1), x = err)]
        if (length(hit)) error_type <- hit[length(hit)]
    }
    data.frame(decision = decision, reason = reason, rejected = rejected,
               equality_wald = equality_wald, equality_p = equality_p,
               wald_1 = wald_1, wald_2 = wald_2, basis = basis, basis_dimension = d,
               valid = valid, failure_stage = if (valid) NA_character_ else stage,
               basis_ready = basis_ready, drm_convergence = convergence,
               optimizer_message = optimizer_message, error_message = err,
               error_type = error_type, error_call = error_call,
               warning_message = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
               stringsAsFactors = FALSE)
}
safe_drm_model_selection <- drm_selection_rule <- safe_drm_selection_rule <- drm_model_selection
