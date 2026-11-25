source("drm_selection.R")
source("data_generating.R")
source("benchmark_aic_bic.R")
source("adaptive_basis.R")
library(dplyr)

run_linear_simulation <- function(true_formula, candidate_formulas, snr = 1,
    n_train = 200L, n_test = 200L, repetitions = 500L, basis = "adaptive",
    kernel_bw = "adaptive", dim_basis = "adaptive", predictor_distribution = "Normal",
    error_distribution = "Normal", test_seed = 4100L, training_seed = 42L, alpha = 0.05, maxit = 2000L,
    save_drm_data = FALSE) {
    if (!inherits(true_formula, "formula") || length(true_formula) != 2L)
        stop("Write true_formula as a one-sided formula, such as ~ 1 + x1.")
    if (!is.list(candidate_formulas) || length(candidate_formulas) != 2L ||
        !setequal(names(candidate_formulas), c("M1", "M2")) ||
        !all(vapply(candidate_formulas, function(f) {
            inherits(f, "formula") && length(f) == 3L && identical(f[[2L]], quote(y))
        }, logical(1))))
        stop("Supply two formulas with response y in a list named M1 and M2.")
    candidate_formulas <- candidate_formulas[c("M1", "M2")]
    if (!is.numeric(snr) || length(snr) != 1L || !is.finite(snr) || snr <= 0)
        stop("snr must be a positive finite number.")

    predictors <- unique(c(all.vars(true_formula[[2L]]),
        all.vars(candidate_formulas$M1[[3L]]), all.vars(candidate_formulas$M2[[3L]])))
    predictors <- setdiff(predictors, ".")
    if (!length(predictors) || "y" %in% predictors)
        stop("Use at least one predictor, and reserve y for the response.")

    start <- proc.time()[["elapsed"]]
    set.seed(test_seed)
    x_test <- generate_predictors(n_test, predictors, predictor_distribution)
    set.seed(training_seed)
    results <- parameters <- vector("list", repetitions)
    if (save_drm_data) drm_data <- vector("list", repetitions)

    for (i in seq_len(repetitions)) {
        data <- generate_experiment(generate_predictors(n_train, predictors, predictor_distribution), x_test,
            true_formula, snr, error_distribution)
        models <- lapply(candidate_formulas, lm, data = data$train)
        complexity <- vapply(models, function(model) model$rank, integer(1))
        drm_input <- list(
            baseline = data$test$y,
            prediction_1 = as.numeric(predict(models$M1, newdata = data$test)),
            prediction_2 = as.numeric(predict(models$M2, newdata = data$test))
        )
        parameters[[i]] <- mutate(data$parameters,
            rmse_M1 = sqrt(mean((drm_input$prediction_1 - drm_input$baseline)^2)),
            rmse_M2 = sqrt(mean((drm_input$prediction_2 - drm_input$baseline)^2)))
        if (save_drm_data) drm_data[[i]] <- drm_input

        drm_start <- proc.time()[["elapsed"]]
        drm <- drm_model_selection(drm_input$baseline, drm_input$prediction_1, drm_input$prediction_2,
            basis, complexity = complexity, alpha = alpha, maxit = maxit,
            kernel_bw = kernel_bw, dim_basis = dim_basis)
        drm$runtime_seconds <- proc.time()[["elapsed"]] - drm_start
        drm$method <- "DRM"

        benchmarks <- benchmark_aic_bic(models)
        reference <- data.frame(method = names(benchmarks), decision = unname(benchmarks), valid = TRUE)
        results[[i]] <- mutate(bind_rows(drm, reference), replicate = i)
    }

    parameters <- bind_rows(parameters) |> mutate(replicate = row_number())
    raw <- bind_rows(results) |> left_join(parameters, by = "replicate") |>
        mutate(snr = snr, n_train = n_train, n_test = n_test, error_sd = sqrt(error_variance))
    drm_raw <- filter(raw, method == "DRM")

    summary <- raw |> group_by(method) |>
        summarise(
            valid_runs = sum(valid), invalid_runs = sum(!valid),
            M1_rate = mean(decision[valid] == "M1"),
            M2_rate = mean(decision[valid] == "M2"),
            tie_rate = mean(decision[valid] == "Tie"),
            .groups = "drop"
        ) |>
        mutate(across(c(M1_rate, M2_rate, tie_rate), ~ replace(.x, is.nan(.x), NA_real_)))
    error_summary <- drm_raw |> filter(!valid) |>
        count(failure_stage, error_type, error_message, name = "runs")

    settings <- list(true_formula = true_formula, candidate_formulas = candidate_formulas, snr = snr,
        n_train = n_train, n_test = n_test, repetitions = repetitions, basis = basis,
        kernel_bw = kernel_bw, dim_basis = dim_basis, predictor_distribution = predictor_distribution,
        error_distribution = error_distribution, test_seed = test_seed, training_seed = training_seed,
        alpha = alpha, maxit = maxit, save_drm_data = save_drm_data)
                  
    result <- list(settings = settings, parameters = parameters, raw = raw, summary = summary,
        error_summary = error_summary,
        runtime = list(total_seconds = proc.time()[["elapsed"]] - start,
                       drm_total_seconds = sum(drm_raw$runtime_seconds),
                       drm_mean_seconds = mean(drm_raw$runtime_seconds)))
    if (save_drm_data) result$drm_data <- drm_data
    result
}
