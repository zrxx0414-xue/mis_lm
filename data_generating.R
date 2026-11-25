standardised_random <- function(n, distribution = "Normal") {
    values <- switch(distribution,
        "Normal" = rnorm(n), "Uniform" = runif(n, -sqrt(3), sqrt(3)),
        "Heavy-tailed t(3)" = rt(n, df = 3) / sqrt(3),
        "Bimodal" = (sample(c(-1, 1), n, replace = TRUE) + rnorm(n, sd = 0.5)) / sqrt(1.25),
        "Right-skewed" = rexp(n) - 1, stop("Unknown distribution: ", distribution))
    as.numeric(scale(values))
}

#generate x
generate_predictors <- function(n, predictors, distribution = "Normal") {
    values <- lapply(predictors, function(name) standardised_random(n, distribution))
    as.data.frame(setNames(values, predictors), check.names = FALSE)
}

generate_experiment <- function(train, test, true_formula, snr, error_distribution) {
    mu <- eval(true_formula[[2L]], train, environment(true_formula))
    test_signal <- eval(true_formula[[2L]], test, environment(true_formula))
    if (!is.numeric(mu) || !is.null(dim(mu)) || length(mu) != nrow(train) || any(!is.finite(mu)) ||
        !is.numeric(test_signal) || !is.null(dim(test_signal)) ||
        length(test_signal) != nrow(test) || any(!is.finite(test_signal)))
        stop("true_formula must yield one finite numeric signal value per row.")

    signal_variance <- var(mu)
    error_variance <- signal_variance / snr
    if (!is.finite(error_variance) || error_variance <= 0)
        stop("The training signal and snr must yield a positive finite error variance.")

    train$y <- mu + sqrt(error_variance) * standardised_random(nrow(train), error_distribution)
    test$y <- test_signal
    parameters <- data.frame(error_variance = error_variance, achieved_snr = signal_variance / error_variance)
    list(train = train, test = test, parameters = parameters)
}
