source("run_simulation.R")

result <- run_linear_simulation(
    true_formula = ~ 1 + x1 + x2 + 0.5 * x1 * x2,
    candidate_formulas = list(
        M1 = y ~ x1 + x2,
        M2 = y ~ x1 + x1:x2
    ),
    snr = 3,
    repetitions = 10L,
    dim_basis = 1L,
    save_drm_data = TRUE
)

print(result$summary)
print(result$error_summary)
print(result$runtime)
