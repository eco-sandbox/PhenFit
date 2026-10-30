# PhenFit: An interactive model of selection on reaction norms in harsh versus nice environments

## The model

- Gaussian stabilizing selection with unit width: `w(z, theta) = exp(-(z - theta)^2 / 2)`.
- Two environments: *harsh* (phenotypic optimum 0, maximum fitness `ratio` < 1) and *nice* (optimum `delta`, maximum fitness 1). The nice environment occurs with frequency `p_nice`.
- **Genotype = intercept** `a` **and plasticity slope** `b` **in [0, 1].** If its cue says "nice" it expresses `z = a + b * delta`, otherwise `z = a`; `b = 1` shifts the phenotype by the full optimum difference, `b = 0` is canalized. The cue is correct with probability `rho` (cue reliability).
- **Cost of plasticity:** fitness is multiplied by `1 - cost * b^2`, so it is paid only through the slope (canalized genotypes pay nothing) and the evolved slope `b*` shrinks as `cost` rises. Fitness in environment e: `Wmax_e * (1 - cost * b^2) * E_cue[exp(-(z - theta_e)^2 / 2)]`.
- **Evolved (plastic) genotype:** the `(a*, b*)` found by grid search (`a` in `[-1, delta + 1]` step 0.005, `b` in `[0, 1]` step 0.01) that maximizes the environment-weighted fitness.
  - *hard selection:* weights harsh `(1 - p) * ratio`, nice `p`
  - *soft selection:* weights harsh `1 - p`, nice `p`
- **Canalized genotype:** the best intercept with `b = 0` (no cost).
- Outputs: reaction norms, fitness-vs-slope curve (peak = `b*`), fitness of each genotype in each environment and on average, relative and absolute gain of plasticity, `b*`, `a*`, the cost paid, the canalized `z*`, and a plain-language verdict.


## Files

- `model.R` - pure, documented functions (`run_model()`, `canalized_optimum()`, plots). Run `Rscript model.R` for a worked example (prints b* for several costs, tables, and writes `example.png`).
- `app.R` - Shiny app with sliders for all parameters (sources `model.R`).

## Run

```r
install.packages(c("shiny", "ggplot2"))
setwd("PhenFit") #adjust if needed
shiny::runApp()            # opens the app
source("model.R"); m <- run_model(p_nice = 0.7, ratio = 0.4, delta = 2, rho = 0.9, cost = 0.05, regime = "hard"); results_table(m)
```
