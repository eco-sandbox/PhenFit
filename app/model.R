# model.R -- Illustrative two-environment model of phenotypic plasticity
# ---------------------------------------------------------------------------
# R port of the interactive model PhenFit
#
# Setup
#  * Trait z; Gaussian stabilizing selection with width omega_harsh / omega_nice
#    (default 1 in both; v4, 2026-10-29: separate widths).
#  * Optimum theta is 0 in the harsh environment and `delta` in the nice one.
#  * Maximum fitness: 1 in the nice environment, `ratio` (<= 1) in the harsh one.
#  * A genotype is an intercept `a` and a plasticity slope `b` in [0, 1].
#    If its cue says "nice" it expresses z = a + b * delta, otherwise z = a.
#    b = 1 shifts the phenotype by the full optimum difference; b = 0 is canalized.
#  * The cue is correct with probability `rho` (0.5 = uninformative).
#  * Cost of plasticity: fitness is multiplied by (1 - cost * b^2). The cost is
#    paid only through the slope, so b = 0 (canalized) is free, and the evolved
#    slope b* shrinks as the cost rises.
#  * Fitness of (a, b) in environment e:
#        W_e = Wmax_e * (1 - cost * b^2) * E_cue[ exp(-(z - theta_e)^2 / (2 omega_e^2)) ]
#    with E_cue over the sensed state (correct w.p. rho).
#  * Selection picks the (a, b) maximizing the weighted sum of the two
#    environments' fitness (grid search: a in [-1, delta + 1] step 0.005,
#    b in [0, 1] step 0.01; first maximum wins). Weights:
#      "hard": frequency x maximum fitness  (harsh: (1-f)*ratio, nice: f)
#      "soft": frequency                    (harsh: 1-f,        nice: f)
#    where f is the frequency of the nice environment (argument f_nice).
#  * "Canalized genotype" = the best intercept with b = 0 (no cost).
#  * "Plastic genotype"   = the evolved (a*, b*); it can never do worse than the
#    canalized one on the selection objective, and equals it when b* = 0.

#' Gaussian stabilizing-selection fitness kernel (unit width, max 1)
#' @param z phenotype
#' @param theta optimum
gauss_fitness <- function(z, theta, omega = 1) exp(-((z - theta)^2) / (2 * omega^2))

#' Cue-averaged fitness shape (before Wmax and cost) in each environment
#' @return list(S_harsh, S_nice) for intercept vector `a` and slope `b`
fitness_shape <- function(a, b, delta, rho, om_h = 1, om_n = 1) {
  zn <- a + b * delta   # phenotype expressed when the cue says "nice"
  list(S_harsh = rho * gauss_fitness(a, 0, om_h)     + (1 - rho) * gauss_fitness(zn, 0, om_h),
       S_nice  = rho * gauss_fitness(zn, delta, om_n) + (1 - rho) * gauss_fitness(a, delta, om_n))
}

#' Evolved genotype by grid search over (a, b)
#' @param f_nice frequency f of the nice environment
#' @param ratio  maximum fitness in harsh / nice environment
#' @param delta  distance between optima (in selection widths)
#' @param rho    cue reliability
#' @param cost   cost coefficient: fitness factor 1 - cost * b^2
#' @param regime "hard" or "soft"
#' @return list(a_star, b_star, z_star (canalized intercept, b = 0),
#'         b_grid, profile (best objective for each b))
solve_slope <- function(f_nice, ratio, delta, rho, cost, regime = c("hard", "soft"),
                        a_step = 0.005, b_step = 0.01, om_h = 1, om_n = 1) {
  regime <- match.arg(regime)
  w_harsh <- if (regime == "hard") (1 - f_nice) * ratio else 1 - f_nice
  w_nice  <- f_nice
  a <- -1 + (0:round((delta + 2) / a_step)) * a_step
  b_grid <- (0:round(1 / b_step)) * b_step
  profile <- numeric(length(b_grid)); a_at <- numeric(length(b_grid))
  for (j in seq_along(b_grid)) {
    b <- b_grid[j]
    S <- fitness_shape(a, b, delta, rho, om_h, om_n)
    v <- (w_harsh * S$S_harsh + w_nice * S$S_nice) * (1 - cost * b^2)
    k <- which.max(v)            # first maximum, as in the JS strict '>' loop
    profile[j] <- v[k]; a_at[j] <- a[k]
  }
  jb <- which.max(profile)       # first maximum over b
  list(a_star = a_at[jb], b_star = b_grid[jb], z_star = a_at[1],
       b_grid = b_grid, profile = profile)
}

#' Best fixed (canalized) phenotype, kept for convenience (b = 0)
canalized_optimum <- function(f_nice, ratio, delta, regime = c("hard", "soft")) {
  solve_slope(f_nice, ratio, delta, rho = 1, cost = 0, regime = regime)$z_star
}

#' Run the  model
#' @param f_nice frequency f of the nice environment (0-1)
#' @param ratio  harsh/nice maximum fitness (0-1]
#' @param delta  optimum difference in widths (>0)
#' @param rho    cue reliability (0.5 = uninformative, 1 = perfect)
#' @param cost   cost coefficient c: fitness is multiplied by 1 - c * b^2
#' @param regime "hard" or "soft" selection
#' @return named list: inputs, z_star (canalized), a_star, b_star, cost_paid,
#'   b_grid, profile, fitness of plastic (pH, pN, pMean) and canalized
#'   (cH, cN, cMean) genotypes, relative / absolute gains, verdict.
run_model <- function(f_nice = 0.7, ratio = 0.4, delta = 2, rho = 0.9,
                          cost = 0.05, regime = c("hard", "soft"),
                          om_h = 1, om_n = 1) {
  regime <- match.arg(regime)
  sol <- solve_slope(f_nice, ratio, delta, rho, cost, regime, om_h = om_h, om_n = om_n)
  a_star <- sol$a_star; b_star <- sol$b_star; z_star <- sol$z_star

  pen <- 1 - cost * b_star^2
  S <- fitness_shape(a_star, b_star, delta, rho, om_h, om_n)
  pH <- ratio * pen * S$S_harsh
  pN <- 1 * pen * S$S_nice
  cH <- ratio * gauss_fitness(z_star, 0, om_h)
  cN <- 1 * gauss_fitness(z_star, delta, om_n)
  pMean <- f_nice * pN + (1 - f_nice) * pH
  cMean <- f_nice * cN + (1 - f_nice) * cH
  relH <- (pH - cH) / cH
  relN <- (pN - cN) / cN
  relMean <- (pMean - cMean) / cMean

  frac <- if (delta > 0) z_star / delta else 0
  frac <- min(max(frac, 0), 1)
  larger <- if (relH > relN) "harsh" else "nice"
  verdict <- sprintf(
    paste0("Under these settings the best fixed phenotype sits %.0f%% of the way ",
           "from the harsh optimum to the nice optimum. The plastic genotype gains ",
           "%.0f%% in the harsh environment and %.0f%% in the nice environment ",
           "relative to the canalized genotype, so the relative advantage of ",
           "plasticity is larger in the %s environment. The evolved plasticity ",
           "slope is b* = %.2f (1 means the phenotype tracks the full optimum shift ",
           "when the cue says so), at a fitness cost of %.1f%%."),
    frac * 100, relH * 100, relN * 100, larger,
    b_star, cost * b_star^2 * 100)

  list(
    inputs = list(f_nice = f_nice, ratio = ratio, delta = delta, rho = rho,
                  cost = cost, regime = regime, om_h = om_h, om_n = om_n),
    z_star = z_star, a_star = a_star, b_star = b_star,
    cost_paid = cost * b_star^2,
    b_grid = sol$b_grid, profile = sol$profile,
    pH = pH, pN = pN, pMean = pMean,
    cH = cH, cN = cN, cMean = cMean,
    relH = relH, relN = relN, relMean = relMean,
    absH = pH - cH, absN = pN - cN, absMean = pMean - cMean,
    larger_relative_gain = larger,
    verdict = verdict
  )
}

#' Results as a data.frame
results_table <- function(m) {
  data.frame(
    Environment = c("Harsh", "Nice", "Population mean"),
    Plastic = c(m$pH, m$pN, m$pMean),
    Canalized = c(m$cH, m$cN, m$cMean),
    Relative_gain = c(m$relH, m$relN, m$relMean),
    Absolute_gain = c(m$absH, m$absN, m$absMean),
    row.names = NULL
  )
}

#' Evolved-plasticity summary (slope, intercept, cost) as a data.frame
plasticity_table <- function(m) {
  data.frame(
    Quantity = c("Evolved slope b* (1 = tracks the optimum shift)",
                 "Fitness cost paid, c\u00b7b*\u00b2",
                 "Intercept a* (phenotype when the cue says harsh)",
                 "Canalized phenotype z*"),
    Value = c(m$b_star, m$cost_paid, m$a_star, m$z_star))
}

# mean phenotype of the evolved plastic genotype in each environment
plastic_means <- function(m) {
  D <- m$inputs$delta; q <- m$inputs$rho
  c(harsh = m$a_star + m$b_star * D * (1 - q),
    nice  = m$a_star + m$b_star * D * q)
}

#' Reaction-norm figure 
#' Dashed line joins the two optima; viridis blue (#3B528B) = evolved plastic
#' genotype's mean phenotype; viridis green (#5DC863) = canalized phenotype z*.
plot_reaction_norms <- function(m) {
  D <- m$inputs$delta; pm <- plastic_means(m)
  plot(NA, xlim = c(0.8, 2.2), ylim = c(-1, D + 1), xaxt = "n",
       xlab = "", ylab = "phenotype z", bty = "l",
       main = sprintf("Reaction norms (evolved slope b* = %.2f)", m$b_star))
  axis(1, at = c(1, 2), labels = c("Harsh (low fitness)", "Nice (high fitness)"))
  segments(1, 0, 2, D, lty = 2, col = "grey50")
  segments(1, m$z_star, 2, m$z_star, col = "#5DC863", lwd = 3)
  segments(1, pm[["harsh"]], 2, pm[["nice"]], col = "#3B528B", lwd = 3)
  points(c(1, 2), c(0, D), pch = 19, cex = 1.6, col = "grey30")
  text(c(1, 2), c(0, D), "optimum", pos = c(2, 4), col = "grey30", xpd = NA)
  legend("topleft", bty = "n", lwd = 3, col = c("#3B528B", "#5DC863"),
         legend = c("Plastic genotype (mean phenotype)", "Canalized genotype"))
}

#' Fitness-vs-slope curve : best weighted fitness for each slope b
plot_slope_profile <- function(m) {
  # y range with 25% headroom so the b* marker and label never leave the panel
  yr <- range(m$profile); span <- max(diff(yr), 1e-6)
  ylim <- c(yr[1] - 0.05 * span, yr[2] + 0.25 * span)
  plot(m$b_grid, m$profile, type = "l", lwd = 3, col = "#3B528B", bty = "l",
       ylim = ylim, xlab = "plasticity slope b",
       ylab = "weighted fitness selection maximizes",
       main = sprintf("Fitness vs. plasticity slope (b* = %.2f)", m$b_star))
  abline(v = m$b_star, lty = 2, col = "grey40")
  points(m$b_star, max(m$profile), pch = 19, cex = 1.4)
  text(m$b_star, yr[2] + 0.15 * span, sprintf("b* = %.2f", m$b_star),
       pos = if (m$b_star > 0.6) 2 else 4)
}

#' Fitness bar chart 
plot_fitness <- function(m) {
  mat <- rbind(Plastic = c(m$pH, m$pN, m$pMean),
               Canalized = c(m$cH, m$cN, m$cMean))
  colnames(mat) <- c("Harsh", "Nice", "Population mean")
  barplot(mat, beside = TRUE, ylim = c(0, 1), col = c("#3B528B", "#5DC863"),
          ylab = "Expected fitness", legend.text = rownames(mat),
          args.legend = list(x = "topleft", bty = "n"))
}

#' ggplot2 versions (used by the Shiny app); require ggplot2
gg_reaction_norms <- function(m) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  D <- m$inputs$delta; pm <- plastic_means(m)
  lines <- data.frame(
    genotype = rep(c("Plastic genotype (mean phenotype)", "Canalized genotype"), each = 2),
    env = rep(c(1, 2), 2),
    z = c(pm[["harsh"]], pm[["nice"]], m$z_star, m$z_star))
  opt <- data.frame(env = c(1, 2), z = c(0, D))
  ggplot2::ggplot() +
    ggplot2::geom_line(data = opt, ggplot2::aes(env, z), linetype = "dashed",
                       colour = "grey50") +
    ggplot2::geom_point(data = opt, ggplot2::aes(env, z), size = 3, colour = "grey30") +
    ggplot2::geom_text(data = opt, ggplot2::aes(env, z, label = "optimum"),
                       hjust = c(1.3, -0.3), colour = "grey30") +
    ggplot2::geom_line(data = lines, ggplot2::aes(env, z, colour = genotype),
                       linewidth = 1.4) +
    ggplot2::scale_colour_manual(values = c("Plastic genotype (mean phenotype)" = "#3B528B",
                                            "Canalized genotype" = "#5DC863")) +
    ggplot2::scale_x_continuous(breaks = c(1, 2), limits = c(0.6, 2.4),
                                labels = c("Harsh (low fitness)", "Nice (high fitness)")) +
    ggplot2::coord_cartesian(ylim = c(-1, D + 1)) +
    ggplot2::labs(x = NULL, y = "phenotype z", colour = NULL,
                  title = sprintf("Reaction norms (evolved slope b* = %.2f)", m$b_star)) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "top")
}

gg_slope_profile <- function(m) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  d <- data.frame(b = m$b_grid, fitness = m$profile)
  # Robust b* label (v2): explicit y-range with 40% headroom above the peak, and
  # the label pinned just inside the top edge of the panel (y = Inf, vjust > 1),
  # so it can never be clipped, whatever the plot height. No clip = "off".
  yr <- range(m$profile); span <- max(diff(yr), 1e-4)
  ylims <- c(yr[1] - 0.05 * span, yr[2] + 0.40 * span)
  right <- m$b_star > 0.6
  lab <- sprintf("b* = %.2f", m$b_star)
  ggplot2::ggplot(d, ggplot2::aes(b, fitness)) +
    ggplot2::geom_line(colour = "#3B528B", linewidth = 1.4) +
    ggplot2::geom_vline(xintercept = m$b_star, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_point(data = data.frame(b = m$b_star, fitness = max(m$profile)),
                        size = 3) +
    ggplot2::annotate("text", x = m$b_star, y = Inf, label = lab,
                      hjust = if (right) 1.1 else -0.1, vjust = 1.6, size = 4.5) +
    ggplot2::scale_x_continuous(limits = c(0, 1),
                                expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::coord_cartesian(ylim = ylims) +
    ggplot2::labs(x = "plasticity slope b",
                  y = "weighted fitness",
                  title = "Fitness vs. plasticity slope (axis does not start at zero)",
                  subtitle = sprintf("Evolved slope %s (dashed line)", lab)) +
    ggplot2::theme_minimal(base_size = 13)
}

gg_fitness <- function(m) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  d <- data.frame(
    env = factor(rep(c("Harsh", "Nice", "Population mean"), 2),
                 levels = c("Harsh", "Nice", "Population mean")),
    genotype = rep(c("Plastic genotype", "Best canalized genotype"), each = 3),
    fitness = c(m$pH, m$pN, m$pMean, m$cH, m$cN, m$cMean))
  ggplot2::ggplot(d, ggplot2::aes(env, fitness, fill = genotype)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", fitness)),
                       position = ggplot2::position_dodge(width = 0.9), vjust = -0.4) +
    ggplot2::scale_fill_manual(values = c("Plastic genotype" = "#3B528B",
                                          "Best canalized genotype" = "#5DC863")) +
    ggplot2::coord_cartesian(ylim = c(0, 1.05)) +
    ggplot2::labs(x = NULL, y = "Expected fitness", fill = NULL) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "top")
}

# ---- Example (runs only when this file is executed directly) --------------
if (sys.nframe() == 0) {
  for (cst in c(0, 0.05, 0.2)) {
    m <- run_model(f_nice = 0.7, ratio = 0.4, delta = 2, rho = 0.9,
                       cost = cst, regime = "hard")
    cat(sprintf("cost = %.2f  ->  b* = %.2f, a* = %.3f, canalized z* = %.3f\n",
                cst, m$b_star, m$a_star, m$z_star))
  }
  m <- run_model(f_nice = 0.7, ratio = 0.4, delta = 2, rho = 0.9,
                     cost = 0.05, regime = "hard")
  print(results_table(m), digits = 3)
  print(plasticity_table(m), digits = 3)
  cat("\n", m$verdict, "\n")
  png("example.png", width = 1800, height = 560, res = 130)
  op <- par(mfrow = c(1, 3)); plot_reaction_norms(m); plot_slope_profile(m)
  plot_fitness(m); par(op)
  dev.off()
  cat("Wrote example.png\n")
}
