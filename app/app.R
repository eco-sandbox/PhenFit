# app.R -- Shiny app for PhenFit
# Run with: shiny::runApp()   (from this folder)
library(shiny)
library(ggplot2)
source("model.R", local = TRUE)   # model + other plots; needs model.R in the same folder

# Fitness-vs-slope plot
# clipped b* label even next to an older model.R: it overrides gg_slope_profile().
gg_slope_profile <- function(m) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  d <- data.frame(b = m$b_grid, fitness = m$profile)
  # Robust b* label: explicit y-range with 40% headroom above the peak, and
  # the label pinned just inside the top edge of the panel (y = Inf, vjust > 1),
  # so it can never be clipped, whatever the plot height. No clip = "off".
  yr <- range(m$profile); span <- max(diff(yr), 1e-4)
  ylims <- c(yr[1] - 0.05 * span, yr[2] + 0.40 * span)
  right <- m$b_star > 0.6
  lab <- sprintf("b* = %.2f", m$b_star)
  ggplot2::ggplot(d, ggplot2::aes(b, fitness)) +
    ggplot2::geom_line(colour = "#3B528B", linewidth = 1.4) +   # viridis blue (plastic genotype)
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


# slider label with a smaller explanatory second line
lab2 <- function(main, sub) HTML(paste0(main,
  "<br><span style='font-weight:normal; font-size:12px; color:#666;'>", sub, "</span>"))

# Every case sets both selection widths (om_n, om_h); notes may contain HTML (<sub>).
cases <- list(
  "1" = list(name = "Baseline",        f_nice = 0.7, ratio = 0.4, delta = 2,   rho = 0.9, cost = 0.05, regime = "hard",
             om_n = 1, om_h = 1,
             note = "Default settings: a moderate slope evolves; the benefit appears in the harsh environment."),
  "2" = list(name = "Unreliable cue",  f_nice = 0.7, ratio = 0.4, delta = 2,   rho = 0.5, cost = 0.05, regime = "hard",
             om_n = 1, om_h = 1,
             note = "Cue reliability 0.5 (uninformative): the slope collapses to zero."),
  "3" = list(name = "High cost",       f_nice = 0.7, ratio = 0.4, delta = 2,   rho = 0.9, cost = 0.5,  regime = "hard",
             om_n = 1, om_h = 1,
             note = "Cost c = 0.5: a high cost pushes the population back to canalization."),
  "4" = list(name = "Harsh common",    f_nice = 0.1, ratio = 0.4, delta = 2,   rho = 0.9, cost = 0.05, regime = "hard",
             om_n = 1, om_h = 1,
             note = "Nice environment rare (0.1): the benefit of plasticity moves to the nice environment."),
  "5" = list(name = "Sink-like harsh", f_nice = 0.5, ratio = 0.2, delta = 2,   rho = 0.9, cost = 0.05, regime = "hard",
             om_n = 1, om_h = 1,
             note = "Harsh/nice fitness 0.2: trait means follow the nice optimum; the gain shows up in the harsh habitat."),
  "6" = list(name = "Small \u0394",    f_nice = 0.7, ratio = 0.4, delta = 0.5, rho = 0.9, cost = 0.05, regime = "hard",
             om_n = 1, om_h = 1,
             note = "Optimum difference 0.5: little plasticity and almost no benefit, despite a large fitness difference."),
  "7" = list(name = "Distant optima",  f_nice = 0.7, ratio = 0.4, delta = 2.5, rho = 0.9, cost = 0.05, regime = "hard",
             om_n = 1, om_h = 1,
             note = "Optimum difference 2.5: fitness vs. slope has two peaks, a near-canalized nice specialist (b \u2248 0.1) and the plastic optimum (b* \u2248 0.82); intermediate slopes do worse than both."),
  "8" = list(name = "Unequal widths",  f_nice = 0.7, ratio = 0.4, delta = 2,   rho = 0.9, cost = 0.05, regime = "hard",
             om_n = 0.5, om_h = 1,
             note = "Narrow nice niche (\u03c9<sub>nice</sub> = 0.5): fitness vs. slope has two peaks, but here the near-canalized nice specialist wins (b* \u2248 0.06); the plastic peak (b \u2248 0.83) is lower. The reverse of case 7.")
)


phenfit_css <- "
.text-col { max-width: 850px; margin: 18px 0 40px 0; font-size: 15px; line-height: 1.55; }
.text-col h3 { margin-top: 28px; }
.text-col figure { margin: 18px 0 24px 0; }
.text-col figure img { width: 100%; max-width: 700px; height: auto; display: block; }
.text-col figcaption { max-width: 700px; font-size: 13px; color: #555; margin-top: 6px; }
.text-col .fig-credit { max-width: 700px; font-size: 11px; color: #888; margin-top: 2px; }
.glossary { width: 100%; border-collapse: collapse; font-size: 14px; }
.glossary th, .glossary td { border-bottom: 1px solid #ddd; padding: 6px 8px; vertical-align: top; text-align: left; }
.glossary th { background: #f5f5f5; border-bottom: 2px solid #bbb; }
.glossary td:first-child { font-weight: 600; width: 22%; }
.glossary td:last-child { width: 26%; color: #444; }
.model-terms td:last-child { width: auto; color: inherit; }
.refs p { margin: 0 0 8px 0; padding-left: 1.5em; text-indent: -1.5em; font-size: 14px; }
/* sliders in viridis blue (#3B528B, the plastic-genotype color) instead of the default Bootstrap blue */
.irs-bar, .irs-bar-edge, .irs-single, .irs-from, .irs-to { background: #3B528B !important; border-color: #3B528B !important; }
.irs-single:before, .irs-from:before, .irs-to:before { border-top-color: #3B528B !important; }
"

# Re-typeset MathJax when a tab is shown (equations in hidden tabs can otherwise
# be laid out at zero width).
phenfit_js <- "
$(document).on('shown.bs.tab', 'a[data-toggle=\"tab\"], a[data-bs-toggle=\"tab\"]', function() {
  if (window.MathJax && MathJax.Hub) { MathJax.Hub.Queue(['Typeset', MathJax.Hub]); }
});
"

fig_block <- function(src, alt, caption, credit = NULL) {
  tags$figure(tags$img(src = src, alt = alt), tags$figcaption(caption),
              if (!is.null(credit)) tags$div(class = "fig-credit", credit))
}

# Turn every DOI URL in a reference into a link that opens in a new tab
ref_par <- function(x) {
  x <- htmltools::htmlEscape(x)
  x <- gsub("(https://doi\\.org/[^ ]+)",
            "<a href=\"\\1\" target=\"_blank\" rel=\"noopener noreferrer\">\\1</a>", x)
  p(HTML(x))
}

model_terms <- list(
  c("Trait", "the measured character, e.g. specific leaf area (SLA)"),
  c("Phenotype (\\(z\\))", "the trait value an individual expresses"),
  c("Genotype (\\(a, b\\))", "the heritable reaction norm; plastic if \\(b > 0\\), canalized if \\(b = 0\\)"),
  c("Optimal trait value (\\(\\theta\\))", "the trait value with the highest fitness in a given environment"),
  c("Cost-free regression slope (\\(\\beta\\))", "the regression of \\(\\theta\\) on \\(s\\), in trait units, i.e. the best slope without costs if fitness were quadratic"),
  c("Intercept (\\(a\\))", "the phenotype a genotype expresses when the cue says harsh (harsh is the reference state)"),
  c("Slope (\\(b\\))", "the strength of plasticity: the change in phenotype when the cue says nice, as a fraction of \\(\\Delta\\)"),
  c("Evolved genotype (\\(a^*, b^*\\))", "the genotype with the highest mean fitness across both environments, found by grid search"),
  c("Harsh/nice maximum-fitness ratio (\\(r\\))", "maximum attainable fitness in the harsh environment relative to the nice one (0–1), matters only under hard selection"),
  c("Distance between optima (\\(\\Delta\\))", "the difference between the optimal trait values in the nice and the harsh environment, in units of the selection width"),
  c("Selection width (\\(\\omega\\))", "the width of the Gaussian fitness curve, at \\(\\theta \\pm \\omega\\) fitness is 61% of its maximum"),
  c("Maximum fitness (\\(W_{max}\\))", "the highest fitness an environment allows: 1 in the nice environment and \\(r\\) in the harsh one"),
  c("Frequency of nice (\\(f\\))", "the probability that a generation meets the nice environment"),
  c("Cue (\\(s\\))", "the environmental signal a plant senses before selection acts: \\(s = 1\\) if it says nice, \\(s = 0\\) if it says harsh"),
  c("Cue reliability (\\(\\rho\\))", "the probability that the cue correctly indicates the environment"),
  c("Cost of plasticity (\\(c\\))", "scales the fitness cost of plasticity, which enters as \\((1 - c b^2)\\)"),
  c("Selection regime", "hard: each environment counts by its frequency and maximum fitness, soft: by frequency alone")
)

model_terms_table <- tags$table(
  class = "glossary model-terms",
  tags$thead(tags$tr(tags$th("Term"), tags$th("Meaning in PhenFit"))),
  tags$tbody(lapply(model_terms, function(r)
    tags$tr(tags$td(r[1]), tags$td(r[2]))))
)

glossary <- list(
  c("Active (anticipatory) plasticity", "the plant senses a cue and a regulated developmental response follows, such as shade avoidance triggered by a low red to far-red ratio, which changes stem elongation before any resource is limiting", "Smith-Gill 1983; Forsman 2015"),
  c("Bet-hedging", "a strategy that lowers fitness variance across unpredictable environments at the cost of a lower mean, e.g. by producing offspring with randomly different phenotypes rather than responding to a cue", "Simons 2011; Starrfelt and Kokko 2012"),
  c("Canalization", "the production of the same phenotype despite environmental (or genetic) variation; a canalized genotype has a flat reaction norm", "Waddington 1953; Van Buskirk and Steiner 2009"),
  c("Costs and limits of plasticity", "a cost is a fitness reduction paid by a plastic genotype even when it produces the same phenotype as a non-plastic one. A limit is a constraint that stops plasticity from producing the optimal phenotype", "DeWitt et al. 1998; Auld et al. 2010; Murren et al. 2015"),
  c("Cryptic genetic variation", "genetic variation that has no phenotypic effect under usual conditions but is expressed under new or stressful ones", "Rutherford and Lindquist 1998; Paaby and Rockman 2014"),
  c("Cue reliability", "how well an environmental cue predicts the conditions under which selection acts", "Tufto 2000; Reed et al. 2010; Chevin and Lande 2015"),
  c("Environmental grain", "fine-grained when individuals experience several environments in their lifetime, coarse-grained when each experiences only one", "Levins 1968; Baythavong 2011"),
  c("Genetic assimilation", "the evolution of an initially plastic, environmentally induced phenotype into one that is produced regardless of the environment", "Waddington 1953; Pigliucci et al. 2006; Ehrenreich and Pfennig 2016"),
  c("Genetic tracking", "the evolution of the mean trait, through selection across generations. It works when autocorrelation is high: the parents' environment predicts the offspring's, so the population can adapt to it over time", "Tufto 2015; Botero et al. 2015"),
  c("Genotype-by-environment interaction (G\u00d7E)", "genetic variation in plasticity, meaning that genotypes differ in reaction-norm slope or shape", "Via and Lande 1985; Scheiner 1993"),
  c("Hard and soft selection", "under hard selection, each environment's contribution to the next generation scales with its productivity. Under soft selection, each contributes a fixed share regardless of productivity", "Van Tienderen 1991; de Jong 1995"),
  c("Jack-of-all-trades / Master-of-some", "plastic genotypes that maintain fitness in unfavorable environments, versus ones that raise fitness in favorable environments", "Richards et al. 2006"),
  c("Passive plasticity", "phenotypic change that follows directly from the environment's effect on growth or physiology (e.g. smaller plants under nutrient shortage) rather than from a regulated, cue-driven response. It is often not adaptive.", "Van Kleunen and Fischer 2005; Forsman 2015"),
  c("Phenotypic plasticity", "the ability of one genotype to produce different phenotypes in response to different environments (in contrast to bet-hedging, where variation arises without an environmental cue)", "Johannsen 1911; Bradshaw 1965; Schlichting 1986; Pigliucci 2005"),
  c("Reaction norm", "the set of phenotypes a genotype expresses across a range of environments, often summarized by its intercept and slope", "Woltereck 1909; Via and Lande 1985; Gomulkiewicz and Kirkpatrick 1992"),
  c("Stabilizing selection", "selection favoring an intermediate optimum phenotype, with fitness falling off on both sides", "Lande and Arnold 1983; Via and Lande 1985")
)

glossary_table <- tags$table(
  class = "glossary",
  tags$thead(tags$tr(tags$th("Term"), tags$th("Definition"), tags$th("Reference"))),
  tags$tbody(lapply(glossary, function(r)
    tags$tr(tags$td(r[1]), tags$td(r[2]), tags$td(r[3]))))
)

background_refs <- c(
  "Auld, Agrawal & Relyea (2010). Re-evaluating the costs and limits of adaptive phenotypic plasticity. Proc R Soc B 277:503-511. https://doi.org/10.1098/rspb.2009.1355",
  "Baythavong (2011). Linking the spatial scale of environmental variation and the evolution of phenotypic plasticity: selection favors adaptive plasticity in fine-grained environments. Am Nat 178:75-87. https://doi.org/10.1086/660281",
  "Botero, Weissing, Wright & Rubenstein (2015). Evolutionary tipping points in the capacity to adapt to environmental change. PNAS 112:184-189. https://doi.org/10.1073/pnas.1408589111",
  "Bradshaw (1965). Evolutionary significance of phenotypic plasticity in plants. Advances in Genetics 13:115-155. https://doi.org/10.1016/S0065-2660(08)60048-6",
  "Chevin & Lande (2015). Evolution of environmental cues for phenotypic plasticity. Evolution 69:2767-2775. https://doi.org/10.1111/evo.12755",
  "de Jong (1995). Phenotypic plasticity as a product of selection in a variable environment. Am Nat 145:493-512. https://doi.org/10.1086/285752",
  "DeWitt, Sih & Wilson (1998). Costs and limits of phenotypic plasticity. Trends Ecol Evol 13:77-81. https://doi.org/10.1016/S0169-5347(97)01274-3",
  "Ehrenreich & Pfennig (2016). Genetic assimilation: a review of its potential proximate causes and evolutionary consequences. Ann Bot 117:769-779. https://doi.org/10.1093/aob/mcv130",
  "Forsman (2015). Rethinking phenotypic plasticity and its consequences for individuals, populations and species. Heredity 115:276-284. https://doi.org/10.1038/hdy.2014.92",
  "Gomulkiewicz & Kirkpatrick (1992). Quantitative genetics and the evolution of reaction norms. Evolution 46:390-411. https://doi.org/10.1111/j.1558-5646.1992.tb02047.x",
  "Holt & Gaines (1992). Analysis of adaptation in heterogeneous landscapes: implications for the evolution of fundamental niches. Evol Ecol 6:433-447. https://doi.org/10.1007/BF02270702",
  "Johannsen (1911). The genotype conception of heredity. Am Nat 45:129-159. https://doi.org/10.1086/279202",
  "Kawecki (1995). Demography of source-sink populations and the evolution of ecological niches. Evol Ecol 9:38-44. https://doi.org/10.1007/BF01237695",
  "Lande (2009). Adaptation to an extraordinary environment by evolution of phenotypic plasticity and genetic assimilation. J Evol Biol 22:1435-1446. https://doi.org/10.1111/j.1420-9101.2009.01754.x",
  "Lande & Arnold (1983). The measurement of selection on correlated characters. Evolution 37:1210-1226. https://doi.org/10.1111/j.1558-5646.1983.tb00236.x",
  "Levins R (1968). Evolution in Changing Environments: Some Theoretical Explorations. Princeton University Press, Princeton, NJ.",
  "Murren, Auld, Callahan, Ghalambor et al. (2015). Constraints on the evolution of phenotypic plasticity: limits and costs of phenotype and plasticity. Heredity 115:293-301. https://doi.org/10.1038/hdy.2015.8",
  "Paaby & Rockman (2014). Cryptic genetic variation: evolution's hidden substrate. Nat Rev Genet 15:247-258. https://doi.org/10.1038/nrg3688",
  "Pigliucci (2005). Evolution of phenotypic plasticity: where are we going now? Trends Ecol Evol 20:481-486. https://doi.org/10.1016/j.tree.2005.06.001",
  "Pigliucci, Murren & Schlichting (2006). Phenotypic plasticity and evolution by genetic assimilation. J Exp Biol 209:2362-2367. https://doi.org/10.1242/jeb.02070",
  "Reed, Waples, Schindler, Hard & Kinnison (2010). Phenotypic plasticity and population viability: the importance of environmental predictability. Proc R Soc B 277:3391-3400. https://doi.org/10.1098/rspb.2010.0771",
  "Richards, Bossdorf, Muth, Gurevitch & Pigliucci (2006). Jack of all trades, master of some? On the role of phenotypic plasticity in plant invasions. Ecol Lett 9:981-993. https://doi.org/10.1111/j.1461-0248.2006.00950.x",
  "Rutherford & Lindquist (1998). Hsp90 as a capacitor for morphological evolution. Nature 396:336-342. https://doi.org/10.1038/24550",
  "Scheiner (1993). Genetics and evolution of phenotypic plasticity. Annu Rev Ecol Syst 24:35-68. https://doi.org/10.1146/annurev.es.24.110193.000343",
  "Schlichting (1986). The evolution of phenotypic plasticity in plants. Annu Rev Ecol Syst 17:667-693. https://doi.org/10.1146/annurev.es.17.110186.003315",
  "Simons (2011). Modes of response to environmental change and the elusive empirical evidence for bet hedging. Proc R Soc B 278:1601-1609. https://doi.org/10.1098/rspb.2011.0176",
  "Smith-Gill (1983). Developmental plasticity: developmental conversion versus phenotypic modulation. Am Zool 23:47-55. https://doi.org/10.1093/icb/23.1.47",
  "Starrfelt & Kokko (2012). Bet-hedging: a triple trade-off between means, variances and correlations. Biol Rev 87:742-755. https://doi.org/10.1111/j.1469-185X.2012.00225.x",
  "Tufto (2000). The evolution of plasticity and nonplastic spatial and temporal adaptations in the presence of imperfect environmental cues. Am Nat 156:121-130. https://doi.org/10.1086/303381",
  "Tufto (2015). Genetic evolution, plasticity, and bet-hedging as adaptive responses to temporally autocorrelated fluctuating selection: a quantitative genetic model. Evolution 69:2034-2049. https://doi.org/10.1111/evo.12716",
  "Van Buskirk & Steiner (2009). The fitness costs of developmental canalization and plasticity. J Evol Biol 22:852-860. https://doi.org/10.1111/j.1420-9101.2009.01685.x",
  "Van Kleunen & Fischer (2005). Constraints on the evolution of adaptive phenotypic plasticity in plants. New Phytol 166:49-60. https://doi.org/10.1111/j.1469-8137.2004.01296.x",
  "Van Tienderen (1991). Evolution of generalists and specialists in spatially heterogeneous environments. Evolution 45:1317-1331. https://doi.org/10.1111/j.1558-5646.1991.tb02638.x",
  "Via & Lande (1985). Genotype-environment interaction and the evolution of phenotypic plasticity. Evolution 39:505-522. https://doi.org/10.1111/j.1558-5646.1985.tb00391.x",
  "Waddington (1953). Genetic assimilation of an acquired character. Evolution 7:118-126. https://doi.org/10.1111/j.1558-5646.1953.tb00070.x",
  "Woltereck (1909). Weitere experimentelle Untersuchungen \u00fcber Artver\u00e4nderung, speziell \u00fcber das Wesen quantitativer Artunterschiede bei Daphniden. Verhandlungen der Deutschen Zoologischen Gesellschaft 19:110-173."
)

tab_background <- withMathJax(div(
  class = "text-col",
  h3("Glossary"),
  glossary_table,

  h3("Model terms"),
  model_terms_table,

  h3("Gaussian stabilizing selection"),
  p("In each environment, fitness is highest at the optimal trait value \\(\\theta\\) and falls off symmetrically on both sides:"),
  p("$$W(z) = W_{max} \\exp\\left[-\\frac{(z-\\theta)^2}{2\\omega^2}\\right]$$"),
  p("\\(W_{max}\\) is the highest fitness the environment allows. \\(\\omega\\) is the width of the fitness curve:",
    "a large \\(\\omega\\) means weak selection, and a small \\(\\omega\\) means strong selection.",
    "A plant one width away from the optimum (\\(z = \\theta \\pm \\omega\\)) keeps 61% of \\(W_{max}\\).",
    "At half a width it keeps 88%, at two widths 14%, and at three widths only 1% (Via and Lande 1985; Lande 2009)."),
  fig_block("fig_gaussian.png", "Gaussian fitness function",
            "Gaussian fitness function. Points mark fitness at 0.5, 1, 2 and 3 widths from the optimum. ©PhenFit"),

  h3("Reaction norms"),
  p("A reaction norm describes how a genotype's phenotype changes with its environment.",
    "With two environments it is a straight line: the intercept \\(a\\) is the phenotype in the reference (harsh) environment,",
    "and the slope \\(b\\) is how much the phenotype changes in the other one.",
    "A flat line (\\(b = 0\\)) is a canalized genotype, and a sloped line (\\(b > 0\\)) is a plastic one."),

  h3("Cues and their reliability"),
  p("Plants can't sense the conditions under which they'll be selected, only a cue that predicts them,",
    "such as a low red to far-red ratio predicting shade. The cue reliability \\(\\rho\\) is the probability that the cue is right.",
    "Without any cost, and if fitness were quadratic, the best slope would be the cost-free regression slope \\(\\beta\\) of the optimal trait value on the cue:"),
  p("$$\\beta = \\frac{\\mathrm{cov}(\\theta, s)}{\\mathrm{var}(s)}$$"),
  p("A perfect cue (\\(\\rho = 1\\)) favors tracking the full distance between the optima.",
    "A cue at chance level (\\(\\rho = 0.5\\)) favors no plasticity at all (Tufto 2000; Chevin and Lande 2015)."),
  fig_block("anim_cue_regression.gif", "Cost-free regression slope \\(\\beta\\) for different cue reliabilities",
            "Cost-free regression slope \\(\\beta\\) for different cue reliabilities \u03c1. Bubble area shows probability. ©PhenFit"),

  h3("Costs and limits of plasticity"),
  p("Being plastic can cost fitness, for example to maintain sensing and regulation.",
    "That cost holds plasticity below the cost-free optimum (DeWitt et al. 1998; Auld et al. 2010)."),

  h3("Hard vs soft selection"),
  p("Under hard selection, a productive environment contributes more offspring to the next generation, so selection there counts more.",
    "Under soft selection, each environment contributes a fixed share whatever its productivity.",
    "Hard selection therefore pulls adaptation toward the nice environment, which is the source-sink logic",
    "(Holt and Gaines 1992; Kawecki 1995; Van Tienderen 1991)."),
  fig_block("Fig_hard_vs_soft_pool.png", "Hard and soft selection",
            "Hard selection: patches contribute offspring by productivity. Soft selection: each patch contributes a fixed share. ©PhenFit"),

  h3("References"),
  div(class = "refs", lapply(background_refs, ref_par))
))

# Preset-case lines: one per case button in `cases` above (checked by the stopifnot below).
# Entries may contain HTML (<sub>); they are rendered with HTML().
case_lines <- c(
  "1" = "Baseline: the default settings (<i>f</i> = 0.7, r = 0.4, \u0394 = 2, \u03c1 = 0.9, c = 0.05, hard selection)",
  "2" = "Unreliable cue: \u03c1 = 0.5",
  "3" = "High cost: c = 0.5",
  "4" = "Harsh common: <i>f</i> = 0.1",
  "5" = "Sink-like harsh: <i>f</i> = 0.5, r = 0.2",
  "6" = "Small \u0394: \u0394 = 0.5",
  "7" = "Distant optima: \u0394 = 2.5 (fitness vs. slope has two peaks)",
  "8" = "Unequal widths: \u03c9<sub>nice</sub> = 0.5 (near-canalized peak wins over the plastic one)"
)
stopifnot(identical(names(case_lines), names(cases)),
          all(vapply(cases, function(cs) is.numeric(cs$om_n) && is.numeric(cs$om_h), logical(1))))

tab_howitworks <- withMathJax(div(
  class = "text-col",
  h3("What it does"),
  p("PhenFit asks how much plasticity selection favors when a population lives in a mix of harsh and nice environments.",
    "It finds the evolved genotype \\((a^*, b^*)\\) and compares it with the best canalized genotype (\\(b = 0\\))."),

  h3("Environments"),
  p("Each generation, a plant meets the nice environment with probability \\(f\\) and the harsh one with probability \\(1 - f\\).",
    "Each generation is drawn independently, so there is no temporal autocorrelation.",
    "The optimal trait value is 0 in the harsh environment and \\(\\Delta\\) in the nice one.",
    "\\(\\Delta\\) is measured in units of the selection width (\\(\\omega = 1\\))."),

  h3("Phenotype"),
  p("A plant receives a binary cue that is correct with probability \\(\\rho\\).",
    "If the cue says nice, it expresses \\(z = a + b\\Delta\\). If the cue says harsh, it expresses \\(z = a\\)."),

  h3("Fitness in each environment"),
  p("$$W_e(a,b) = W_{max,e}\\,(1 - c b^2)\\; E\\left[\\exp\\left(-\\frac{(z-\\theta_e)^2}{2\\omega_e^2}\\right)\\right]$$"),
  p("This applies Gaussian selection separately in each environment and averages over both cue outcomes.",
    "The factor \\((1 - cb^2)\\) is the cost of plasticity. Because the cost rises with \\(b^2\\),",
    "a little plasticity is almost free, but strong plasticity becomes expensive."),
  fig_block("anim_b_sweep.gif", "Fitness landscape in the harsh and the nice environment",
            HTML(paste("Mean fitness of the plastic genotype peaks at an intermediate slope (b* = 0.58), because plasticity pays in the harsh environment but costs fitness and misfires when the cue is wrong. ©PhenFit"))),
  h3("Mean fitness"),
  p("$$\\bar{W} = \\sum_e w_e W_e$$"),
  p("Under hard selection, \\(w_e\\) is the environment's frequency \u00d7 its \\(W_{max}\\).",
    "Under soft selection, \\(w_e\\) is the frequency only."),

  h3("Finding b*"),
  p("PhenFit searches a fine grid of \\((a, b)\\) values for the genotype with the highest \\(\\bar{W}\\).",
    "It runs the same search with \\(b\\) fixed at 0 to find the best canalized genotype."),

  h3("Controls"),
  tags$ul(
    tags$li("Frequency of the nice environment (\\(f\\)): 0.05\u20130.95 (default 0.7)"),
    tags$li("Harsh/nice maximum-fitness ratio (r): 0.1\u20131 (default 0.4)"),
    tags$li("Distance between optima (\\(\\Delta\\)): 0.5\u20133 (default 2)"),
    tags$li("Cue reliability (\\(\\rho\\)): 0.5\u20131 (default 0.9)"),
    tags$li("Cost of plasticity (\\(c\\)): 0\u20130.5 (default 0.05)"),
    tags$li("Selection width in the nice and in the harsh environment (\\(\\omega_{nice}\\), \\(\\omega_{harsh}\\)): 0.25\u20133 each (default 1)"),
    tags$li("Selection regime: hard or soft selection (default hard)")
  ),

  h3("Outputs"),
  p("Top panel: reaction norms of the plastic and canalized genotypes against the two optima.",
    "Middle panel: best mean fitness as a function of \\(b\\), with \\(b^*\\) marked.",
    "Bottom panel: expected fitness of the plastic and the best canalized genotype in each environment and the population mean.",
    "Below the panels, tables list these fitness values with the relative and absolute gains of plasticity",
    "and a summary of the evolved genotype, followed by a short summary."),

  h3("Preset cases"),
  p("Eight preset scenarios load parameter combinations:"),
  tags$ol(lapply(unname(case_lines), function(x) tags$li(HTML(x)))),

  h3("Assumptions and limits"),
  p("PhenFit is illustrative, not predictive. It models a single trait, two discrete environments, a binary cue,",
    "a quadratic cost and no temporal autocorrelation, with no demography, gene flow, genetic variance or evolutionary dynamics.",
    "It finds the optimal genotype but doesn't track how a population gets there."),
  
  h3("Attribution"),
  p("PhenFit (v1.9 / 10.2026) © 2026 [Michel J] BSD-2 / CC BY 4.0."),
  
))

ui <- fluidPage(
  tags$head(tags$style(HTML(phenfit_css)), tags$script(HTML(phenfit_js))),
  titlePanel("PhenFit: Finding the fittest genotype"),
  p("Simplified model of selection on plastic and canalized genotypes in harsh vs. nice environments (v1.9 / 10.2026)"),
  tabsetPanel(
    id = "main_tabs",
    tabPanel(
      "Model",
      sidebarLayout(
        sidebarPanel(
          tags$b("Example cases"),
          div(style = "margin: 6px 0 4px 0;",
              lapply(names(cases), function(k)
                actionButton(paste0("case_", k), paste(k, cases[[k]]$name),
                             class = "btn-sm", style = "margin: 2px;"))),
          helpText(uiOutput("case_note", inline = TRUE)),
          hr(),
          sliderInput("f_nice", HTML("Frequency of nice environment (<i>f</i>)"), 0.05, 0.95, 0.7, step = 0.05),
          sliderInput("ratio", "Harsh/nice maximum-fitness ratio (r)", 0.1, 1, 0.4, step = 0.1),
          sliderInput("delta", lab2("Optimum difference \u0394", "in units of width \u03c9 = 1"), 0.5, 3, 2, step = 0.5),
          sliderInput("rho", lab2("Cue reliability \u03c1", "0.5 = uninformative"), 0.5, 1, 0.9, step = 0.05),
          sliderInput("om_n", lab2("Selection width in nice environment (\u03c9<sub>nice</sub>)", "larger width = weaker selection"), 0.25, 3, 1, step = 0.25),
          sliderInput("om_h", lab2("Selection width in harsh environment (\u03c9<sub>harsh</sub>)", "smaller width = stronger selection"), 0.25, 3, 1, step = 0.25),
          sliderInput("cost", lab2("Cost of plasticity c", "fitness \u00d7 (1 \u2212 c\u00b7b\u00b2)"), 0, 0.5, 0.05, step = 0.05),
          radioButtons("regime", "Selection regime",
                       choices = c("Hard: frequency \u00d7 maximum fitness" = "hard",
                                   "Soft: frequency only" = "soft"),
                       selected = "hard")
        ),
        mainPanel(
          plotOutput("rn", height = 320),
          plotOutput("slope", height = 260),
          plotOutput("fit", height = 280),
          tableOutput("tab"),
          tableOutput("ptab"),
          textOutput("verdict")
        )
      )
    ),
    tabPanel("Theoretical background", tab_background),
    tabPanel("How PhenFit works", tab_howitworks)
  )
)

server <- function(input, output, session) {
  case_note <- reactiveVal("Click a case to load settings, change sliders to explore.")
  lapply(names(cases), function(k) {
    observeEvent(input[[paste0("case_", k)]], {
      cs <- cases[[k]]
      updateSliderInput(session, "f_nice", value = cs$f_nice)
      updateSliderInput(session, "ratio",  value = cs$ratio)
      updateSliderInput(session, "delta",  value = cs$delta)
      updateSliderInput(session, "rho",    value = cs$rho)
      updateSliderInput(session, "cost",   value = cs$cost)
      updateSliderInput(session, "om_n",   value = cs$om_n)   # every case sets both widths,
      updateSliderInput(session, "om_h",   value = cs$om_h)   # so case 1 after case 8 restores 1
      updateRadioButtons(session, "regime", selected = cs$regime)
      case_note(paste0("Case ", k, " (", cs$name, "): ", cs$note))
    })
  })
  output$case_note <- renderUI(HTML(case_note()))
  m <- reactive(run_model(input$f_nice, input$ratio, input$delta,
                              input$rho, input$cost, input$regime,
                              om_h = input$om_h, om_n = input$om_n))
  output$rn <- renderPlot(gg_reaction_norms(m()))
  output$slope <- renderPlot(gg_slope_profile(m()))
  output$fit <- renderPlot(gg_fitness(m()))
  output$ptab <- renderTable({
    d <- plasticity_table(m())
    d$Value <- sprintf("%.3f", d$Value)
    d
  })
  output$tab <- renderTable({
    d <- results_table(m())
    d$Relative_gain <- sprintf("%.0f%%", d$Relative_gain * 100)
    d$Plastic <- sprintf("%.2f", d$Plastic)
    d$Canalized <- sprintf("%.2f", d$Canalized)
    d$Absolute_gain <- sprintf("%.2f", d$Absolute_gain)
    names(d) <- c("Environment", "Plastic", "Canalized", "Relative gain", "Absolute gain")
    d
  })
  output$verdict <- renderText(m()$verdict)
}

shinyApp(ui, server)
