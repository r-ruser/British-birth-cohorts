options(stringsAsFactors = FALSE)
Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(survey)
  library(scales)
})
options(survey.lonely.psu = "adjust")

root <- "E:/公共数据库/P1"
v2 <- file.path(root, "life_course_v2_20260928")
v3 <- file.path(root, "life_course_v3_20260928")
sev <- file.path(root, "severe_asthma_relevant_analysis_20260928")
out <- file.path(root, "manuscript_figures_20260929")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_table <- function(folder, name) read.csv(file.path(folder, name), check.names = FALSE)
write_source <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE)
blue <- "#527C9C"
rust <- "#B27F67"
grey <- "#707A82"
theme_set(theme_classic(base_size = 8.5, base_family = "Arial") +
            theme(axis.line = element_line(linewidth = .35),
                  axis.ticks = element_line(linewidth = .35),
                  strip.background = element_blank(),
                  strip.text = element_text(face = "bold"),
                  legend.position = "bottom",
                  plot.tag = element_text(face = "bold", size = 11)))
save_figure <- function(plot, name, width_mm = 183, height_mm = 105) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  base <- file.path(out, name)
  ggsave(paste0(base, ".pdf"), plot, width = width_in, height = height_in,
         device = cairo_pdf, bg = "white")
  ggsave(paste0(base, ".svg"), plot, width = width_in, height = height_in,
         device = svglite::svglite, bg = "white")
  ggsave(paste0(base, ".png"), plot, width = width_in, height = height_in,
         dpi = 500, device = ragg::agg_png, bg = "white")
  ggsave(paste0(base, ".tiff"), plot, width = width_in, height = height_in,
         dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
}

# Figure 1 already matches the manuscript caption. Preserve the original exports.
for (ext in c("png", "pdf", "svg")) {
  from <- file.path(v3, "figures", paste0("01_cohort_equivalence.", ext))
  stopifnot(file.exists(from))
  file.copy(from, file.path(out, paste0("Figure_1_cohort_equivalence.", ext)), overwrite = TRUE)
}

# Figure 2: fit the published survey model once more to obtain the covariance of
# each main exposure coefficient and its interaction with previous symptoms.
long <- readRDS(file.path(v2, "analysis_longitudinal_frames.rds"))
z <- droplevels(long$MCS$monotone[long$MCS$monotone$age > 3, ])
terms <- c("male", "low_birthweight", "ses", "pregnancy_continued_smoke",
           "respondent_asthma", "atopy3", "bmi3_z", "ever_asthma3")
labels <- c("Male sex", "Low birthweight", "Early socioeconomic disadvantage",
            "Prenatal smoking proxy", "Parental asthma", "Atopy at age 3",
            "BMI at age 3 (per SD)", "Ever-asthma report by age 3")
names(labels) <- terms
design <- svydesign(ids = ~family, strata = ~design_stratum, weights = ~w,
                    data = z, nest = TRUE)
formula <- as.formula(paste("y ~ previous * (wave +", paste(terms, collapse = " + "), ")"))
fit <- svyglm(formula, design, family = quasipoisson("log"))
b <- coef(fit)
V <- vcov(fit)
stage <- read_table(v3, "stage_interactions.csv")
stage <- stage[stage$method == "stage_specific_wave", ]
records <- list()
for (term in terms) {
  interaction <- intersect(c(paste0("previous:", term), paste0(term, ":previous")), names(b))
  stopifnot(length(interaction) == 1)
  a <- b[term]
  c <- b[interaction]
  se_a <- sqrt(V[term, term])
  se_p <- sqrt(V[term, term] + V[interaction, interaction] + 2 * V[term, interaction])
  se_r <- sqrt(V[interaction, interaction])
  records[[term]] <- rbind(
    data.frame(term, factor = labels[term], measure = "Emergence", estimate = exp(a),
               lower = exp(a - 1.96 * se_a), upper = exp(a + 1.96 * se_a)),
    data.frame(term, factor = labels[term], measure = "Persistence", estimate = exp(a + c),
               lower = exp(a + c - 1.96 * se_p), upper = exp(a + c + 1.96 * se_p)),
    data.frame(term, factor = labels[term], measure = "Stage ROR", estimate = exp(-c),
               lower = exp(-c - 1.96 * se_r), upper = exp(-c + 1.96 * se_r))
  )
  published <- stage[stage$term == term, ]
  stopifnot(nrow(published) == 1,
            abs(exp(a) - published$onset_RR) < 1e-8,
            abs(exp(a + c) - published$persistence_RR) < 1e-8,
            abs(exp(-c) - published$stage_ROR) < 1e-8)
}
stage_data <- do.call(rbind, records)
write_source(stage_data, "Figure_2_source.csv")
stage_data$factor <- factor(stage_data$factor, levels = rev(unname(labels)))
stage_data$measure <- factor(stage_data$measure,
                             levels = c("Emergence", "Persistence", "Stage ROR"))
p2 <- ggplot(stage_data, aes(x = estimate, y = factor)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey65", linewidth = .35) +
  geom_errorbar(aes(xmin = lower, xmax = upper, colour = measure),
                orientation = "y", width = .12, linewidth = .5) +
  geom_point(aes(colour = measure), size = 1.8) +
  facet_grid(. ~ measure) +
  scale_colour_manual(values = c("Emergence" = blue, "Persistence" = rust,
                                 "Stage ROR" = grey), guide = "none") +
  scale_x_log10(breaks = c(.75, 1, 1.5, 2, 3), labels = c("0.75", "1", "1.5", "2", "3")) +
  labs(x = "Adjusted RR or stage ratio of RRs (95% CI)", y = NULL) +
  theme(panel.spacing.x = grid::unit(5, "mm"),
        axis.text.y = element_text(size = 7.4))
save_figure(p2, "Figure_2_transition_stages", 183, 100)

# Figure 3A: standardized MCS risk by childhood symptom history.
risk <- read_table(sev, "MCS_prior_history_standardized_risks.csv")
contrast <- read_table(sev, "04_childhood_burden_to_severe_wheeze.csv")
contrast <- contrast[contrast$term == "repeated_vs_isolated", ]
stopifnot(nrow(contrast) == 1)
risk$prior_group <- factor(risk$prior_group,
                           levels = c("never", "isolated", "repeated"),
                           labels = c("0", "1", expression("\u22652")))
write_source(risk, "Figure_3A_source.csv")
p3a <- ggplot(risk, aes(prior_group, adjusted_risk)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .12,
                linewidth = .6, colour = blue) +
  geom_point(size = 2.6, colour = blue) +
  scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, .26),
                     breaks = seq(0, .25, .05)) +
  labs(x = "Positive wheeze waves at ages 3, 5, 7",
       y = "Age-11 severe-wheeze proxy risk (95% CI)", title = "MCS") +
  annotate("text", x = 1.05, y = .245, hjust = 0, size = 2.55,
           label = sprintf("Repeated vs isolated: RR %.2f (%.2f\u2013%.2f)",
                           contrast$estimate, contrast$lower, contrast$upper))

# Figure 3B: average repeated adult reporting only; the six-wave endpoint belongs
# to Figure 4 and is deliberately excluded here.
adult <- read_table(v3, "new_model_results.csv")
adult <- adult[adult$method == "IPCW" &
                 adult$endpoint == "child_phenotype_adult_two_recent_windows" &
                 adult$term %in% c("phenotypeisolated", "phenotyperepeated") &
                 adult$cohort %in% c("NCDS", "BCS70"), ]
stopifnot(nrow(adult) == 4)
write_source(adult, "Figure_3B_source.csv")
adult$comparison <- factor(adult$term,
                           levels = c("phenotyperepeated", "phenotypeisolated"),
                           labels = c("2 vs 0", "1 vs 0"))
adult$cohort <- factor(adult$cohort, levels = c("NCDS", "BCS70"))
p3b <- ggplot(adult, aes(RR, comparison)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey65", linewidth = .35) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                width = .12, linewidth = .55, colour = blue) +
  geom_point(size = 2.2, colour = blue) +
  facet_grid(cohort ~ ., switch = "y") +
  scale_x_log10(limits = c(1, 8), breaks = c(1, 2, 4, 8)) +
  labs(x = "Adult repeated-report RR (95% CI)", y = "Childhood positive waves") +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0),
        panel.spacing.y = grid::unit(5, "mm"))
fig3 <- (p3a | p3b) + plot_layout(widths = c(1, 1.07)) +
  plot_annotation(tag_levels = "A")
save_figure(fig3, "Figure_3_childhood_to_respiratory_burden", 183, 100)

# Figure 4A: conditional weighted category distribution in those with all six
# adult outcome waves. Figure 4B: FEV1/FVC contrasts only, as captioned.
burden <- read_table(v3, "child_to_adult_burden_table.csv")
burden <- burden[burden$definition == "two_recent_windows" &
                   burden$cohort %in% c("NCDS", "BCS70"), ]
stopifnot(nrow(burden) == 18)
write_source(burden, "Figure_4A_source.csv")
burden$cohort <- factor(burden$cohort, levels = c("NCDS", "BCS70"))
burden$phenotype <- factor(burden$phenotype,
                           levels = c("never", "isolated", "repeated"),
                           labels = c("0", "1", "2"))
burden$adult_category <- factor(burden$adult_category,
                                levels = c("repeated", "isolated", "none"),
                                labels = c("\u22652", "1", "0"))
p4a <- ggplot(burden, aes(phenotype, weighted_probability, fill = adult_category)) +
  geom_col(width = .68) +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c("\u22652" = blue, "1" = "#B8CBD7", "0" = "#E5EBEF"),
                    breaks = c("0", "1", "\u22652")) +
  scale_y_continuous(labels = label_percent(), breaks = seq(0, 1, .25), expand = c(0, 0)) +
  labs(x = "Childhood positive waves", y = "Weighted adult distribution",
       fill = "Positive adult waves") +
  theme(legend.position = "bottom", legend.key.size = grid::unit(3, "mm"))

lung <- read_table(sev, "07_NCDS_lung_function_models.csv")
lung <- lung[lung$endpoint == "ratio" &
               lung$analysis == "strict_QC_adjusted_linear" &
               lung$term %in% c("child_groupisolated", "child_grouprepeated"), ]
stopifnot(nrow(lung) == 2)
write_source(lung, "Figure_4B_source.csv")
lung$comparison <- factor(lung$term,
                          levels = c("child_grouprepeated", "child_groupisolated"),
                          labels = c("2 vs 0", "1 vs 0"))
p4b <- ggplot(lung, aes(estimate, comparison)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey65", linewidth = .35) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                width = .12, linewidth = .55, colour = blue) +
  geom_point(size = 2.4, colour = blue) +
  scale_x_continuous(limits = c(-.12, .01), breaks = c(-.10, -.05, 0)) +
  labs(x = "Adjusted FEV1/FVC difference (95% CI)",
       y = "Childhood positive waves", title = "NCDS, age 44\u201346")
fig4 <- (p4a | p4b) + plot_layout(widths = c(1.35, .85)) +
  plot_annotation(tag_levels = "A")
save_figure(fig4, "Figure_4_adult_burden_and_lung_function", 183, 96)

writeLines(c(
  "Figure 1 copied without data changes from the existing V3 manuscript figure.",
  "Figure 2: refit of the published MCS survey-weighted, family-clustered stage-specific-wave model; 95% Wald CIs from its full covariance matrix; eight source point estimates verified against stage_interactions.csv.",
  "Figure 3A: MCS standardized age-11 high-symptom-proxy risks and model-conditional 95% CIs; Figure 3B: NCDS and BCS70 IPCW average repeated adult reporting RRs and 95% CIs.",
  "Figure 4A: weighted adult outcome distribution conditional on all six adult waves observed; Figure 4B: NCDS strict-QC age-44-46 FEV1/FVC coefficient contrasts and HC3 robust 95% CIs.",
  "All source CSVs and the R rendering script are in this folder. Figures 2-4 are exported as PDF, SVG, PNG 500 dpi, and TIFF 600 dpi."
), file.path(out, "figure_methods_and_sources.txt"))
cat("Manuscript figures 1-4 ready in:", out, "\n")
