options(stringsAsFactors = FALSE)
Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scales)
})
root <- "E:/公共数据库/P1"
out <- file.path(root, "manuscript_figures_20261001")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
previous <- file.path(root, "manuscript_figures_20260929")
for (figure in c("Figure_1_cohort_equivalence", "Figure_2_transition_stages",
                 "Figure_4_adult_burden_and_lung_function")) {
  for (ext in c("pdf", "svg", "png", "tiff")) {
    source <- file.path(previous, paste0(figure, ".", ext))
    if (file.exists(source)) {
      stopifnot(file.copy(source, file.path(out, basename(source)), overwrite = TRUE))
    }
  }
}
blue <- "#527C9C"
theme_set(theme_classic(base_size = 8.5, base_family = "Arial") +
            theme(axis.line = element_line(linewidth = .35),
                  axis.ticks = element_line(linewidth = .35),
                  strip.background = element_blank(),
                  strip.text = element_text(face = "bold"),
                  plot.tag = element_text(face = "bold", size = 11)))

risk <- read.csv(file.path(root, "severe_asthma_relevant_analysis_20260928",
                           "MCS_prior_history_standardized_risks.csv"))
contrast <- read.csv(file.path(root, "severe_asthma_relevant_analysis_20260928",
                               "04_childhood_burden_to_severe_wheeze.csv"))
contrast <- contrast[contrast$term == "repeated_vs_isolated", ]
risk$prior_group <- factor(risk$prior_group,
                           levels = c("never", "isolated", "repeated"),
                           labels = c("0", "1", "≥2"))
p_a <- ggplot(risk, aes(prior_group, adjusted_risk)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .12,
                linewidth = .6, colour = blue) +
  geom_point(size = 2.6, colour = blue) +
  scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0, .26),
                     breaks = seq(0, .25, .05)) +
  labs(x = "Positive wheeze waves at ages 3, 5, 7",
       y = "Age-11 severe-wheeze proxy risk (95% CI)", title = "MCS") +
  annotate("text", x = 1.05, y = .245, hjust = 0, size = 2.55,
           label = sprintf("Repeated vs isolated: RR %.2f (%.2f–%.2f)",
                           contrast$estimate, contrast$lower, contrast$upper))

old <- read.csv(file.path(root, "life_course_v3_20260928", "new_model_results.csv"))
ncds <- old[old$cohort == "NCDS" & old$method == "IPCW" &
              old$endpoint == "child_phenotype_adult_two_recent_windows" &
              old$term %in% c("phenotypeisolated", "phenotyperepeated"),
            c("cohort", "term", "RR", "lower", "upper")]
extension <- read.csv(file.path(root, "bcs70_age5_age51_extension_20261001",
                                "adult_GEE_6_vs_7_waves.csv"))
bcs <- extension[extension$definition == "two_recent_windows" &
                   extension$adult_waves == 7 & extension$method == "IPCW" &
                   extension$term %in% c("phenotypeisolated", "phenotyperepeated"),
                 c("term", "RR", "lower", "upper")]
bcs$cohort <- "BCS70"
adult <- rbind(ncds, bcs[, names(ncds)])
stopifnot(nrow(adult) == 4L)
write.csv(adult, file.path(out, "Figure_3B_source.csv"), row.names = FALSE)
adult$comparison <- factor(adult$term,
                           levels = c("phenotyperepeated", "phenotypeisolated"),
                           labels = c("2 vs 0", "1 vs 0"))
adult$cohort <- factor(adult$cohort, levels = c("NCDS", "BCS70"),
                       labels = c("NCDS, 23–62 y", "BCS70, 26–51 y"))
p_b <- ggplot(adult, aes(RR, comparison)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey65", linewidth = .35) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y",
                width = .12, linewidth = .55, colour = blue) +
  geom_point(size = 2.2, colour = blue) +
  facet_grid(cohort ~ ., switch = "y") +
  scale_x_log10(limits = c(1, 8), breaks = c(1, 2, 4, 8)) +
  labs(x = "Adult repeated-report RR (95% CI)", y = "Childhood positive waves") +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0),
        panel.spacing.y = grid::unit(5, "mm"))

figure <- (p_a | p_b) + plot_layout(widths = c(1, 1.07)) +
  plot_annotation(tag_levels = "A")
base <- file.path(out, "Figure_3_childhood_to_respiratory_burden_age51")
ggsave(paste0(base, ".pdf"), figure, width = 183/25.4, height = 100/25.4,
       device = cairo_pdf, bg = "white")
ggsave(paste0(base, ".svg"), figure, width = 183/25.4, height = 100/25.4,
       device = svglite::svglite, bg = "white")
ggsave(paste0(base, ".png"), figure, width = 183/25.4, height = 100/25.4,
       dpi = 500, device = ragg::agg_png, bg = "white")
ggsave(paste0(base, ".tiff"), figure, width = 183/25.4, height = 100/25.4,
       dpi = 600, device = ragg::agg_tiff, compression = "lzw", bg = "white")
cat("Updated Figure 3 saved in", out, "\n")
