# British birth cohorts: asthma life-course analysis code

This repository contains the R scripts used for the local analyses of NCDS, BCS70, MCS, and Next Steps. It includes source audits, cohort construction, harmonisation, longitudinal trajectory and RCS analyses, parent-behaviour models, the formal life-course analyses, V2/V3 sensitivity analyses, quality checks, and figure-generation code.

**Only code is published.** Cohort microdata, derived RDS files, tables, figures, reports, logs, and fitted model objects remain outside this repository.

## Script groups

- `build_full_cohort_rds.R`, `audit_*`, `inspect_*`, `trajectory_*`: source construction and variable audits.
- `augment_bcs70_full_rds.R`: audited addition of the UKDA 2699 age-5 and UKDA 9347 age-51 files to the existing BCS70 full RDS, with a preserved backup.
- `four_cohort_trajectory_rcs.R`, `mcs_*`, `plot_four_cohort_rcs.R`, `plot_mcs_groups.R`: exploratory trajectories and restricted cubic splines.
- `formal_*`, `plot_formal_results.R`: formal MCS, adult-cohort, and Next Steps analyses.
- `parent_behavior_*`, `plot_parent_behavior.R`: parent-behaviour harmonisation and association models.
- `v2_*`: transition and adult repeated-burden models with weighting and checks.
- `v3_child_audit.R`, `v3_models.R`, `v3_mi.R`, `v3_finalize_pool.R`, `v3_quality.R`, `v3_figures.R`: childhood phenotype, adult burden, equivalence, stage interaction, marginal risk, MI + IPCW, and figures.
- `analysis_severe_asthma_relevant.R`, `plot_severe_asthma_relevant.R`: high-symptom, hospitalisation and lung-function analyses and their source figures.
- `render_manuscript_figures.R`: manuscript Figures 1–4, including newly calculated stage-specific confidence intervals and composite Figures 3–4.

## Reproducing the analysis

The scripts are snapshots of the code executed locally. They contain absolute paths to the original **local** RDS inputs and output staging folders. Change those paths to directories containing your authorised cohort files before running. The repository does not provide licensed cohort data or reproduce the analysis without those inputs. `v3_mi.R` reads `v3_models.R` from the current working directory, so run it from `R/` or adjust that reference.

Typical V3 order after the upstream source and V2 datasets have been prepared: `v3_child_audit.R` → `v3_models.R` → `v3_mi.R` → `v3_quality.R` → `v3_figures.R`. `v3_finalize_pool.R` is an idempotent correction/check for pooled MI diagnostics.

The scripts use R packages including `survey`, `sandwich`, `geepack`, `mice`, and `ggplot2`; individual files declare their other dependencies.
