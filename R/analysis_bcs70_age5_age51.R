options(stringsAsFactors = FALSE)
Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8")
suppressPackageStartupMessages({
  library(haven)
  library(geepack)
  library(sandwich)
})

root <- "E:/公共数据库/P1"
out <- file.path(root, "bcs70_age5_age51_extension_20261001")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
full <- readRDS("E:/公共数据库/Four_Cohorts_Full_RDS_20260925/BCS70_full.rds")
get_source <- function(name) {
  idx <- which(basename(full$source_inventory$source_entry) == name &
                 full$source_inventory$format == ".dta")
  stopifnot(length(idx) == 1L)
  full$datasets[[idx]]
}
age5 <- get_source("f699b.dta")
age51 <- get_source("bcs11_age51_main.dta")
rm(full)
gc()
base <- readRDS(file.path(root, "life_course_v3_20260928", "v3_cohort_datasets.rds"))$BCS70
stopifnot(!anyDuplicated(base$id), !anyDuplicated(age5$bcsid),
          !anyDuplicated(age51$bcsid))
num <- function(x) suppressWarnings(as.numeric(unclass(x)))
match5 <- match(base$id, as.character(age5$bcsid))
match51 <- match(base$id, as.character(age51$bcsid))
v5 <- num(age5$e087)[match5]
v51 <- num(age51$b11khlpb01)[match51]
v51doctor <- num(age51$b11asthd)[match51]
base$ch5ever <- ifelse(v5 %in% c(1, 3), 1,
                        ifelse(v5 == 2, 0, NA_real_))
base$y51 <- ifelse(v51 == 1, 1, ifelse(v51 == 0, 0, NA_real_))
base$y51doctor <- ifelse(base$y51 == 0, 0,
                         ifelse(base$y51 == 1 & v51doctor == 1, 1,
                                ifelse(base$y51 == 1 & v51doctor == 2, 0, NA_real_)))
base$ch5_admission <- {
  z <- num(age5$e092)[match5]
  ifelse(base$ch5ever == 0, 0,
         ifelse(base$ch5ever == 1 & z >= 0, as.integer(z > 0), NA_real_))
}
common <- c("male", "low_birthweight", "ses", "pregnancy_continued_smoke")
ages6 <- c(26, 30, 34, 38, 42, 46)
ages7 <- c(ages6, 51)

count_status <- data.frame(
  measure = c("birth_cohort", "age5_file_rows", "age5_matched",
              "age5_ever_wheeze_valid", "age5_ever_wheeze_yes",
              "age5_wheeze_hospital_valid", "age5_wheeze_hospital_yes",
              "age51_file_rows", "age51_matched", "age51_item_asked_valid",
              "age51_asthma_wheezy_bronchitis_yes", "age51_doctor_confirmed_yes",
              "age51_not_asked_fieldwork_stage"),
  n = c(nrow(base), nrow(age5), sum(!is.na(match5)),
        sum(!is.na(base$ch5ever)), sum(base$ch5ever == 1, na.rm = TRUE),
        sum(!is.na(base$ch5_admission)), sum(base$ch5_admission == 1, na.rm = TRUE),
        nrow(age51), sum(!is.na(match51)), sum(!is.na(base$y51)),
        sum(base$y51 == 1, na.rm = TRUE), sum(base$y51doctor == 1, na.rm = TRUE),
        sum(v51 == -3, na.rm = TRUE))
)
write.csv(count_status, file.path(out, "source_linkage_and_outcome_counts.csv"), row.names = FALSE)

ipcw_long <- function(d, ages, vars) {
  req <- vars
  d <- d[complete.cases(d[req]), ]
  d$base <- d$child_selection_weight
  active <- seq_len(nrow(d))
  cw <- rep(1, nrow(d))
  last <- hist <- rep(0, nrow(d))
  waves <- list()
  all_waves <- list()
  for (j in seq_along(ages)) {
    age <- ages[j]
    y <- d[[paste0("y", age)]]
    active <- active[!is.na(active)]
    r <- as.integer(!is.na(y[active]))
    q <- d[active, ]
    q$r <- r
    q$last <- last[active]
    q$history <- hist[active] / max(j - 1, 1)
    model <- glm(reformulate(c(vars, if (j > 1) c("last", "history")), "r"),
                 q, family = binomial())
    pd <- pmin(.98, pmax(.02, predict(model, type = "response")))
    pn <- rep(mean(r), length(r))
    cw[active] <- cw[active] * pn / pd
    keep <- active[r == 1]
    z <- d[keep, ]
    z$age <- age
    z$y <- y[keep]
    z$previous <- last[keep]
    z$ipcw_raw <- cw[keep]
    cuts <- quantile(z$ipcw_raw, c(.01, .99))
    z$ipcw <- pmin(cuts[2], pmax(cuts[1], z$ipcw_raw))
    z$w <- z$base * z$ipcw
    waves[[j]] <- z
    last[keep] <- y[keep]
    hist[keep] <- hist[keep] + y[keep]
    active <- keep
    u <- d[!is.na(y), ]
    u$age <- age
    u$y <- y[!is.na(y)]
    u$w <- u$base
    all_waves[[j]] <- u
  }
  mono <- do.call(rbind, waves)
  mono$wave <- factor(mono$age, levels = ages)
  all <- do.call(rbind, all_waves)
  all$wave <- factor(all$age, levels = ages)
  list(monotone = mono, all = all, baseline = d)
}

fit_gee <- function(z, terms) {
  z <- droplevels(z[order(z$id, z$age), ])
  z$gid <- match(z$id, unique(z$id))
  f <- geeglm(reformulate(c("wave", common, "phenotype"), "y"),
              data = z, id = gid, weights = w,
              family = poisson("log"), corstr = "independence")
  beta <- coef(f)
  V <- vcov(f)
  s <- sqrt(diag(V))
  out <- data.frame(term = terms, RR = exp(beta[terms]),
                    lower = exp(beta[terms] - 1.96 * s[terms]),
                    upper = exp(beta[terms] + 1.96 * s[terms]),
                    p = 2 * pnorm(-abs(beta[terms] / s[terms])),
                    row.names = NULL)
  delta <- beta["phenotyperepeated"] - beta["phenotypeisolated"]
  se <- sqrt(V["phenotyperepeated", "phenotyperepeated"] +
               V["phenotypeisolated", "phenotypeisolated"] -
               2 * V["phenotyperepeated", "phenotypeisolated"])
  contrast <- data.frame(term = "repeated_vs_isolated", RR = exp(delta),
                         lower = exp(delta - 1.96 * se),
                         upper = exp(delta + 1.96 * se),
                         p = 2 * pnorm(-abs(delta / se)))
  rbind(out, contrast)
}

results <- list()
flows <- list()
weight_diagnostics <- list()
for (definition in c("two_recent_windows", "three_mixed_windows")) {
  chvars <- if (definition == "two_recent_windows") c("chfirst", "ch16")
            else c("ch5ever", "chfirst", "ch16")
  b <- base[complete.cases(base[common]), ]
  b$child_complete <- as.integer(complete.cases(b[chvars]))
  selection <- glm(reformulate(common, "child_complete"), b, family = binomial())
  ps <- pmin(.98, pmax(.02, predict(selection, type = "response")))
  b$child_selection_weight <- mean(b$child_complete) / ps
  b <- b[b$child_complete == 1, ]
  cuts <- quantile(b$child_selection_weight, c(.01, .99))
  b$child_selection_weight <- pmin(cuts[2], pmax(cuts[1], b$child_selection_weight))
  n_positive <- rowSums(b[chvars])
  b$phenotype <- factor(ifelse(n_positive == 0, "never",
                               ifelse(n_positive == 1, "isolated", "repeated")),
                        levels = c("never", "isolated", "repeated"))
  for (years in c(6L, 7L)) {
    ages <- if (years == 6L) ages6 else ages7
    l <- ipcw_long(b, ages, c(common, "phenotype"))
    for (age in ages) {
      observed <- l$monotone[l$monotone$age == age, ]
      w <- observed$ipcw
      weight_diagnostics[[paste(definition, years, age)]] <- data.frame(
        definition = definition, adult_waves = years, age = age,
        observed_monotone = nrow(observed),
        ipcw_min = min(w), ipcw_p01 = unname(quantile(w, .01)),
        ipcw_median = median(w), ipcw_p99 = unname(quantile(w, .99)),
        ipcw_max = max(w),
        effective_sample_size = sum(observed$w)^2 / sum(observed$w^2)
      )
    }
    for (method in c("IPCW", "all_available")) {
      z <- if (method == "IPCW") l$monotone else l$all
      estimates <- fit_gee(z, c("phenotypeisolated", "phenotyperepeated"))
      estimates$definition <- definition
      estimates$adult_waves <- years
      estimates$method <- method
      estimates$people <- length(unique(z$id))
      estimates$person_waves <- nrow(z)
      results[[paste(definition, years, method)]] <- estimates
    }
    flows[[paste(definition, years)]] <- data.frame(
      definition = definition, adult_waves = years,
      baseline_core_complete = sum(complete.cases(base[common])),
      childhood_complete = nrow(b),
      childhood_never = sum(b$phenotype == "never"),
      childhood_isolated = sum(b$phenotype == "isolated"),
      childhood_repeated = sum(b$phenotype == "repeated"),
      age51_observed_in_childhood_complete = sum(!is.na(b$y51)),
      any_adult_observation = length(unique(l$all$id)),
      monotone_all_waves = sum(l$monotone$age == max(ages)),
      stringsAsFactors = FALSE
    )
    cat("DONE", definition, years, "waves:", nrow(b), "childhood complete;",
        sum(l$monotone$age == max(ages)), "all-wave responders\n")
  }
}
result_table <- do.call(rbind, results)
flow_table <- do.call(rbind, flows)
write.csv(result_table, file.path(out, "adult_GEE_6_vs_7_waves.csv"), row.names = FALSE)
write.csv(flow_table, file.path(out, "extension_sample_flow.csv"), row.names = FALSE)
write.csv(do.call(rbind, weight_diagnostics),
          file.path(out, "IPCW_weight_diagnostics.csv"), row.names = FALSE)

# Verify the original six-wave, two-childhood-window primary model reproduces
# the published V3 point estimates before interpreting the extension.
published <- read.csv(file.path(root, "life_course_v3_20260928", "new_model_results.csv"))
published <- published[published$cohort == "BCS70" & published$method == "IPCW" &
                         published$endpoint == "child_phenotype_adult_two_recent_windows" &
                         published$term %in% c("phenotypeisolated", "phenotyperepeated"), ]
reproduced <- result_table[result_table$definition == "two_recent_windows" &
                             result_table$adult_waves == 6 &
                             result_table$method == "IPCW" &
                             result_table$term %in% published$term, ]
published <- published[match(reproduced$term, published$term), ]
stopifnot(nrow(reproduced) == 2L,
          all(abs(reproduced$RR - published$RR) < 1e-7))
writeLines(c("Original BCS70 six-wave, two-childhood-window IPCW model reproduced exactly.",
             "Age-5 wheeze is an ever-history variable; the three-window phenotype mixes recall periods and is secondary.",
             "Age-51 asthma/wheezy bronchitis is a since-last-interview report; 677 fieldwork-stage not-asked responses remain missing.",
             "Adult seven-wave GEE extends the reporting-burden estimand through age 51; all-available models are response-pattern sensitivity analyses.",
             "Neither age-5 wheeze-related admission nor age-51 doctor diagnosis establishes clinical severe asthma."),
           file.path(out, "analysis_notes.txt"))
cat("BCS70 age-5 / age-51 extension finished\n")
