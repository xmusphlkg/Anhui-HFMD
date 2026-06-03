rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(MASS)
  library(readr)
  library(tidyr)
  library(purrr)
})

if (getRversion() >= '2.15.1') {
  utils::globalVariables(c(
    'virus', 'virus_std', 'datafiles_all', 'datafiles_lo', 'datafiles_virus',
    '发病日期', 'non_missing_typed', 'reported_cases', 'typing_fraction',
    'method', 'counterfactual', 'lower_95', 'upper_95', 'actual_cases',
    'predicted_cases', 'difference', 'n_rescaled', 'endpoint'
  ))
}

if (dir.exists('../outcome')) {
  setwd('../outcome')
} else if (dir.exists('outcome')) {
  setwd('outcome')
} else {
  stop('Cannot find outcome/ directory from current working directory.')
}

data_dir <- if (dir.exists('../data_public')) '../data_public' else '../data'
load(file.path(data_dir, 'datafiles_all.RData'))
load(file.path(data_dir, 'datafile_virus.RData'))
load(file.path(data_dir, 'datafile_lo.RData'))

# -----------------------------
# Helpers
# -----------------------------

normalize_virus <- function(x) {
  dplyr::case_when(
    x %in% c('Cox A16', 'CV-A16') ~ 'CV-A16',
    x %in% c('EV71') ~ 'EV71',
    is.na(x) ~ NA_character_,
    TRUE ~ 'Others'
  )
}

get_date_col <- function(df) {
  if ('date' %in% names(df)) return('date')
  if ('发病日期' %in% names(df)) return('发病日期')
  stop('No recognized date column in dataset.')
}

get_severe_col <- function(df) {
  if ('severe' %in% names(df)) return('severe')
  if ('重症患者' %in% names(df)) return('重症患者')
  stop('No recognized severe column in dataset.')
}

prepare_monthly_counts <- function(df, end_date, virus_type = NULL, start_date = as.Date('2008-01-01')) {
  date_col <- get_date_col(df)

  dat <- df %>%
    mutate(date = as.Date(.data[[date_col]])) %>%
    filter(date >= start_date, date < end_date)

  if (!is.null(virus_type)) {
    dat <- dat %>%
      mutate(virus_std = normalize_virus(virus)) %>%
      filter(virus_std == virus_type)
  }

  dat %>%
    mutate(date = as.Date(format(date, '%Y-%m-01'))) %>%
    count(date, name = 'n') %>%
    arrange(date) %>%
    mutate(
      year = year(date),
      month = factor(month(date), levels = 1:12),
      time_index = row_number()
    )
}

annual_typing_fraction <- function() {
  all_year <- datafiles_all %>%
    mutate(date = as.Date(`发病日期`), year = year(date)) %>%
    filter(year >= 2008, year <= 2023) %>%
    count(year, name = 'reported_cases')

  typed_year <- datafiles_virus %>%
    mutate(date = as.Date(date), year = year(date)) %>%
    filter(year >= 2008, year <= 2023) %>%
    group_by(year) %>%
    summarise(non_missing_typed = n(), .groups = 'drop')

  all_year %>%
    left_join(typed_year, by = 'year') %>%
    mutate(
      non_missing_typed = coalesce(non_missing_typed, 0L),
      typing_fraction = pmax(non_missing_typed / reported_cases, 1e-6)
    ) %>%
    dplyr::select(year, typing_fraction)
}

fit_nb_counterfactual <- function(monthly_df, intervention_date, label_method, offset_col = NULL, trend = c('linear', 'quadratic')) {
  trend <- match.arg(trend)
  data_pre <- monthly_df %>% filter(date < intervention_date)

  if (nrow(data_pre) < 24) {
    stop('Pre-intervention window is too short (<24 months).')
  }

  if (is.null(offset_col)) {
    model_formula <- if (trend == 'quadratic') {
      n ~ poly(time_index, 2, raw = TRUE) + month
    } else {
      n ~ time_index + month
    }
    fit <- MASS::glm.nb(model_formula, data = data_pre)
    pred <- predict(fit, newdata = monthly_df, type = 'link', se.fit = TRUE)
  } else {
    trend_term <- if (trend == 'quadratic') {
      'poly(time_index, 2, raw = TRUE)'
    } else {
      'time_index'
    }
    offset_formula <- as.formula(
      paste0('n ~ ', trend_term, ' + month + offset(log(', offset_col, '))')
    )
    fit <- MASS::glm.nb(
      offset_formula,
      data = data_pre
    )
    pred <- predict(fit, newdata = monthly_df, type = 'link', se.fit = TRUE)
  }

  fit_link <- as.numeric(pred$fit)
  se_link <- as.numeric(pred$se.fit)

  cf <- exp(fit_link)
  lcl <- exp(fit_link - 1.96 * se_link)
  ucl <- exp(fit_link + 1.96 * se_link)

  post <- monthly_df %>%
    mutate(
      method = label_method,
      counterfactual = if_else(date >= intervention_date, cf, NA_real_),
      lower_95 = if_else(date >= intervention_date, lcl, NA_real_),
      upper_95 = if_else(date >= intervention_date, ucl, NA_real_)
    )

  pearson_res <- residuals(fit, type = 'pearson')
  disp_ratio <- sum(pearson_res^2, na.rm = TRUE) / fit$df.residual
  lb <- Box.test(pearson_res, lag = 12, type = 'Ljung-Box')

  diag_row <- tibble(
    method = label_method,
    theta = fit$theta,
    overdispersion_ratio = disp_ratio,
    aic = AIC(fit),
    ljung_box_p = lb$p.value,
    n_train_months = nrow(data_pre)
  )

  list(post = post, diagnostics = diag_row)
}

fit_nb_counterfactual_window <- function(monthly_df, intervention_date, training_end_date, label_method, trend = c('linear', 'quadratic')) {
  trend <- match.arg(trend)
  data_pre <- monthly_df %>%
    filter(date < intervention_date, date < training_end_date)

  if (nrow(data_pre) < 24) {
    stop('Pre-intervention window is too short (<24 months).')
  }

  model_formula <- if (trend == 'quadratic') {
    n ~ poly(time_index, 2, raw = TRUE) + month
  } else {
    n ~ time_index + month
  }

  fit <- MASS::glm.nb(model_formula, data = data_pre)
  pred <- predict(fit, newdata = monthly_df, type = 'link', se.fit = TRUE)

  fit_link <- as.numeric(pred$fit)
  se_link <- as.numeric(pred$se.fit)

  cf <- exp(fit_link)
  lcl <- exp(fit_link - 1.96 * se_link)
  ucl <- exp(fit_link + 1.96 * se_link)

  post <- monthly_df %>%
    mutate(
      method = label_method,
      counterfactual = if_else(date >= intervention_date, cf, NA_real_),
      lower_95 = if_else(date >= intervention_date, lcl, NA_real_),
      upper_95 = if_else(date >= intervention_date, ucl, NA_real_)
    )

  pearson_res <- residuals(fit, type = 'pearson')
  disp_ratio <- sum(pearson_res^2, na.rm = TRUE) / fit$df.residual
  lb <- Box.test(pearson_res, lag = 12, type = 'Ljung-Box')

  diag_row <- tibble(
    method = label_method,
    theta = fit$theta,
    overdispersion_ratio = disp_ratio,
    aic = AIC(fit),
    ljung_box_p = lb$p.value,
    n_train_months = nrow(data_pre),
    training_end_date = training_end_date
  )

  list(post = post, diagnostics = diag_row)
}

summarize_post_year <- function(df, intervention_date) {
  df %>%
    filter(date >= intervention_date) %>%
    mutate(year = year(date)) %>%
    group_by(method, year) %>%
    summarise(
      actual_cases = sum(n, na.rm = TRUE),
      predicted_cases = sum(counterfactual, na.rm = TRUE),
      lower_95 = sum(lower_95, na.rm = TRUE),
      upper_95 = sum(upper_95, na.rm = TRUE),
      difference = actual_cases - predicted_cases,
      percent_change = if_else(predicted_cases > 0, (difference / predicted_cases) * 100, NA_real_),
      .groups = 'drop'
    )
}

fit_apc_glm <- function(df, family_type = c('poisson', 'nb')) {
  family_type <- match.arg(family_type)

  if (family_type == 'poisson') {
    fit <- glm(n ~ year_c, family = poisson(), data = df)
    vc <- vcov(fit)
    if (requireNamespace('sandwich', quietly = TRUE)) {
      vc <- sandwich::vcovHC(fit, type = 'HC0')
    }
    beta <- coef(fit)['year_c']
    se <- sqrt(diag(vc))['year_c']
    model_name <- 'Poisson'
    theta <- NA_real_
  } else {
    fit <- MASS::glm.nb(n ~ year_c, data = df)
    beta <- coef(fit)['year_c']
    se <- sqrt(diag(vcov(fit)))['year_c']
    model_name <- 'Negative binomial'
    theta <- fit$theta
  }

  tibble(
    model = model_name,
    beta_year = as.numeric(beta),
    se_year = as.numeric(se),
    apc_pct = (exp(beta) - 1) * 100,
    apc_lcl_pct = (exp(beta - 1.96 * se) - 1) * 100,
    apc_ucl_pct = (exp(beta + 1.96 * se) - 1) * 100,
    p_value = 2 * pnorm(abs(beta / se), lower.tail = FALSE),
    theta = theta
  )
}

run_typed_sensitivity <- function(endpoint_name, virus_type, scenario, end_date, intervention_date) {
  vaccine_scenario <- identical(scenario, 'EV71 vaccine impact')
  analysis_start <- if (vaccine_scenario) as.Date('2010-01-01') else as.Date('2008-01-01')
  trend <- if (vaccine_scenario) 'quadratic' else 'linear'

  monthly_raw <- prepare_monthly_counts(datafiles_virus, end_date = end_date, virus_type = virus_type, start_date = analysis_start)
  typing_frac <- annual_typing_fraction()

  monthly <- monthly_raw %>%
    left_join(typing_frac, by = 'year') %>%
    mutate(
      typing_fraction = coalesce(typing_fraction, 1e-6),
      n_rescaled = n / typing_fraction
    )

  baseline <- fit_nb_counterfactual(monthly, intervention_date, 'baseline_typed', trend = trend)

  offset_fit <- fit_nb_counterfactual(
    monthly,
    intervention_date,
    'typing_fraction_offset',
    offset_col = 'typing_fraction',
    trend = trend
  )

  monthly_rescaled <- monthly %>% mutate(n = n_rescaled)
  rescaled_fit <- fit_nb_counterfactual(monthly_rescaled, intervention_date, 'annual_rescaled_typed', trend = trend)

  all_post <- bind_rows(
    baseline$post,
    offset_fit$post,
    rescaled_fit$post
  )

  out_year <- summarize_post_year(all_post, intervention_date) %>%
    mutate(
      scenario = scenario,
      endpoint = endpoint_name
    ) %>%
    dplyr::select(scenario, endpoint, everything())

  out_diag <- bind_rows(
    baseline$diagnostics,
    offset_fit$diagnostics,
    rescaled_fit$diagnostics
  ) %>%
    mutate(
      scenario = scenario,
      endpoint = endpoint_name
    ) %>%
    dplyr::select(scenario, endpoint, everything())

  list(effect = out_year, diag = out_diag)
}

# -----------------------------
# 1) Typing-fraction sensitivity for typed ITS
# -----------------------------

sens_ev71_vax <- run_typed_sensitivity(
  endpoint_name = 'EV71 typed cases',
  virus_type = 'EV71',
  scenario = 'EV71 vaccine impact',
  end_date = as.Date('2020-01-01'),
  intervention_date = as.Date('2016-06-01')
)

sens_cva16_vax <- run_typed_sensitivity(
  endpoint_name = 'CV-A16 typed cases',
  virus_type = 'CV-A16',
  scenario = 'EV71 vaccine impact',
  end_date = as.Date('2020-01-01'),
  intervention_date = as.Date('2016-06-01')
)

sens_other_vax <- run_typed_sensitivity(
  endpoint_name = 'Other typed cases',
  virus_type = 'Others',
  scenario = 'EV71 vaccine impact',
  end_date = as.Date('2020-01-01'),
  intervention_date = as.Date('2016-06-01')
)

sens_ev71_covid <- run_typed_sensitivity(
  endpoint_name = 'EV71 typed cases',
  virus_type = 'EV71',
  scenario = 'COVID-19 impact',
  end_date = as.Date('2024-01-01'),
  intervention_date = as.Date('2020-01-01')
)

sens_cva16_covid <- run_typed_sensitivity(
  endpoint_name = 'CV-A16 typed cases',
  virus_type = 'CV-A16',
  scenario = 'COVID-19 impact',
  end_date = as.Date('2024-01-01'),
  intervention_date = as.Date('2020-01-01')
)

sens_other_covid <- run_typed_sensitivity(
  endpoint_name = 'Other typed cases',
  virus_type = 'Others',
  scenario = 'COVID-19 impact',
  end_date = as.Date('2024-01-01'),
  intervention_date = as.Date('2020-01-01')
)

its_sensitivity_effect <- bind_rows(
  sens_ev71_vax$effect,
  sens_cva16_vax$effect,
  sens_other_vax$effect,
  sens_ev71_covid$effect,
  sens_cva16_covid$effect,
  sens_other_covid$effect
)

its_sensitivity_diag <- bind_rows(
  sens_ev71_vax$diag,
  sens_cva16_vax$diag,
  sens_other_vax$diag,
  sens_ev71_covid$diag,
  sens_cva16_covid$diag,
  sens_other_covid$diag
)

# -----------------------------
# 2) Severity-driven sampling bias check
# -----------------------------

sev_col <- get_severe_col(datafiles_lo)

typed <- datafiles_lo %>%
  mutate(
    date = as.Date(date),
    year = year(date),
    virus_std = normalize_virus(virus),
    severe_std = case_when(
      .data[[sev_col]] %in% c('Y', '是', 'Yes', '1', 1) ~ 'Severe',
      .data[[sev_col]] %in% c('N', '否', 'No', '0', 0) ~ 'Mild',
      TRUE ~ NA_character_
    )
  ) %>%
  filter(year >= 2008, year <= 2023, !is.na(virus_std), !is.na(severe_std))

sev_year <- typed %>%
  group_by(year, severe_std) %>%
  summarise(
    n_typed = n(),
    n_ev71 = sum(virus_std == 'EV71', na.rm = TRUE),
    ev71_prop = n_ev71 / n_typed,
    odds_ev71 = n_ev71 / pmax(n_typed - n_ev71, 1),
    .groups = 'drop'
  )

sev_or <- sev_year %>%
  dplyr::select(year, severe_std, odds_ev71) %>%
  pivot_wider(names_from = severe_std, values_from = odds_ev71) %>%
  mutate(
    severe_vs_mild_or = Severe / Mild
  ) %>%
  dplyr::select(year, severe_vs_mild_or)

typed2 <- typed %>%
  mutate(
    ev71_binary = if_else(virus_std == 'EV71', 1L, 0L),
    year_c = year - 2008,
    severe_binary = if_else(severe_std == 'Severe', 1L, 0L)
  )

fit_bias <- glm(ev71_binary ~ year_c * severe_binary, family = binomial(), data = typed2)
bias_coef <- summary(fit_bias)$coefficients

bias_coef_df <- as.data.frame(bias_coef)
bias_coef_df$term <- rownames(bias_coef_df)
rownames(bias_coef_df) <- NULL
names(bias_coef_df) <- c('estimate', 'std_error', 'z_value', 'p_value', 'term')

# -----------------------------
# 3) County-year typed sample size for dominance confidence
# -----------------------------

county_year <- datafiles_lo %>%
  mutate(
    date = as.Date(date),
    year = year(date),
    virus_std = normalize_virus(virus)
  ) %>%
  filter(year >= 2008, year <= 2023, !is.na(virus_std))

county_virus <- county_year %>%
  group_by(name, year, virus_std) %>%
  summarise(typed_cases = n(), .groups = 'drop')

county_rank <- county_virus %>%
  group_by(name, year) %>%
  arrange(desc(typed_cases), .by_group = TRUE) %>%
  mutate(rk = row_number()) %>%
  ungroup()

county_conf <- county_rank %>%
  group_by(name, year) %>%
  summarise(
    typed_total = sum(typed_cases),
    top_virus = if_else(sum(typed_cases == max(typed_cases)) > 1, 'Mixed', first(virus_std[typed_cases == max(typed_cases)])),
    top_cases = max(typed_cases),
    second_cases = if_else(any(rk == 2), typed_cases[rk == 2][1], 0L),
    top_share = top_cases / typed_total,
    dominance_gap = top_cases - second_cases,
    .groups = 'drop'
  )

# -----------------------------
# 4) GLM-based trend sensitivity for typed serotypes
# -----------------------------

annual_serotype <- datafiles_virus %>%
  mutate(
    date = as.Date(date),
    year = year(date),
    virus_std = normalize_virus(virus)
  ) %>%
  filter(year >= 2008, year <= 2023, !is.na(virus_std)) %>%
  count(year, virus_std, name = 'n') %>%
  mutate(year_c = year - min(year))

glm_apc_sensitivity <- annual_serotype %>%
  group_by(virus_std) %>%
  group_modify(~ bind_rows(
    fit_apc_glm(.x, 'poisson'),
    fit_apc_glm(.x, 'nb')
  )) %>%
  ungroup() %>%
  rename(serotype = virus_std)

# -----------------------------
# 5) COVID-period sensitivity using a pre-vaccine baseline
# -----------------------------

run_covid_pre_vaccine_sensitivity <- function(endpoint_name, virus_type = NULL) {
  monthly_raw <- prepare_monthly_counts(
    if (is.null(virus_type)) datafiles_all else datafiles_virus,
    end_date = as.Date('2024-01-01'),
    virus_type = virus_type,
    start_date = as.Date('2008-01-01')
  )

  fit <- fit_nb_counterfactual_window(
    monthly_df = monthly_raw,
    intervention_date = as.Date('2020-01-01'),
    training_end_date = as.Date('2016-01-01'),
    label_method = 'pre_vaccine_baseline_2008_2015',
    trend = 'linear'
  )

  effect <- summarize_post_year(fit$post, as.Date('2020-01-01')) %>%
    mutate(
      scenario = 'COVID-19 sensitivity (pre-vaccine baseline)',
      endpoint = endpoint_name
    ) %>%
    dplyr::select(scenario, endpoint, method, year, actual_cases, predicted_cases, lower_95, upper_95, difference, percent_change)

  list(effect = effect, diag = fit$diagnostics)
}

covid_all_pre_vax <- run_covid_pre_vaccine_sensitivity(
  endpoint_name = 'All reported HFMD cases',
  virus_type = NULL
)

covid_ev71_pre_vax <- run_covid_pre_vaccine_sensitivity(
  endpoint_name = 'Typed EV71 cases',
  virus_type = 'EV71'
)

covid_cva16_pre_vax <- run_covid_pre_vaccine_sensitivity(
  endpoint_name = 'Typed CV-A16 cases',
  virus_type = 'CV-A16'
)

covid_other_pre_vax <- run_covid_pre_vaccine_sensitivity(
  endpoint_name = 'Typed other-enterovirus cases',
  virus_type = 'Others'
)

covid_pre_vax_effect <- bind_rows(
  covid_all_pre_vax$effect,
  covid_ev71_pre_vax$effect,
  covid_cva16_pre_vax$effect,
  covid_other_pre_vax$effect
)

# -----------------------------
# Write outputs
# -----------------------------

dir.create('supplementary_tables', showWarnings = FALSE, recursive = TRUE)

write_csv(its_sensitivity_effect, 'supplementary_tables/its_typing_fraction_sensitivity.csv')
write_csv(its_sensitivity_diag, 'supplementary_tables/its_nb_diagnostics_extended.csv')
write_csv(sev_year, 'supplementary_tables/severe_mild_ev71_yearly.csv')
write_csv(sev_or, 'supplementary_tables/severe_vs_mild_ev71_or_yearly.csv')
write_csv(bias_coef_df, 'supplementary_tables/severe_sampling_bias_logit.csv')
write_csv(county_conf, 'supplementary_tables/supplementary_data_D7_county_map_masking.csv')
write_csv(glm_apc_sensitivity, 'supplementary_tables/joinpoint_glm_apc_sensitivity.csv')
write_csv(covid_pre_vax_effect, 'supplementary_tables/covid_pre_vaccine_baseline_sensitivity.csv')

cat('Wrote supplementary_tables/its_typing_fraction_sensitivity.csv\n')
cat('Wrote supplementary_tables/its_nb_diagnostics_extended.csv\n')
cat('Wrote supplementary_tables/severe_mild_ev71_yearly.csv\n')
cat('Wrote supplementary_tables/severe_vs_mild_ev71_or_yearly.csv\n')
cat('Wrote supplementary_tables/severe_sampling_bias_logit.csv\n')
cat('Wrote supplementary_tables/supplementary_data_D7_county_map_masking.csv\n')
cat('Wrote supplementary_tables/joinpoint_glm_apc_sensitivity.csv\n')
cat('Wrote supplementary_tables/covid_pre_vaccine_baseline_sensitivity.csv\n')
