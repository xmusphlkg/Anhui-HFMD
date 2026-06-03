rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(MASS)
  library(patchwork)
  library(readr)
  library(scales)
})

if (dir.exists("../outcome")) {
  setwd("../outcome")
} else if (dir.exists("outcome")) {
  setwd("outcome")
} else {
  stop("Cannot find outcome/ directory from current working directory.")
}

data_dir <- if (dir.exists("../data_public")) "../data_public" else "../data"
load(file.path(data_dir, "datafile_virus.RData"))
load(file.path(data_dir, "datafiles_all.RData"))

dir.create("supplementary_tables", showWarnings = FALSE, recursive = TRUE)

normalize_virus <- function(x) {
  dplyr::case_when(
    x %in% c("CV-A16", "Cox A16") ~ "CV-A16",
    x == "EV71" ~ "EV71",
    is.na(x) ~ NA_character_,
    TRUE ~ "Others"
  )
}

get_date_col <- function(df) {
  if ("date" %in% names(df)) {
    return("date")
  }
  if ("发病日期" %in% names(df)) {
    return("发病日期")
  }
  stop("No recognized date column in dataset.")
}

prepare_monthly_counts <- function(df, start_date, end_date, virus_type = NULL) {
  date_col <- get_date_col(df)

  dat <- df %>%
    mutate(date = as.Date(.data[[date_col]])) %>%
    filter(date >= start_date, date < end_date)

  if (!is.null(virus_type)) {
    dat <- dat %>%
      mutate(virus_std = normalize_virus(virus)) %>%
      filter(virus_std == virus_type)
  }

  monthly_seq <- seq.Date(
    from = as.Date(format(start_date, "%Y-%m-01")),
    to = as.Date(format(end_date - days(1), "%Y-%m-01")),
    by = "month"
  )

  dat %>%
    mutate(date = as.Date(format(date, "%Y-%m-01"))) %>%
    count(date, name = "n") %>%
    right_join(tibble(date = monthly_seq), by = "date") %>%
    mutate(n = coalesce(n, 0L)) %>%
    arrange(date) %>%
    mutate(
      month = factor(month(date), levels = 1:12),
      time_index = row_number()
    )
}

fit_nb_counterfactual <- function(monthly_df, intervention_date, trend = c("linear", "quadratic")) {
  trend <- match.arg(trend)

  data_pre <- monthly_df %>%
    filter(date < intervention_date)

  if (nrow(data_pre) < 24) {
    stop("Pre-intervention window is too short for a stable seasonal model.")
  }

  model_formula <- if (trend == "quadratic") {
    n ~ poly(time_index, 2, raw = TRUE) + month
  } else {
    n ~ time_index + month
  }

  fit <- MASS::glm.nb(model_formula, data = data_pre)
  pred <- predict(fit, newdata = monthly_df, type = "link", se.fit = TRUE)

  fit_link <- as.numeric(pred$fit)
  se_link <- as.numeric(pred$se.fit)

  set.seed(20260424)
  n_sim <- 1500L
  eta_sim <- matrix(
    rnorm(length(fit_link) * n_sim, mean = rep(fit_link, n_sim), sd = rep(se_link, n_sim)),
    nrow = length(fit_link),
    ncol = n_sim
  )
  mu_sim <- exp(eta_sim)
  sim_matrix <- matrix(
    rnbinom(length(mu_sim), mu = as.vector(mu_sim), size = fit$theta),
    nrow = length(fit_link),
    ncol = n_sim
  )

  out <- monthly_df %>%
    mutate(
      counterfactual = if_else(date >= intervention_date, exp(fit_link), NA_real_),
      lower_95 = if_else(date >= intervention_date, apply(sim_matrix, 1, quantile, probs = 0.025), NA_real_),
      upper_95 = if_else(date >= intervention_date, apply(sim_matrix, 1, quantile, probs = 0.975), NA_real_),
      difference = n - counterfactual,
      change = case_when(
        is.na(difference) ~ NA_character_,
        difference >= 0 ~ "Increase",
        TRUE ~ "Decrease"
      ),
      train_start = min(data_pre$date),
      train_end = max(data_pre$date),
      trend = trend
    )

  attr(out, "sim_matrix") <- sim_matrix
  out
}

summarize_post_period <- function(df, scenario, endpoint, model_spec, intervention_date) {
  start_year <- if (month(intervention_date) == 1) year(intervention_date) else year(intervention_date) + 1
  sim_matrix <- attr(df, "sim_matrix")

  post_index <- which(df$date >= intervention_date & year(df$date) >= start_year)
  year_index <- year(df$date[post_index])
  year_levels <- sort(unique(year_index))

  annual_pi <- lapply(year_levels, function(year_value) {
    row_ids <- post_index[year_index == year_value]
    sim_year <- colSums(sim_matrix[row_ids, , drop = FALSE])
    tibble(
      year = year_value,
      lower_95 = quantile(sim_year, probs = 0.025, na.rm = TRUE),
      upper_95 = quantile(sim_year, probs = 0.975, na.rm = TRUE)
    )
  }) %>%
    bind_rows()

  df %>%
    filter(date >= intervention_date) %>%
    mutate(year = year(date)) %>%
    filter(year >= start_year) %>%
    group_by(year) %>%
    summarise(
      actual_cases = sum(n, na.rm = TRUE),
      predicted_cases = sum(counterfactual, na.rm = TRUE),
      difference = actual_cases - predicted_cases,
      percent_change = if_else(predicted_cases > 0, difference / predicted_cases * 100, NA_real_),
      train_start = first(train_start),
      train_end = first(train_end),
      trend = first(trend),
      .groups = "drop"
    ) %>%
    left_join(annual_pi, by = "year") %>%
    mutate(
      scenario = scenario,
      endpoint = endpoint,
      model_spec = model_spec,
      intervention_date = intervention_date
    ) %>%
    dplyr::select(
      scenario,
      endpoint,
      model_spec,
      intervention_date,
      train_start,
      train_end,
      trend,
      year,
      actual_cases,
      predicted_cases,
      lower_95,
      upper_95,
      difference,
      percent_change
    )
}

plot_counterfactual <- function(df, intervention_date, plot_title_main, plot_title_diff, x_limits, y_limits = NULL) {
  if (!inherits(x_limits, "Date")) {
    x_limits <- as.Date(x_limits, origin = "1970-01-01")
  }

  plot_main <- ggplot(df, aes(x = date)) +
    geom_ribbon(
      data = df %>% filter(date >= intervention_date),
      aes(ymin = lower_95, ymax = upper_95),
      fill = "#F18B00FF",
      alpha = 0.18
    ) +
    geom_line(aes(y = n, color = "Observed"), linewidth = 0.5) +
    geom_line(
      data = df %>% filter(date >= intervention_date),
      aes(y = counterfactual, color = "Predicted"),
      linewidth = 0.6,
      linetype = "22"
    ) +
    geom_vline(xintercept = as.numeric(intervention_date), color = "black", linewidth = 0.4) +
    scale_color_manual(
      values = c("Observed" = "#3CC8C0FF", "Predicted" = "#F18B00FF"),
      limits = c("Observed", "Predicted"),
      drop = FALSE
    ) +
    scale_x_date(
      limits = x_limits,
      expand = c(0, 0),
      date_breaks = "2 years",
      date_labels = "%Y"
    ) +
    scale_y_continuous(
      labels = label_comma(),
      limits = y_limits,
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(
      title = plot_title_main,
      x = NULL,
      y = "Number of cases",
      color = NULL
    ) +
    theme_bw(base_family = "Times New Roman") +
    theme(
      panel.grid = element_blank(),
      text = element_text(color = "black"),
      axis.text = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      plot.title.position = "plot"
    )

  plot_diff <- ggplot(df %>% filter(date >= intervention_date), aes(x = date, y = difference, fill = change)) +
    geom_col(alpha = 0.45, show.legend = FALSE) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
    scale_fill_manual(values = c("Increase" = "#C0392B", "Decrease" = "#27AE60")) +
    scale_x_date(
      limits = x_limits,
      expand = c(0, 0),
      date_breaks = "2 years",
      date_labels = "%Y"
    ) +
    scale_y_continuous(labels = label_comma()) +
    labs(
      title = plot_title_diff,
      x = NULL,
      y = "Difference"
    ) +
    theme_bw(base_family = "Times New Roman") +
    theme(
      panel.grid = element_blank(),
      text = element_text(color = "black"),
      axis.text = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.title.position = "plot"
    )

  list(main = plot_main, diff = plot_diff)
}

run_counterfactual <- function(df, endpoint, scenario, start_date, end_date, intervention_date, trend, virus_type = NULL) {
  monthly <- prepare_monthly_counts(
    df = df,
    start_date = start_date,
    end_date = end_date,
    virus_type = virus_type
  )

  fit_nb_counterfactual(
    monthly_df = monthly,
    intervention_date = intervention_date,
    trend = trend
  )
}

intervention_ev71 <- as.Date("2016-06-01")
intervention_covid <- as.Date("2020-01-01")

vax_all_primary <- run_counterfactual(
  df = datafiles_all,
  endpoint = "All HFMD cases",
  scenario = "EV71 vaccine period",
  start_date = as.Date("2013-01-01"),
  end_date = as.Date("2020-01-01"),
  intervention_date = intervention_ev71,
  trend = "linear"
)

vax_ev71_primary <- run_counterfactual(
  df = datafiles_virus,
  endpoint = "Typed EV71 cases",
  scenario = "EV71 vaccine period",
  start_date = as.Date("2013-01-01"),
  end_date = as.Date("2020-01-01"),
  intervention_date = intervention_ev71,
  trend = "linear",
  virus_type = "EV71"
)

vax_cva16_primary <- run_counterfactual(
  df = datafiles_virus,
  endpoint = "Typed CV-A16 cases",
  scenario = "EV71 vaccine period",
  start_date = as.Date("2013-01-01"),
  end_date = as.Date("2020-01-01"),
  intervention_date = intervention_ev71,
  trend = "linear",
  virus_type = "CV-A16"
)

vax_others_primary <- run_counterfactual(
  df = datafiles_virus,
  endpoint = "Typed other-enterovirus cases",
  scenario = "EV71 vaccine period",
  start_date = as.Date("2013-01-01"),
  end_date = as.Date("2020-01-01"),
  intervention_date = intervention_ev71,
  trend = "linear",
  virus_type = "Others"
)

vax_all_sensitivity <- run_counterfactual(
  df = datafiles_all,
  endpoint = "All HFMD cases",
  scenario = "EV71 vaccine period",
  start_date = as.Date("2010-01-01"),
  end_date = as.Date("2020-01-01"),
  intervention_date = intervention_ev71,
  trend = "quadratic"
)

vax_ev71_sensitivity <- run_counterfactual(
  df = datafiles_virus,
  endpoint = "Typed EV71 cases",
  scenario = "EV71 vaccine period",
  start_date = as.Date("2010-01-01"),
  end_date = as.Date("2020-01-01"),
  intervention_date = intervention_ev71,
  trend = "quadratic",
  virus_type = "EV71"
)

vax_cva16_sensitivity <- run_counterfactual(
  df = datafiles_virus,
  endpoint = "Typed CV-A16 cases",
  scenario = "EV71 vaccine period",
  start_date = as.Date("2010-01-01"),
  end_date = as.Date("2020-01-01"),
  intervention_date = intervention_ev71,
  trend = "quadratic",
  virus_type = "CV-A16"
)

vax_others_sensitivity <- run_counterfactual(
  df = datafiles_virus,
  endpoint = "Typed other-enterovirus cases",
  scenario = "EV71 vaccine period",
  start_date = as.Date("2010-01-01"),
  end_date = as.Date("2020-01-01"),
  intervention_date = intervention_ev71,
  trend = "quadratic",
  virus_type = "Others"
)

covid_all_primary <- run_counterfactual(
  df = datafiles_all,
  endpoint = "All HFMD cases",
  scenario = "COVID-19 period",
  start_date = as.Date("2008-01-01"),
  end_date = as.Date("2024-01-01"),
  intervention_date = intervention_covid,
  trend = "linear"
)

covid_ev71_primary <- run_counterfactual(
  df = datafiles_virus,
  endpoint = "Typed EV71 cases",
  scenario = "COVID-19 period",
  start_date = as.Date("2008-01-01"),
  end_date = as.Date("2024-01-01"),
  intervention_date = intervention_covid,
  trend = "linear",
  virus_type = "EV71"
)

covid_cva16_primary <- run_counterfactual(
  df = datafiles_virus,
  endpoint = "Typed CV-A16 cases",
  scenario = "COVID-19 period",
  start_date = as.Date("2008-01-01"),
  end_date = as.Date("2024-01-01"),
  intervention_date = intervention_covid,
  trend = "linear",
  virus_type = "CV-A16"
)

covid_others_primary <- run_counterfactual(
  df = datafiles_virus,
  endpoint = "Typed other-enterovirus cases",
  scenario = "COVID-19 period",
  start_date = as.Date("2008-01-01"),
  end_date = as.Date("2024-01-01"),
  intervention_date = intervention_covid,
  trend = "linear",
  virus_type = "Others"
)

all_cases_limit <- range(
  c(
    vax_all_primary$n,
    vax_all_primary$counterfactual,
    covid_all_primary$n,
    covid_all_primary$counterfactual
  ),
  finite = TRUE
)

typed_ev71_limit <- range(
  c(
    vax_ev71_primary$n,
    vax_ev71_primary$counterfactual,
    covid_ev71_primary$n,
    covid_ev71_primary$counterfactual
  ),
  finite = TRUE
)

plots_vax_all <- plot_counterfactual(
  df = vax_all_primary,
  intervention_date = intervention_ev71,
  plot_title_main = "A: All reported HFMD cases, vaccine-period counterfactual",
  plot_title_diff = "B: Vaccine-period difference for all reported cases",
  x_limits = as.Date(c("2013-01-01", "2020-01-01")),
  y_limits = all_cases_limit
)

plots_vax_ev71 <- plot_counterfactual(
  df = vax_ev71_primary,
  intervention_date = intervention_ev71,
  plot_title_main = "C: Typed EV71 cases, vaccine-period counterfactual",
  plot_title_diff = "D: Vaccine-period difference for typed EV71 cases",
  x_limits = as.Date(c("2013-01-01", "2020-01-01")),
  y_limits = typed_ev71_limit
)

plots_covid_all <- plot_counterfactual(
  df = covid_all_primary,
  intervention_date = intervention_covid,
  plot_title_main = "E: All reported HFMD cases, COVID-19-period counterfactual",
  plot_title_diff = "F: COVID-19-period difference for all reported cases",
  x_limits = as.Date(c("2008-01-01", "2024-01-01")),
  y_limits = all_cases_limit
)

plots_covid_ev71 <- plot_counterfactual(
  df = covid_ev71_primary,
  intervention_date = intervention_covid,
  plot_title_main = "G: Typed EV71 cases, COVID-19-period counterfactual",
  plot_title_diff = "H: COVID-19-period difference for typed EV71 cases",
  x_limits = as.Date(c("2008-01-01", "2024-01-01")),
  y_limits = typed_ev71_limit
)

fig_main <- (
  plots_vax_all$main + plots_vax_all$diff +
  plots_vax_ev71$main + plots_vax_ev71$diff
) / (
  plots_covid_all$main + plots_covid_all$diff +
  plots_covid_ev71$main + plots_covid_ev71$diff
) +
  plot_layout(widths = c(2, 1, 2, 1), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "fig5.pdf",
  plot = fig_main,
  width = 18,
  height = 10,
  device = cairo_pdf,
  family = "Times New Roman"
)

ggsave(
  "fig5.png",
  plot = fig_main,
  width = 18,
  height = 10,
  dpi = 320
)

plots_vax_cva16 <- plot_counterfactual(
  df = vax_cva16_primary,
  intervention_date = intervention_ev71,
  plot_title_main = "A: Typed CV-A16 cases, vaccine-period counterfactual",
  plot_title_diff = "B: Vaccine-period difference for typed CV-A16 cases",
  x_limits = as.Date(c("2013-01-01", "2020-01-01"))
)

plots_vax_others <- plot_counterfactual(
  df = vax_others_primary,
  intervention_date = intervention_ev71,
  plot_title_main = "C: Typed other-enterovirus cases, vaccine-period counterfactual",
  plot_title_diff = "D: Vaccine-period difference for typed other-enterovirus cases",
  x_limits = as.Date(c("2013-01-01", "2020-01-01"))
)

plots_covid_cva16 <- plot_counterfactual(
  df = covid_cva16_primary,
  intervention_date = intervention_covid,
  plot_title_main = "E: Typed CV-A16 cases, COVID-19-period counterfactual",
  plot_title_diff = "F: COVID-19-period difference for typed CV-A16 cases",
  x_limits = as.Date(c("2008-01-01", "2024-01-01"))
)

plots_covid_others <- plot_counterfactual(
  df = covid_others_primary,
  intervention_date = intervention_covid,
  plot_title_main = "G: Typed other-enterovirus cases, COVID-19-period counterfactual",
  plot_title_diff = "H: COVID-19-period difference for typed other-enterovirus cases",
  x_limits = as.Date(c("2008-01-01", "2024-01-01"))
)

fig_secondary <- (
  plots_vax_cva16$main + plots_vax_cva16$diff +
  plots_vax_others$main + plots_vax_others$diff
) / (
  plots_covid_cva16$main + plots_covid_cva16$diff +
  plots_covid_others$main + plots_covid_others$diff
) +
  plot_layout(widths = c(2, 1, 2, 1), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "figS7_typed_secondary_counterfactual.pdf",
  plot = fig_secondary,
  width = 18,
  height = 10,
  device = cairo_pdf,
  family = "Times New Roman"
)

ggsave(
  "figS7_typed_secondary_counterfactual.png",
  plot = fig_secondary,
  width = 18,
  height = 10,
  dpi = 320
)

its_summary <- bind_rows(
  summarize_post_period(vax_all_primary, "EV71 vaccine period", "All reported HFMD cases", "primary_train_2013_2016", intervention_ev71),
  summarize_post_period(vax_ev71_primary, "EV71 vaccine period", "Typed EV71 cases", "primary_train_2013_2016", intervention_ev71),
  summarize_post_period(vax_cva16_primary, "EV71 vaccine period", "Typed CV-A16 cases", "primary_train_2013_2016", intervention_ev71),
  summarize_post_period(vax_others_primary, "EV71 vaccine period", "Typed other-enterovirus cases", "primary_train_2013_2016", intervention_ev71),
  summarize_post_period(covid_all_primary, "COVID-19 period", "All reported HFMD cases", "primary_train_2008_2019", intervention_covid),
  summarize_post_period(covid_ev71_primary, "COVID-19 period", "Typed EV71 cases", "primary_train_2008_2019", intervention_covid),
  summarize_post_period(covid_cva16_primary, "COVID-19 period", "Typed CV-A16 cases", "primary_train_2008_2019", intervention_covid),
  summarize_post_period(covid_others_primary, "COVID-19 period", "Typed other-enterovirus cases", "primary_train_2008_2019", intervention_covid)
)

training_window_sensitivity <- bind_rows(
  summarize_post_period(vax_all_primary, "EV71 vaccine period", "All reported HFMD cases", "primary_train_2013_2016", intervention_ev71),
  summarize_post_period(vax_all_sensitivity, "EV71 vaccine period", "All reported HFMD cases", "sensitivity_train_2010_2016", intervention_ev71),
  summarize_post_period(vax_ev71_primary, "EV71 vaccine period", "Typed EV71 cases", "primary_train_2013_2016", intervention_ev71),
  summarize_post_period(vax_ev71_sensitivity, "EV71 vaccine period", "Typed EV71 cases", "sensitivity_train_2010_2016", intervention_ev71),
  summarize_post_period(vax_cva16_primary, "EV71 vaccine period", "Typed CV-A16 cases", "primary_train_2013_2016", intervention_ev71),
  summarize_post_period(vax_cva16_sensitivity, "EV71 vaccine period", "Typed CV-A16 cases", "sensitivity_train_2010_2016", intervention_ev71),
  summarize_post_period(vax_others_primary, "EV71 vaccine period", "Typed other-enterovirus cases", "primary_train_2013_2016", intervention_ev71),
  summarize_post_period(vax_others_sensitivity, "EV71 vaccine period", "Typed other-enterovirus cases", "sensitivity_train_2010_2016", intervention_ev71)
)

write_csv(its_summary, "supplementary_tables/tableS4_primary_its_counterfactual.csv")
write_csv(training_window_sensitivity, "supplementary_tables/tableS8_vaccine_training_window_sensitivity.csv")
