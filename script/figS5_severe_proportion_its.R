rm(list = ls())

library(dplyr)
library(ggplot2)
library(lubridate)
library(patchwork)
library(readr)

if (dir.exists('../outcome')) {
  setwd('../outcome')
} else if (dir.exists('outcome')) {
  setwd('outcome')
} else {
  stop('Cannot find outcome/ directory from current working directory.')
}

data_dir <- if (dir.exists('../data_public')) '../data_public' else '../data'
load(file.path(data_dir, 'datafiles_all.RData'))

scientific_10 <- function(x) {
  ifelse(x == 0, 0, parse(text = gsub("[+]", "", gsub("e", "%*%10^", scales::scientific_format()(x)))))
}

prepare_monthly_severe <- function(data, end_date, start_date_filter = as.Date('2008-01-01')) {
  if (!('发病日期' %in% names(data))) stop('Expected column 发病日期 in datafiles_all')
  if (!('重症患者' %in% names(data))) stop('Expected column 重症患者 in datafiles_all')

  df <- data %>%
    mutate(
      date = as.Date(.data[['发病日期']]),
      severe_flag = .data[['重症患者']] == '是'
    ) %>%
    filter(date < end_date, date >= start_date_filter)

  monthly <- df %>%
    mutate(date = as.Date(format(date, '%Y-%m-01'))) %>%
    group_by(date) %>%
    summarise(
      total = n(),
      severe = sum(severe_flag, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(date) %>%
    mutate(
      month = factor(month(date), levels = 1:12),
      time_index = row_number(),
      prop = severe / total
    )

  monthly
}

build_and_predict_binom_its <- function(monthly, intervention_start_date, trend = c('linear', 'quadratic')) {
  trend <- match.arg(trend)
  pre <- monthly %>% filter(date < intervention_start_date)
  if (nrow(pre) < 24) stop('Pre-intervention window is too short for a stable seasonal model.')

  model_formula <- if (trend == 'quadratic') {
    cbind(severe, total - severe) ~ poly(time_index, 2, raw = TRUE) + month
  } else {
    cbind(severe, total - severe) ~ time_index + month
  }

  model <- glm(model_formula, data = pre, family = binomial())

  pred <- predict(model, newdata = monthly, type = 'link', se.fit = TRUE)
  fit <- as.numeric(pred$fit)
  se <- as.numeric(pred$se.fit)

  invlogit <- function(x) 1 / (1 + exp(-x))
  cf_mean <- invlogit(fit)
  cf_lcl <- invlogit(fit - 1.96 * se)
  cf_ucl <- invlogit(fit + 1.96 * se)

  is_post <- monthly$date >= intervention_start_date

  combined <- data.frame(
    Date = monthly$date,
    Total = monthly$total,
    Severe = monthly$severe,
    Actual = monthly$prop,
    Counterfactual = ifelse(is_post, cf_mean, NA_real_),
    Lower_95 = ifelse(is_post, cf_lcl, NA_real_),
    Upper_95 = ifelse(is_post, cf_ucl, NA_real_)
  ) %>%
    mutate(
      Diff = Actual - Counterfactual,
      Change = ifelse(Diff >= 0, 'Increase', 'Decrease')
    )

  list(model = model, combined_data = combined)
}

plot_severe_its <- function(combined_data, intervention_date, title_A, title_B, x_axis_limits) {
  min_date_plot <- x_axis_limits[1]
  max_date_plot <- x_axis_limits[2]

  plot_A <- ggplot(combined_data, aes(x = Date)) +
    geom_line(aes(y = Actual, color = 'Observed')) +
    geom_line(aes(y = Counterfactual, color = 'Predicted')) +
    geom_vline(xintercept = as.numeric(intervention_date), color = 'black') +
    scale_color_manual(values = c('Observed' = "#3CC8C0FF", 'Predicted' = "#F18B00FF"), drop = FALSE) +
    labs(title = title_A, y = 'Severe proportion', x = NULL, color = NULL) +
    theme_bw(base_family = 'Times New Roman') +
    theme(
      text = element_text(color = 'black', size = 13.5, family = 'Times New Roman'),
      panel.grid = element_blank(),
      axis.text = element_text(size = 12, color = 'black', family = 'Times New Roman'),
      axis.title.y = element_text(size = 14, face = 'bold', family = 'Times New Roman'),
      plot.title = element_text(size = 16, face = 'bold', family = 'Times New Roman'),
      plot.title.position = 'plot'
    ) +
    scale_x_date(
      breaks = seq(min_date_plot, max_date_plot, by = '2 years'),
      limits = c(min_date_plot, max_date_plot),
      expand = c(0, 0),
      date_labels = '%Y'
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0, 0))

  plot_B <- ggplot(combined_data, aes(x = Date, y = Diff, fill = Change)) +
    geom_col(alpha = 0.3, show.legend = FALSE) +
    geom_hline(yintercept = 0, color = 'black') +
    scale_fill_manual(values = c('Increase' = "#C0392B", 'Decrease' = "#27AE60"), drop = FALSE) +
    labs(title = title_B, y = 'Absolute difference (proportion)', x = NULL) +
    theme_bw(base_family = 'Times New Roman') +
    theme(
      text = element_text(color = 'black', size = 13.5, family = 'Times New Roman'),
      panel.grid = element_blank(),
      axis.text = element_text(size = 12, color = 'black', family = 'Times New Roman'),
      axis.title.y = element_text(size = 14, face = 'bold', family = 'Times New Roman'),
      plot.title = element_text(size = 16, face = 'bold', family = 'Times New Roman'),
      plot.title.position = 'plot'
    ) +
    scale_x_date(
      breaks = seq(min_date_plot, max_date_plot, by = '2 years'),
      limits = c(min_date_plot, max_date_plot),
      expand = c(0, 0),
      date_labels = '%Y'
    )

  list(plot_A = plot_A, plot_B = plot_B)
}

# Vaccine impact (evaluate up to 2020-01)
intervention_ev71 <- as.Date('2016-06-01')
monthly_vaccine <- prepare_monthly_severe(datafiles_all, end_date = as.Date('2020-01-01'), start_date_filter = as.Date('2010-01-01'))
fit_vaccine <- build_and_predict_binom_its(monthly_vaccine, intervention_ev71, trend = 'quadratic')
plots_vaccine <- plot_severe_its(
  fit_vaccine$combined_data,
  intervention_ev71,
  title_A = 'A: EV71 vaccine',
  title_B = 'B: Difference',
  x_axis_limits = c(as.Date('2008-01-01'), as.Date('2020-01-01'))
)

# COVID impact (evaluate up to 2024-01)
intervention_covid <- as.Date('2020-01-01')
monthly_covid <- prepare_monthly_severe(datafiles_all, end_date = as.Date('2024-01-01'))
fit_covid <- build_and_predict_binom_its(monthly_covid, intervention_covid, trend = 'linear')
plots_covid <- plot_severe_its(
  fit_covid$combined_data,
  intervention_covid,
  title_A = 'C: COVID-19',
  title_B = 'D: Difference',
  x_axis_limits = c(as.Date('2008-01-01'), as.Date('2024-01-01'))
)

fig <- (plots_vaccine$plot_A + plots_vaccine$plot_B) /
  (plots_covid$plot_A + plots_covid$plot_B) +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom')

# Save

ggsave('figS5_severe_proportion_its.pdf', fig, width = 14, height = 10, device = cairo_pdf, family = 'Times New Roman')
ggsave('figS5_severe_proportion_its.png', fig, width = 14, height = 10)

# Summary table (yearly)
summary_yearly <- bind_rows(
  fit_vaccine$combined_data %>%
    filter(year(Date) >= 2017, year(Date) <= 2019) %>%
    mutate(year = year(Date), scenario = 'EV71 vaccine', period = '2017-2019') %>%
    group_by(scenario, period, year) %>%
    summarise(
      actual_severe = sum(Severe, na.rm = TRUE),
      actual_total = sum(Total, na.rm = TRUE),
      actual_prop = actual_severe / actual_total,
      expected_prop = mean(Counterfactual, na.rm = TRUE),
      expected_prop_lcl = mean(Lower_95, na.rm = TRUE),
      expected_prop_ucl = mean(Upper_95, na.rm = TRUE),
      .groups = 'drop'
    ),
  fit_covid$combined_data %>%
    filter(year(Date) >= 2020, year(Date) <= 2023) %>%
    mutate(year = year(Date), scenario = 'COVID-19', period = '2020-2023') %>%
    group_by(scenario, period, year) %>%
    summarise(
      actual_severe = sum(Severe, na.rm = TRUE),
      actual_total = sum(Total, na.rm = TRUE),
      actual_prop = actual_severe / actual_total,
      expected_prop = mean(Counterfactual, na.rm = TRUE),
      expected_prop_lcl = mean(Lower_95, na.rm = TRUE),
      expected_prop_ucl = mean(Upper_95, na.rm = TRUE),
      .groups = 'drop'
    )
)

dir.create('supplementary_tables', showWarnings = FALSE, recursive = TRUE)
write_csv(summary_yearly, 'supplementary_tables/severe_proportion_its_summary.csv')
