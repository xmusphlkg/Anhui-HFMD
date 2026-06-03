#####################################
## @Description: 
## @version: 
## @Author: Li Kangguo
## @Date: 2026-03-10 17:36:55
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-03-10 17:37:01
#####################################
rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(forecast)
  library(lubridate)
})

data_dir <- if (dir.exists('data_public')) 'data_public' else 'data'
load(file.path(data_dir, 'datafile_virus.RData'))
load(file.path(data_dir, 'datafiles_all.RData'))

cat("Loaded data objects.\n")
cat("Virus labels in datafile_virus.RData:", paste(sort(unique(as.character(datafiles_virus$virus))), collapse = " | "), "\n")

prepare_time_series_data <- function(data, virus_type = NULL, end_date, start_date_filter = as.Date("2008-01-01")) {
  date_col_name <- if ("发病日期" %in% names(data)) "发病日期" else "date"

  data <- data |>
    mutate(date = as.Date(.data[[date_col_name]]))

  data_filtered <- data |>
    filter(date < end_date, date >= start_date_filter)

  if (!is.null(virus_type)) {
    if (virus_type == "CV-A16") {
      data_filtered <- data_filtered |>
        filter(virus %in% c("CV-A16", "Cox A16"))
    } else if (virus_type == "Others") {
      data_filtered <- data_filtered |>
        filter(!is.na(virus), !(virus %in% c("EV71", "CV-A16", "Cox A16")))
    } else {
      data_filtered <- data_filtered |>
        filter(virus == virus_type)
    }
  }

  processed_data <- data_filtered |>
    transmute(date = as.Date(format(date, "%Y-%m-01"))) |>
    count(date, name = "n") |>
    arrange(date)

  ts_object <- ts(
    processed_data$n,
    start = c(year(min(processed_data$date)), month(min(processed_data$date))),
    frequency = 12
  )

  list(processed_data = processed_data, ts_object = ts_object)
}

format_arima_order <- function(model) {
  arma <- model$arma
  sprintf("(%d,%d,%d)(%d,%d,%d)[12]",
          arma[1], arma[6], arma[2], arma[3], arma[7], arma[4])
}

fit_counterfactual <- function(dataset, scenario, endpoint, virus_type = NULL, end_date, intervention_date) {
  cat("Fitting:", scenario, "-", endpoint, "\n")
  prepared <- prepare_time_series_data(dataset, virus_type = virus_type, end_date = end_date)
  processed_data <- prepared$processed_data

  if (nrow(processed_data) == 0) {
    stop(paste("No observations available for", scenario, endpoint))
  }

  pre_data <- processed_data |>
    filter(date < intervention_date)

  ts_pre <- ts(
    pre_data$n,
    start = c(year(min(pre_data$date)), month(min(pre_data$date))),
    frequency = 12
  )

  model <- auto.arima(
    ts_pre,
    seasonal = TRUE,
    stepwise = TRUE,
    approximation = TRUE,
    max.p = 3,
    max.q = 3,
    max.P = 2,
    max.Q = 2,
    max.order = 8,
    lambda = NULL
  )
  residuals_arima <- residuals(model)
  ljung_box <- Box.test(residuals_arima, lag = 12, type = "Ljung-Box")
  shapiro <- if (length(residuals_arima) >= 3 && length(residuals_arima) <= 5000) shapiro.test(residuals_arima) else NULL

  fitted_values <- fitted(model)
  train_actual <- as.numeric(ts_pre)
  train_mape <- mean(abs((train_actual - fitted_values) / pmax(train_actual, 1)), na.rm = TRUE) * 100

  forecast_horizon <- nrow(processed_data) - nrow(pre_data)
  forecast_obj <- forecast(model, h = forecast_horizon)

  combined_data <- tibble(
    date = processed_data$date,
    actual = processed_data$n,
    counterfactual = c(rep(NA_real_, nrow(pre_data)), pmax(as.numeric(forecast_obj$mean), 0)),
    lower_95 = c(rep(NA_real_, nrow(pre_data)), pmax(as.numeric(forecast_obj$lower[, 2]), 0)),
    upper_95 = c(rep(NA_real_, nrow(pre_data)), pmax(as.numeric(forecast_obj$upper[, 2]), 0))
  ) |>
    mutate(
      absolute_effect = actual - counterfactual,
      relative_effect_pct = ifelse(!is.na(counterfactual) & counterfactual > 0,
                                   (actual - counterfactual) / counterfactual * 100,
                                   NA_real_)
    )

  post_period <- combined_data |>
    filter(date >= intervention_date)

  summary_row <- tibble(
    scenario = scenario,
    endpoint = endpoint,
    model_order = format_arima_order(model),
    training_start = format(min(pre_data$date), "%Y-%m"),
    training_end = format(max(pre_data$date), "%Y-%m"),
    forecast_start = format(min(post_period$date), "%Y-%m"),
    forecast_end = format(max(post_period$date), "%Y-%m"),
    train_mape_pct = round(train_mape, 2),
    aic = round(AIC(model), 2),
    ljung_box_statistic = round(as.numeric(ljung_box$statistic), 3),
    ljung_box_p = signif(ljung_box$p.value, 3),
    shapiro_p = if (!is.null(shapiro)) signif(shapiro$p.value, 3) else NA_real_,
    max_relative_increase_pct = round(max(post_period$relative_effect_pct, na.rm = TRUE), 2),
    min_relative_change_pct = round(min(post_period$relative_effect_pct, na.rm = TRUE), 2),
    max_absolute_increase = round(max(post_period$absolute_effect, na.rm = TRUE), 2),
    max_absolute_decrease = round(min(post_period$absolute_effect, na.rm = TRUE), 2)
  )

  cat("Finished:", scenario, "-", endpoint, "\n")

  list(summary = summary_row, details = combined_data)
}

vaccine_summaries <- bind_rows(
  fit_counterfactual(datafiles_all, "EV71 vaccine rollout", "All HFMD cases", NULL, as.Date("2020-01-01"), as.Date("2016-06-01"))$summary,
  fit_counterfactual(datafiles_virus, "EV71 vaccine rollout", "EV71 typed cases", "EV71", as.Date("2020-01-01"), as.Date("2016-06-01"))$summary,
  fit_counterfactual(datafiles_virus, "EV71 vaccine rollout", "CV-A16 typed cases", "CV-A16", as.Date("2020-01-01"), as.Date("2016-06-01"))$summary,
  fit_counterfactual(datafiles_virus, "EV71 vaccine rollout", "Other typed cases", "Others", as.Date("2020-01-01"), as.Date("2016-06-01"))$summary
)

covid_summaries <- bind_rows(
  fit_counterfactual(datafiles_all, "COVID-19 period", "All HFMD cases", NULL, as.Date("2024-01-01"), as.Date("2020-01-01"))$summary,
  fit_counterfactual(datafiles_virus, "COVID-19 period", "EV71 typed cases", "EV71", as.Date("2024-01-01"), as.Date("2020-01-01"))$summary,
  fit_counterfactual(datafiles_virus, "COVID-19 period", "CV-A16 typed cases", "CV-A16", as.Date("2024-01-01"), as.Date("2020-01-01"))$summary,
  fit_counterfactual(datafiles_virus, "COVID-19 period", "Other typed cases", "Others", as.Date("2024-01-01"), as.Date("2020-01-01"))$summary
)

annual_total <- datafiles_all |>
  mutate(date = as.Date(`发病日期`), year = year(date)) |>
  filter(year >= 2008, year <= 2023) |>
  count(year, name = "reported_cases")

annual_typed <- datafiles_virus |>
  mutate(date = as.Date(date), year = year(date)) |>
  filter(year >= 2008, year <= 2023) |>
  group_by(year) |>
  summarise(
    records_with_serotype_field = n(),
    non_missing_typed_cases = n(),
    ev71_typed_cases = sum(virus == "EV71", na.rm = TRUE),
    cva16_typed_cases = sum(virus == "CV-A16" | virus == "Cox A16", na.rm = TRUE),
    other_typed_cases = sum(!is.na(virus) & !(virus %in% c("EV71", "CV-A16", "Cox A16")), na.rm = TRUE),
    .groups = "drop"
  )

typing_coverage <- annual_total |>
  left_join(annual_typed, by = "year") |>
  mutate(
    non_missing_typed_cases = coalesce(non_missing_typed_cases, 0L),
    records_with_serotype_field = coalesce(records_with_serotype_field, 0L),
    ev71_typed_cases = coalesce(ev71_typed_cases, 0L),
    cva16_typed_cases = coalesce(cva16_typed_cases, 0L),
    other_typed_cases = coalesce(other_typed_cases, 0L),
    typed_fraction_pct = round(non_missing_typed_cases / reported_cases * 100, 3),
    ev71_share_pct = round(ev71_typed_cases / pmax(non_missing_typed_cases, 1) * 100, 2),
    cva16_share_pct = round(cva16_typed_cases / pmax(non_missing_typed_cases, 1) * 100, 2),
    other_share_pct = round(other_typed_cases / pmax(non_missing_typed_cases, 1) * 100, 2)
  )

dir.create("outcome/supplementary_tables", showWarnings = FALSE, recursive = TRUE)

write.csv(bind_rows(vaccine_summaries, covid_summaries), "outcome/supplementary_tables/sarima_model_diagnostics.csv", row.names = FALSE)
write.csv(typing_coverage, "outcome/supplementary_tables/annual_typing_coverage.csv", row.names = FALSE)

cat("Wrote outcome/supplementary_tables/sarima_model_diagnostics.csv\n")
cat("Wrote outcome/supplementary_tables/annual_typing_coverage.csv\n")