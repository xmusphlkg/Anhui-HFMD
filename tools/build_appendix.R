rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
repo_root <- normalizePath(file.path(script_dir, ".."))

template_path <- file.path(repo_root, "appendix_temp.md")
output_path <- file.path(repo_root, "appendix.md")
tables_dir <- file.path(repo_root, "outcome", "supplementary_tables")

normalize_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  replacements <- c(
    "（" = "(",
    "）" = ")",
    "％" = "%",
    "＜" = "<",
    "＞" = ">",
    "，" = ",",
    "\u00A0" = " "
  )
  for (pattern in names(replacements)) {
    x <- gsub(pattern, replacements[[pattern]], x, fixed = TRUE)
  }
  x <- trimws(gsub("[[:space:]]+", " ", x))
  x
}

escape_md <- function(x) {
  x <- normalize_text(x)
  x <- gsub("\\|", "\\\\|", x, perl = TRUE)
  x <- gsub("\r?\n", "<br>", x, perl = TRUE)
  x
}

render_md_table <- function(df) {
  df[] <- lapply(df, escape_md)
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  paste(c(header, separator, rows), collapse = "\n")
}

fmt_p <- function(x) {
  ifelse(
    is.na(x),
    "",
    ifelse(
      x < 0.001,
      "<.001",
      ifelse(
        x > 0.99,
        ">.99",
        sub("^0", "", ifelse(x < 0.01, sprintf("%.3f", x), sprintf("%.2f", x)))
      )
    )
  )
}

fmt_num <- function(x, digits = 0) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits, drop0trailing = FALSE))
}

fmt_int <- function(x) {
  fmt_num(round(x), digits = 0)
}

fmt_prop <- function(x) {
  sprintf("%.6f", x)
}

build_table_s1 <- function() {
  raw <- readxl::read_excel(
    file.path(repo_root, "outcome", "tableS1.xlsx"),
    col_names = FALSE,
    .name_repair = "minimal"
  )
  years <- normalize_text(unlist(raw[2, 2:15]))
  body <- as.data.frame(raw[3:22, 1:15], stringsAsFactors = FALSE)
  names(body) <- c("Characteristic", years)
  body[] <- lapply(body, normalize_text)
  body$Characteristic <- dplyr::recode(
    body$Characteristic,
    "Cox A16" = "CV-A16",
    "Other" = "Others",
    .default = body$Characteristic
  )
  body <- body[rowSums(body != "") > 0, , drop = FALSE]

  coverage <- readr::read_csv(file.path(tables_dir, "annual_typing_coverage.csv"), show_col_types = FALSE) %>%
    filter(year >= 2010, year <= 2023)

  fmt_count_pct <- function(count, pct) {
    sprintf("%d(%.2f%%)", round(count), pct)
  }

  virus_header_idx <- which(body$Characteristic == "Virus")
  virus_rows <- data.frame(
    Characteristic = c("Total", "EV71", "CV-A16", "Others"),
    stringsAsFactors = FALSE
  )
  for (yr in coverage$year) {
    yr_chr <- as.character(yr)
    row_dat <- coverage %>% filter(year == yr)
    virus_rows[[yr_chr]] <- c(
      sprintf("%d(100%%)", round(row_dat$non_missing_typed_cases)),
      fmt_count_pct(row_dat$ev71_typed_cases, row_dat$ev71_share_pct),
      fmt_count_pct(row_dat$cva16_typed_cases, row_dat$cva16_share_pct),
      fmt_count_pct(row_dat$other_typed_cases, row_dat$other_share_pct)
    )
  }

  body[(virus_header_idx + 1):(virus_header_idx + 4), ] <- virus_rows
  body
}

build_table_s2 <- function() {
  readxl::read_excel(file.path(repo_root, "outcome", "tableS2.xlsx")) %>%
    transmute(
      Name = name,
      Radius = fmt_num(RADIUS, 2),
      `Start date` = START_DATE,
      `End date` = END_DATE,
      `P value` = fmt_p(P_VALUE),
      Observed = fmt_int(OBSERVED),
      Expected = fmt_num(EXPECTED, 2),
      `Relative risk` = fmt_num(REL_RISK, 2),
      Year = fmt_int(year),
      LLR = fmt_num(LLR, 2)
    )
}

build_table_s3 <- function() {
  readxl::read_excel(file.path(repo_root, "outcome", "tableS2 moran.xlsx")) %>%
    transmute(
      Year = year,
      `Moran's I` = fmt_num(moran_i, 2),
      `Z score` = fmt_num(z_score, 2),
      `P value` = fmt_p(p_value)
    )
}

endpoint_display <- c(
  "All reported HFMD cases" = "All reported HFMD",
  "Typed EV71 cases" = "Typed EV71",
  "Typed CV-A16 cases" = "Typed CV-A16",
  "Typed other-enterovirus cases" = "Typed other enteroviruses"
)

build_table_s4 <- function() {
  readr::read_csv(file.path(tables_dir, "tableS4_primary_its_counterfactual.csv"), show_col_types = FALSE) %>%
    mutate(
      scenario_label = if_else(scenario == "EV71 vaccine period", "Vaccine-era", "COVID-19-period"),
      period = if_else(scenario == "EV71 vaccine period", "2017-2019", "2020-2023"),
      endpoint_label = endpoint_display[endpoint]
    ) %>%
    transmute(
      Model = paste0(endpoint_label, " (", scenario_label, ")"),
      Period = period,
      Year = fmt_int(year),
      `Actual cases` = fmt_int(actual_cases),
      `Predicted cases` = fmt_int(predicted_cases),
      `Lower 95% PI` = fmt_int(lower_95),
      `Upper 95% PI` = fmt_int(upper_95),
      Difference = fmt_int(difference),
      `Percent change (%)` = fmt_num(percent_change, 2)
    )
}

build_table_s5 <- function() {
  readr::read_csv(file.path(tables_dir, "annual_typing_coverage.csv"), show_col_types = FALSE) %>%
    transmute(
      Year = fmt_int(year),
      `Reported cases` = fmt_int(reported_cases),
      `Records with serotype field` = fmt_int(records_with_serotype_field),
      `Typed cases` = fmt_int(non_missing_typed_cases),
      `EV71 typed cases` = fmt_int(ev71_typed_cases),
      `CV-A16 typed cases` = fmt_int(cva16_typed_cases),
      `Other typed cases` = fmt_int(other_typed_cases),
      `Typed fraction (%)` = fmt_num(typed_fraction_pct, 2),
      `EV71 share (%)` = fmt_num(ev71_share_pct, 2),
      `CV-A16 share (%)` = fmt_num(cva16_share_pct, 2),
      `Other share (%)` = fmt_num(other_share_pct, 2)
    )
}

build_table_s6 <- function() {
  readr::read_csv(file.path(tables_dir, "severe_proportion_its_summary.csv"), show_col_types = FALSE) %>%
    transmute(
      Scenario = scenario,
      Period = period,
      Year = fmt_int(year),
      `Actual severe cases` = fmt_int(actual_severe),
      `Actual total cases` = fmt_int(actual_total),
      `Observed severe proportion` = fmt_prop(actual_prop),
      `Expected severe proportion` = fmt_prop(expected_prop),
      `Expected proportion LCL` = fmt_prop(expected_prop_lcl),
      `Expected proportion UCL` = fmt_prop(expected_prop_ucl)
    )
}

build_table_s7 <- function() {
  yearly <- readr::read_csv(file.path(tables_dir, "severe_death_cfr_yearly.csv"), show_col_types = FALSE)
  overall <- readr::read_csv(file.path(tables_dir, "severe_death_cfr_overall.csv"), show_col_types = FALSE)

  bind_rows(
    yearly %>%
      transmute(
        Year = as.character(year),
        `Reported cases` = reported_cases,
        `Severe cases` = severe_cases,
        Deaths = deaths,
        `Deaths per 100,000 reported cases` = deaths_per_100000_reported_cases
      ),
    overall %>%
      transmute(
        Year = period,
        `Reported cases` = reported_cases,
        `Severe cases` = severe_cases,
        Deaths = deaths,
        `Deaths per 100,000 reported cases` = deaths_per_100000_reported_cases
      )
  ) %>%
    transmute(
      Year = Year,
      `Reported cases` = fmt_int(`Reported cases`),
      `Severe cases` = fmt_int(`Severe cases`),
      Deaths = fmt_int(Deaths),
      `Deaths per 100,000 reported cases` = fmt_num(`Deaths per 100,000 reported cases`, 2)
    )
}

build_table_s8 <- function() {
  readr::read_csv(file.path(tables_dir, "typed_composition_multinom.csv"), show_col_types = FALSE) %>%
    mutate(
      term = dplyr::recode(
        term,
        "Intercept" = "Intercept",
        "year" = "Year",
        "gendermale" = "Male (vs female)",
        "age_group2-3" = "Age 2-3 years (vs <1 year)",
        "age_group4+" = "Age 4+ years (vs <1 year)",
        "severeY" = "Severe (vs mild)",
        .default = term
      )
    ) %>%
    transmute(
      Outcome = outcome,
      Term = term,
      Estimate = fmt_num(estimate, 2),
      `Standard error` = fmt_num(std_error, 2),
      RR = fmt_num(rr, 2),
      `RR LCL` = fmt_num(rr_lcl, 2),
      `RR UCL` = fmt_num(rr_ucl, 2)
    )
}

build_table_s9 <- function() {
  model_display <- c(
    "primary_train_2013_2016" = "Primary 2013-2016 linear",
    "sensitivity_train_2010_2016" = "Sensitivity 2010-2016 quadratic"
  )

  readr::read_csv(file.path(tables_dir, "tableS8_vaccine_training_window_sensitivity.csv"), show_col_types = FALSE) %>%
    filter(scenario == "EV71 vaccine period") %>%
    mutate(
      endpoint_label = endpoint_display[endpoint],
      model_label = model_display[model_spec]
    ) %>%
    transmute(
      `Endpoint and model` = paste(endpoint_label, model_label, sep = " - "),
      Year = fmt_int(year),
      `Actual cases` = fmt_int(actual_cases),
      `Predicted cases` = fmt_int(predicted_cases),
      `Lower 95% PI` = fmt_int(lower_95),
      `Upper 95% PI` = fmt_int(upper_95),
      Difference = fmt_int(difference),
      `Percent change (%)` = fmt_num(percent_change, 2)
    )
}

build_table_s10 <- function() {
  readr::read_csv(file.path(tables_dir, "covid_pre_vaccine_baseline_sensitivity.csv"), show_col_types = FALSE) %>%
    mutate(endpoint_label = endpoint_display[endpoint]) %>%
    transmute(
      `Endpoint and model` = paste(endpoint_label, "Pre-vaccine 2008-2015 linear", sep = " - "),
      Year = fmt_int(year),
      `Actual cases` = fmt_int(actual_cases),
      `Predicted cases` = fmt_int(predicted_cases),
      `Lower 95% PI` = fmt_int(lower_95),
      `Upper 95% PI` = fmt_int(upper_95),
      Difference = fmt_int(difference),
      `Percent change (%)` = fmt_num(percent_change, 2)
    )
}

appendix_template <- readChar(template_path, file.info(template_path)$size, useBytes = TRUE)

table_map <- list(
  "{{table_content_s1}}" = render_md_table(build_table_s1()),
  "{{table_content_s2}}" = render_md_table(build_table_s2()),
  "{{table_content_s3}}" = render_md_table(build_table_s3()),
  "{{table_content_s4}}" = render_md_table(build_table_s4()),
  "{{table_content_s5}}" = render_md_table(build_table_s5()),
  "{{table_content_s6}}" = render_md_table(build_table_s6()),
  "{{table_content_s7}}" = render_md_table(build_table_s7()),
  "{{table_content_s8}}" = render_md_table(build_table_s8()),
  "{{table_content_s9}}" = render_md_table(build_table_s9()),
  "{{table_content_s10}}" = render_md_table(build_table_s10())
)

for (token in names(table_map)) {
  appendix_template <- gsub(token, table_map[[token]], appendix_template, fixed = TRUE)
}

writeLines(appendix_template, output_path, useBytes = TRUE)
cat("Wrote ", output_path, "\n", sep = "")
