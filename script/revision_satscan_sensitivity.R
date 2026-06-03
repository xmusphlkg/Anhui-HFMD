rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(lubridate)
  library(readr)
})

if (getRversion() >= '2.15.1') {
  utils::globalVariables(c('RADIUS', 'START_DATE', 'END_DATE', 'start_date'))
}

if (dir.exists('../outcome')) {
  setwd('../outcome')
} else if (dir.exists('outcome')) {
  setwd('outcome')
} else {
  stop('Cannot find outcome/ directory from current working directory.')
}

data_dir <- if (dir.exists('../data_public')) '../data_public' else '../data'

# Expected inputs:
# - data/spt.col.dbf: base SaTScan result (typically 15% max spatial size)
# - optional additional files for sensitivity runs, e.g.:
#   data/spt_10.col.dbf and data/spt_20.col.dbf
#   or data/spt10.col.dbf and data/spt20.col.dbf

candidate_files <- tibble::tribble(
  ~cluster_size_pct, ~path,
  10, file.path(data_dir, 'spt_10.col.dbf'),
  10, file.path(data_dir, 'spt10.col.dbf'),
  15, file.path(data_dir, 'spt.col.dbf'),
  15, file.path(data_dir, 'spt_15.col.dbf'),
  20, file.path(data_dir, 'spt_20.col.dbf'),
  20, file.path(data_dir, 'spt20.col.dbf')
)

existing <- candidate_files %>%
  filter(file.exists(path)) %>%
  group_by(cluster_size_pct) %>%
  slice_head(n = 1) %>%
  ungroup()

if (nrow(existing) == 0) {
  stop('No SaTScan DBF files found. Please place DBF outputs under data/.')
}

read_satscan <- function(file_path, cluster_size_pct) {
  dat <- sf::st_read(file_path, quiet = TRUE) %>%
    as.data.frame()

  if (!all(c('RADIUS', 'REL_RISK', 'P_VALUE', 'START_DATE', 'END_DATE', 'OBSERVED', 'EXPECTED', 'LLR') %in% names(dat))) {
    stop(paste('Missing required columns in', file_path))
  }

  dat %>%
    filter(RADIUS != 0) %>%
    mutate(
      cluster_size_pct = cluster_size_pct,
      start_date = suppressWarnings(ymd(START_DATE)),
      end_date = suppressWarnings(ymd(END_DATE)),
      year = year(start_date)
    )
}

all_clusters <- purrr::map2_dfr(existing$path, existing$cluster_size_pct, read_satscan)

summary_table <- all_clusters %>%
  group_by(cluster_size_pct) %>%
  summarise(
    n_clusters = n(),
    n_significant_p005 = sum(P_VALUE < 0.05, na.rm = TRUE),
    median_rr = median(REL_RISK, na.rm = TRUE),
    iqr_rr = IQR(REL_RISK, na.rm = TRUE),
    max_rr = max(REL_RISK, na.rm = TRUE),
    max_llr = max(LLR, na.rm = TRUE),
    median_radius = median(RADIUS, na.rm = TRUE),
    earliest_start = as.character(min(start_date, na.rm = TRUE)),
    latest_end = as.character(max(end_date, na.rm = TRUE)),
    monte_carlo_replications = NA_integer_,
    multiple_testing_correction = 'SaTScan default Monte Carlo inference',
    .groups = 'drop'
  ) %>%
  arrange(cluster_size_pct)

cluster_detail <- all_clusters %>%
  transmute(
    cluster_size_pct,
    start_date,
    end_date,
    p_value = P_VALUE,
    observed = OBSERVED,
    expected = EXPECTED,
    rel_risk = REL_RISK,
    llr = LLR,
    radius = RADIUS
  ) %>%
  arrange(cluster_size_pct, start_date)

dir.create('supplementary_tables', showWarnings = FALSE, recursive = TRUE)
write_csv(summary_table, 'supplementary_tables/satscan_cluster_size_sensitivity_summary.csv')
write_csv(cluster_detail, 'supplementary_tables/satscan_cluster_size_sensitivity_detail.csv')

cat('Wrote supplementary_tables/satscan_cluster_size_sensitivity_summary.csv\n')
cat('Wrote supplementary_tables/satscan_cluster_size_sensitivity_detail.csv\n')
cat('Note: Fill monte_carlo_replications manually from the SaTScan .prm or run log if needed.\n')
