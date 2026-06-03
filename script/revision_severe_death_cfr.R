library(dplyr)
library(lubridate)
library(readr)

if (dir.exists("../outcome")) {
  setwd("../outcome")
} else if (dir.exists("outcome")) {
  setwd("outcome")
} else {
  stop("Cannot find outcome/ directory from current working directory.")
}

load("../data/datafiles_all.RData")

yes_value <- intToUtf8(0x662F)

severe_death_yearly <- datafiles_all %>%
  mutate(
    onset_date = as.Date(.data[["发病日期"]]),
    year = year(onset_date),
    death_date = trimws(as.character(.data[["死亡日期"]])),
    severe_value = trimws(as.character(.data[["重症患者"]])),
    death = !is.na(death_date) & !(death_date %in% c("", ".", "NA")),
    severe = severe_value == yes_value
  ) %>%
  filter(year >= 2008, year <= 2023) %>%
  group_by(year) %>%
  summarise(
    reported_cases = n(),
    severe_cases = sum(severe, na.rm = TRUE),
    deaths = sum(death, na.rm = TRUE),
    cfr_pct = deaths / reported_cases * 100,
    deaths_per_100000_reported_cases = deaths / reported_cases * 100000,
    .groups = "drop"
  )

severe_death_overall <- severe_death_yearly %>%
  summarise(
    period = "2008-2023",
    reported_cases = sum(reported_cases),
    severe_cases = sum(severe_cases),
    deaths = sum(deaths),
    cfr_pct = deaths / reported_cases * 100,
    deaths_per_100000_reported_cases = deaths / reported_cases * 100000
  )

dir.create("supplementary_tables", showWarnings = FALSE, recursive = TRUE)
write_csv(severe_death_yearly, "supplementary_tables/severe_death_cfr_yearly.csv")
write_csv(severe_death_overall, "supplementary_tables/severe_death_cfr_overall.csv")

cat("Wrote supplementary_tables/severe_death_cfr_yearly.csv\n")
cat("Wrote supplementary_tables/severe_death_cfr_overall.csv\n")
