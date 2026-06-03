rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(ggplot2)
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
load(file.path(data_dir, "clean2.RData"))

dir.create("supplementary_tables", showWarnings = FALSE, recursive = TRUE)

age_shift_df <- datafiles_2 %>%
  mutate(
    year = year(Date),
    age_num = suppressWarnings(as.numeric(if_else(age == "＜1", "0", as.character(age)))),
    age_group = case_when(
      age_num <= 1 ~ "<1",
      age_num %in% 2:3 ~ "2-3",
      age_num >= 4 ~ "4+",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(year >= 2008, year <= 2023)

annual_age_composition <- age_shift_df %>%
  filter(!is.na(age_group)) %>%
  count(year, age_group, name = "cases") %>%
  group_by(year) %>%
  mutate(prop_pct = cases / sum(cases) * 100) %>%
  ungroup() %>%
  mutate(age_group = factor(age_group, levels = c("<1", "2-3", "4+")))

period_age_comparison <- bind_rows(
  age_shift_df %>%
    filter(year %in% 2017:2019, !is.na(age_group)) %>%
    count(age_group, name = "cases") %>%
    mutate(period = "2017-2019"),
  age_shift_df %>%
    filter(year %in% 2021:2023, !is.na(age_group)) %>%
    count(age_group, name = "cases") %>%
    mutate(period = "2021-2023")
) %>%
  group_by(period) %>%
  mutate(prop_pct = cases / sum(cases) * 100) %>%
  ungroup() %>%
  mutate(age_group = factor(age_group, levels = c("<1", "2-3", "4+")))

chisq_input <- age_shift_df %>%
  filter(year %in% c(2017:2019, 2021:2023), !is.na(age_group)) %>%
  mutate(period = if_else(year <= 2019, "2017-2019", "2021-2023"))

chisq_fit <- chisq.test(table(chisq_input$period, chisq_input$age_group))

older_model_df <- chisq_input %>%
  transmute(
    post_period = if_else(period == "2021-2023", 1L, 0L),
    older = if_else(age_group == "4+", 1L, 0L)
  )

older_fit <- glm(older ~ post_period, family = binomial(), data = older_model_df)
coef_mat <- coef(summary(older_fit))
beta <- coef_mat["post_period", "Estimate"]
se <- coef_mat["post_period", "Std. Error"]

age_shift_statistics <- tibble(
  comparison = "2021-2023 versus 2017-2019",
  chi_square = unname(chisq_fit$statistic),
  chi_square_df = unname(chisq_fit$parameter),
  chi_square_p = chisq_fit$p.value,
  older_group_or = exp(beta),
  older_group_lcl = exp(beta - 1.96 * se),
  older_group_ucl = exp(beta + 1.96 * se),
  older_group_p = coef_mat["post_period", "Pr(>|z|)"]
)

write_csv(annual_age_composition, "supplementary_tables/supplementary_data_D9_annual_age_composition.csv")
write_csv(period_age_comparison, "supplementary_tables/supplementary_data_D9_period_summary.csv")
write_csv(age_shift_statistics, "supplementary_tables/supplementary_data_D9_statistics.csv")

fig_age_shift <- ggplot(annual_age_composition, aes(x = year, y = prop_pct / 100, fill = age_group)) +
  geom_col(color = "white", linewidth = 0.2) +
  geom_vline(xintercept = c(2016.5, 2019.5), color = "black", linewidth = 0.35, linetype = "22") +
  scale_fill_manual(
    values = c("<1" = "#F2EBBBFF", "2-3" = "#3CC8C0FF", "4+" = "#04578CFF"),
    drop = FALSE
  ) +
  scale_x_continuous(breaks = seq(2008, 2023, by = 1), expand = expansion(add = c(0.3, 0.3))) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
  labs(
    title = "Annual age composition of reported HFMD cases with known age",
    x = NULL,
    y = "Share of reported cases",
    fill = "Age group (years)"
  ) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    panel.grid = element_blank(),
    text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    plot.title.position = "plot"
  )

ggsave(
  "figS6_age_shift.pdf",
  plot = fig_age_shift,
  width = 11,
  height = 5.8,
  device = cairo_pdf,
  family = "Times New Roman"
)

ggsave(
  "figS6_age_shift.png",
  plot = fig_age_shift,
  width = 11,
  height = 5.8,
  dpi = 320
)
