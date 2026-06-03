#####################################
## @Description: 
## @version: 
## @Author: Li Kangguo
## @Date: 2026-03-13 20:30:15
## @LastEditors: Li Kangguo
## @LastEditTime: 2026-03-13 20:30:22
#####################################
rm(list = ls())

library(dplyr)
library(lubridate)
library(readr)
library(nnet)

if (dir.exists('../outcome')) {
  setwd('../outcome')
} else if (dir.exists('outcome')) {
  setwd('outcome')
} else {
  stop('Cannot find outcome/ directory from current working directory.')
}

data_dir <- if (dir.exists('../data_public')) '../data_public' else '../data'
load(file.path(data_dir, 'datafile_virus.RData'))

# Prepare typed dataset

typed <- datafiles_virus %>%
  mutate(
    year = year(date),
    year_c = year - 2008,
    age_num = case_when(age == '＜1' ~ 0, TRUE ~ suppressWarnings(as.numeric(age))),
    age_group = case_when(
      is.na(age_num) ~ NA_character_,
      age_num <= 1 ~ '<1',
      age_num %in% 2:3 ~ '2-3',
      age_num >= 4 ~ '4+'
    ),
    gender = factor(gender),
    severe = factor(severe, levels = c('N', 'Y')),
    virus = case_when(
      virus %in% c('CV-A16', 'Cox A16') ~ 'CV-A16',
      virus == 'EV71' ~ 'EV71',
      is.na(virus) ~ NA_character_,
      TRUE ~ 'Others'
    ),
    virus = factor(virus, levels = c('Others', 'EV71', 'CV-A16'))
  ) %>%
  filter(year >= 2008, year <= 2023, !is.na(age_group), !is.na(gender), !is.na(severe), !is.na(virus))

# Multinomial model: baseline = Others
# Coefficients are log-relative-risk ratios vs baseline
model <- nnet::multinom(virus ~ year_c + gender + age_group + severe, data = typed, trace = FALSE)

s <- summary(model)
coef_mat <- s$coefficients
se_mat <- s$standard.errors

coef_df <- as.data.frame(as.table(coef_mat)) %>%
  rename(outcome = Var1, term = Var2, estimate = Freq) %>%
  left_join(
    as.data.frame(as.table(se_mat)) %>%
      rename(outcome = Var1, term = Var2, std_error = Freq),
    by = c('outcome', 'term')
  ) %>%
  mutate(
    rr = exp(estimate),
    rr_lcl = exp(estimate - 1.96 * std_error),
    rr_ucl = exp(estimate + 1.96 * std_error)
  ) %>%
  arrange(outcome, term)

# Clean term names
coef_df <- coef_df %>%
  mutate(
    term = case_when(
      term == '(Intercept)' ~ 'Intercept',
      term == 'year_c' ~ 'year',
      TRUE ~ term
    )
  )

# Save

dir.create('supplementary_tables', showWarnings = FALSE, recursive = TRUE)
write_csv(coef_df, 'supplementary_tables/typed_composition_multinom.csv')

# Also save basic yearly composition for context
composition_yearly <- typed %>%
  count(year, virus) %>%
  group_by(year) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

write_csv(composition_yearly, 'supplementary_tables/typed_composition_yearly.csv')
