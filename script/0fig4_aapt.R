#joinpoint 结果文件导入
rm(list=ls())

library(tidyverse)
library(readxl)
library(grid)
library(patchwork)
library(nih.joinpoint)
library(paletteer)

if (dir.exists('../outcome')) {
     setwd('../outcome')
} else if (dir.exists('outcome')) {
     setwd('outcome')
} else {
     stop('Cannot find outcome/ directory from current working directory.')
}
data_dir <- if (dir.exists('../data_public')) '../data_public' else '../data'
load(file.path(data_dir, 'clean2.RData'))

# population denominators (city_pop.xlsx stores population in units of /10,000)
pop_city <- read_excel(file.path(data_dir, 'city_pop.xlsx'))
pop_city_long <- pop_city |>
     pivot_longer(cols = -name, names_to = 'year', values_to = 'population') |>
     mutate(year = as.integer(year)) |>
     filter(year >= 2008, year <= 2023)

pop_province <- pop_city_long |>
     group_by(year) |>
     summarise(population = sum(population, na.rm = TRUE), .groups = 'drop')

# visual ------------------------------------------------------------------

scientific_10 <- function(x) {
     ifelse(x == 0, 0, parse(text = gsub("[+]", "", gsub("e", "%*%10^", scales::scientific_format()(x)))))
}

## visual aapc
plot_apc <- function(jp_model, data, show.legend = F) {
     df_jp_apc <- jp_model$apc |>
          # round result
          mutate(across(c(apc, apc_95_lcl, apc_95_ucl), ~round(., 2)),
                 p_value = as.numeric(p_value),
                 p_value_label = case_when(p_value < 0.001 ~ '***',
                                           p_value < 0.01 ~ '**',
                                           p_value < 0.05 ~ '*',
                                           TRUE ~ ''),
                 legend = paste0(segment_start, '~', segment_end, '\n',
                                 apc, '(', apc_95_lcl, '~', apc_95_ucl, ')', p_value_label))
     
     df_jp_model <- jp_model$data_export
     
     # get breaks of y axis
     breaks <- pretty(c(0, data$val, df_jp_model$model))
     
     # set colors
     colors <- tail(paletteer_d("futurevisions::atomic_clock"), n = nrow(df_jp_apc))
     colors <- colors[order(order(df_jp_apc$apc))]
     
     # browser()
     
     fig <- ggplot(data)+
          geom_vline(data = df_jp_apc,
                     mapping = aes(xintercept = segment_end),
                     alpha = 0.5,
                     color = 'grey50')+
          geom_rect(data = df_jp_apc,
                    aes(xmin = segment_start, xmax = segment_end,
                        ymin = min(breaks), ymax = max(breaks),
                        fill = legend),
                    alpha = 0.5)+
          geom_point(mapping = aes(x = year, y = val, color = 'Observed'),
                     show.legend = show.legend) +
          geom_line(data = df_jp_model,
                    mapping = aes(x = year, y = model, color = 'Fitted'),
                    show.legend = show.legend) +
          scale_x_continuous(breaks = seq(2008, 2023, 2),
                             expand = expansion(add = c(0.1, 0.1))) +
          scale_y_continuous(limit = range(breaks),
                             breaks = breaks,
                             labels = scientific_10,
                             expand = expansion(mult = c(0, 0))) +
          scale_color_manual(name = 'Type',
                             values = c('Observed' = "#3CC8C0FF", 'Fitted' = "#04578CFF")) +
          scale_fill_manual(name = 'APC (95% CI)',
                            values = colors)+
          theme_bw()+
          theme(plot.title.position = 'plot',
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.text = element_text(size = 12, color = "black",family = 'Times New Roman'),
                axis.title = element_text(size = 14, face = "bold",family = 'Times New Roman'),
                plot.title = element_text(size = 16, face = "bold", family = 'Times New Roman'),
                plot.background = element_rect(fill = "transparent", color = NA),
                legend.text = element_text(size = 10, family = 'Times New Roman'),
                legend.title = element_text(size = 12, family = 'Times New Roman'),
                legend.position="bottom",
                legend.box="vertical", 
                legend.box.just = 'top',
                legend.justification.bottom = 'center',
                legend.title.position = 'top',
                legend.margin = margin(t = 10, unit = "pt"),
                legend.key.spacing.y = unit(0.35, 'cm'))+
          guides(fill = guide_legend(ncol = 4, byrow = TRUE, nrow = 1, order = 1, override.aes = list(shape = NA, linetype = 0)),
                 color = guide_legend(order = 2, override.aes = list(linetype = c(1, 0), shape = c(NA, 16))))
     
     fig
}

# jp model ----------------------------------------------------------------

## joinpoint setting for rate
run_opt_number = run_options(model="ln",
                             model_selection_method = 'permutation test',
                             ci_method = 'parametric',
                             dependent_variable_type = 'count',
                             max_joinpoints = 3,
                             n_cores=parallel::detectCores())

export_opt = export_options(aapc_full_range  = TRUE,
                            export_aapc = TRUE)

## all age group --------------------------------------------------------

data_group_0 <- datafiles_2 |> 
     filter(Date < as.Date("2024-01-01")) |> 
     mutate(year = year(Date)) |> 
     group_by(year) |>
     summarise(cases = n(), .groups = 'drop') |>
     left_join(pop_province, by = 'year') |>
     mutate(val = (cases / population) * 10)

## build joinpoint model for number
model_number_incidence <- joinpoint(data_group_0,
                                    year,
                                    val,
                                    run_opt = run_opt_number,
                                    export_opt = export_opt)

# fig1
fig_a <- plot_apc(model_number_incidence, data_group_0)+
     guides(fill = guide_legend(ncol=4, byrow=TRUE))+
     labs(x=NULL,
          y='Incidence rate (per 100,000)',
          fill="APC (95% CI)",
          title='A: All cases')

## gender group --------------------------------------------------------

data_group_1 <-  datafiles_2 |> 
     filter(Date < as.Date("2024-01-01")) |> 
     mutate(year = year(Date)) |> 
     group_by(year, gender) |>
     summarise(cases = n(), .groups = 'drop') |>
     rename(group = gender) |>
     left_join(pop_province, by = 'year') |>
     mutate(val = (cases / population) * 10) |>
     ungroup()

## build joinpoint model for number
model_number_gender <- joinpoint(data_group_1,
                                 year,
                                 val,
                                 by = group,
                                 run_opt = run_opt_number,
                                 export_opt = export_opt)

model_number_male <- map(model_number_gender, ~ .x %>%
                              filter(group == 'male'))

model_number_female <- map(model_number_gender, ~ .x %>%
                                filter(group == 'female'))

# fig2
fig_b <- plot_apc(model_number_male, filter(data_group_1, group == 'male'))+
     labs(x=NULL,
          y=NULL,
          fill="APC (95% CI)",
          title='B: Male')


# fig3
fig_c <- plot_apc(model_number_female, filter(data_group_1, group == 'female'))+
     labs(x=NULL,
          y=NULL,
          fill="APC (95% CI)",
          title='C: Female')

## age group --------------------------------------------------------

data_group_2 <-  datafiles_2 |> 
     mutate(year = year(Date),
            age = case_when(age == '＜1' ~ 0,
                            TRUE ~ as.numeric(age)),
            age = case_when(age <= 1 ~ '<1',
                            age %in% 2:3 ~ '2-3',
                            age > 3 ~ '4+')) |>
     filter(Date < as.Date("2024-01-01")) |> 
     group_by(year, age) |>
     count() |> 
     rename(val = n, group = age) |> 
     ungroup()

## build joinpoint model for number
model_number_age <- joinpoint(data_group_2,
                              year,
                              val,
                              by = group,
                              run_opt = run_opt_number,
                              export_opt = export_opt)

model_number_age_0 <- map(model_number_age, ~ .x %>%
                               filter(group == '<1'))

model_number_age_1 <- map(model_number_age, ~ .x %>%
                               filter(group == '2-3'))

model_number_age_2 <- map(model_number_age, ~ .x %>%
                               filter(group == '4+'))

# fig4
fig_d <- plot_apc(model_number_age_0, filter(data_group_2, group == '<1'))+
     labs(x=NULL,
          y='Number of cases',
          fill="APC (95% CI)",
          title='D: Age ≤ 1')

# fig5
fig_e <- plot_apc(model_number_age_1, filter(data_group_2, group == '2-3'))+
     labs(x=NULL,
          y=NULL,
          fill="APC (95% CI)",
          title='E: Age 2-3')

# fig6
fig_f <- plot_apc(model_number_age_2, filter(data_group_2, group == '4+'))+
     labs(x=NULL,
          y=NULL,
          fill="APC (95% CI)",
          title='F: Age 4+')

## severe --------------------------------------------------------

data_group_3 <-  datafiles_2 |> 
     mutate(year = year(Date),
            virus = case_when(virus == 'EV71' ~ 'EV71',
                              virus == 'Cox A16' ~ 'CV-A16',
                              virus == "其他肠道病毒" ~ 'Other',
                              TRUE ~ virus)) |>
     filter(Date < as.Date("2024-01-01")) |>
     group_by(year, virus) |>
     count() |>
     rename(val = n, group = virus) |> 
     ungroup()

# build joinpoint model for number
model_number_virus <- joinpoint(data_group_3,
                                year,
                                val,
                                by = group,
                                run_opt = run_opt_number,
                                export_opt = export_opt)

model_number_ev71 <- map(model_number_virus, ~ .x %>%
                              filter(group == 'EV71'))

model_number_cv_a16 <- map(model_number_virus, ~ .x %>%
                                 filter(group == 'CV-A16'))

model_number_other <- map(model_number_virus, ~ .x %>%
                                 filter(group == 'Other'))

# fig7
fig_g <- plot_apc(model_number_ev71, filter(data_group_3, group == 'EV71'))+
     labs(x=NULL,
          y='Number of cases',
          fill="APC (95% CI)",
          title='G: EV71')+
     guides(fill = guide_legend(ncol = 3)) 

# fig8
fig_h <- plot_apc(model_number_cv_a16, filter(data_group_3, group == 'CV-A16'), show.legend = F)+
     labs(x=NULL,
          y=NULL,
          fill="APC (95% CI)",
          title='H: CV-A16')


# fig9
fig_i <- plot_apc(model_number_other, filter(data_group_3, group == 'Other'))+
     labs(x=NULL,
          y=NULL,
          fill="APC (95% CI)",
          title='I: Other')

# combine plots ----------------------------------------------------------
combined_plot <- 
     (fig_a + fig_b + fig_c) / 
     (fig_d + fig_e + fig_f) / 
     (fig_g + fig_h + fig_i)

ggsave(
  'fig4_apc.pdf', 
  plot = combined_plot,
  width = 16, 
  height = 10, 
  device = cairo_pdf,
  family = "Times New Roman"
)

ggsave(
     'fig4_apc.png', 
     plot = combined_plot,
     width = 16, 
     height = 10
)
