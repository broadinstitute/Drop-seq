library(tidyverse)

input_path <- '/Users/kshakir/Downloads/output.tsv'
display_read_num <- 1016

df <-
  input_path %>%
  read_delim(
    delim = "\t",
    col_types = cols(
      READ_NUM = col_number(),
      POSITION = col_number(),
      BASES = col_character()
    )
  )

sub_df <-
  df %>%
  filter(READ_NUM == display_read_num) %>%
  pivot_longer(
    cols = c('AVG', 'BASELINE', 'CAG', 'CCG', 'CCA'),
    names_to = 'TRIPLET',
    values_to = 'SCORE'
  ) %>%
  filter(BASES == TRIPLET)

sub_df_max <- max(sub_df$SCORE)

sub_df_plot <-
  sub_df %>%
  ggplot(aes(x = POSITION, y = SCORE, shape = TRIPLET)) +
  geom_point(aes(color = TRIPLET)) +
  ylim(0, sub_df_max) +
  theme_light()

show(sub_df_plot)
