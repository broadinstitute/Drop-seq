library(tidyverse)

input_path <- '/Users/kshakir/Downloads/output.tsv'

display_read_nums <- c(
  # Reads with longer Q2 tails
  9538,
  1245,
  991,
  3,
  6291,
  7167,
  2401,
  4308,
  3117,
  5632,
  # Other example "good" reads
  1010,
  1012,
  1016
)

df <-
  input_path %>%
    read_delim(
      delim = "\t",
      col_types = cols(
        READ_NUM = col_number(),
        TRIPLET_NUM = col_number(),
        BASES = col_character()
      )
    )

sub_df <-
  df %>%
    filter(READ_NUM %in% display_read_nums) %>%
    pivot_longer(
      cols = !c('READ_NUM', 'TRIPLET_NUM', 'BASES'),
      names_to = 'TRIPLET',
      values_to = 'SCORE'
    ) %>%
    filter(BASES == TRIPLET)

# Convert accuracy log likelihood (0 to negative inifinity)
# to error likelihood (0.0 to 1.0)
# then to phred (0 to positive infitity)
#sub_df$SCORE <- log10(1 - (10^sub_df$SCORE)) * -10.0

# If needed, the reverse
#sub_df$SCORE <- log10(1 - (10^(sub_df$SCORE / -10.0)))

sub_df_min <- min(0, sub_df$SCORE)
sub_df_max <- max(0, sub_df$SCORE)

sub_df_plot <-
  sub_df %>%
    ggplot(aes(x = TRIPLET_NUM, y = SCORE, shape = TRIPLET)) +
    geom_point(aes(color = TRIPLET)) +
    ylim(sub_df_min, sub_df_max) +
    facet_grid(rows = vars(READ_NUM)) +
    theme_light()

show(sub_df_plot)
