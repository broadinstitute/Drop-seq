library(tidyverse)
library(gridExtra)

input_path <- '/Users/kshakir/Downloads/htt_repeats.glutamine.txt'

df <-
  input_path %>%
    read_delim(
      delim = '\t',
      col_names = 'NUM_READS',
      col_types = cols(
        NUM_READS = col_integer()
      )
    )
df$COUNT_GLUTAMINE_REPEATS <- as.integer(as.integer(row.names(df)) - 1)
df <- df %>% filter(COUNT_GLUTAMINE_REPEATS <= 75)

num_observances <- format(sum(df$NUM_READS), nsmall = 0, big.mark = ",")

hist_plot <-
  df %>%
    ggplot(aes(x = COUNT_GLUTAMINE_REPEATS, y = NUM_READS)) +
    geom_col() +
    labs(
      title = paste0('Count of glutamine repeats in ', num_observances, ' of 10,000 reads'),
      subtitle = '250bp reads; reads selected that aligned CAG/CAA repeats between CCTTCGAGTCCCTCAAGTCCTTC and CCGCCACCG',
      caption = '/broad/mccarroll_huntingtons_storage/Nora/MiSeqs/CAG_MiSeq_3p_sc_HOXy_091621_JVYG3/get.broadinstitute.org/pkgs/SN0233471/1_JVYG3.1.AGGCAGAA.unmapped.2.fastq.gz'
    ) +
    theme_light()

show(hist_plot)
