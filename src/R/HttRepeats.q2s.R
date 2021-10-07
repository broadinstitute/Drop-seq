library(tidyverse)

input_path <- '/Users/kshakir/Downloads/htt_repeats.q2s.txt'

df <-
  input_path %>%
    read_delim(
      delim = '\t',
      col_names = 'COUNT',
      col_types = cols(
        COUNT = col_integer()
      )
    )
df$Q2_TAIL_LENGTH <- as.integer(as.integer(row.names(df)) - 1)
df <- df %>% filter(Q2_TAIL_LENGTH <= 250)

count_no_tail <- format(df$COUNT[1], big.mark = ',')

hist_plot <-
  df %>%
    ggplot(aes(x = Q2_TAIL_LENGTH, y = COUNT)) +
    geom_point() +
    scale_y_log10() +
    labs(
      title = 'Count of consecutive Q2 scores on the tail of the first 10,000 reads',
      subtitle = paste0('250bp reads; ', count_no_tail, ' containing no Q2 tail'),
      caption = '/broad/mccarroll_huntingtons_storage/Nora/MiSeqs/CAG_MiSeq_3p_sc_HOXy_091621_JVYG3/get.broadinstitute.org/pkgs/SN0233471/1_JVYG3.1.AGGCAGAA.unmapped.2.fastq.gz'
    ) +
    theme_light()

show(hist_plot)
