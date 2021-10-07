library(tidyverse)
library(gridExtra)

rxns_path <- '/Users/kshakir/Downloads/rxns.tsv'
repeats_dir <- '/Users/kshakir/Downloads/glutamine_repeats_2021_10_07_a/'
plots_dir <- '/Users/kshakir/Downloads/'

read_glutamine_repeats_df <- function(rxns_df, idx_rxn) {
  glutamine_repeats_path <- paste0(repeats_dir, idx_rxn, '.tsv')

  glutamine_repeats_df <-
    glutamine_repeats_path %>%
      read_delim(
        delim = '\t',
        col_types = cols(
          PREFIX_OFFSET = col_integer(),
          PREFIX_MISMATCHES = col_integer(),
          SUFFIX_MISMATCHES = col_integer(),
          GLUTAMINE_COUNT = col_integer(),
          Q2_START = col_integer(),
          COUNT = col_integer(),
          FIRST = col_integer()
        )
      )

  rxns_df <-
    rxns_path %>%
      read_delim(
        delim = '\t',
        col_types = cols(
          ORDER = col_integer()
        )
      )

  rxn_df <- rxns_df %>%
    filter(
      RXN_IDX_0 == idx_rxn |
        RXN_IDX_1 == idx_rxn |
        RXN_IDX_2 == idx_rxn
    )
  sample_name <- rxn_df$SAMPLE_NAME
  sequencer <- rxn_df$SEQUENCER
  if (!is.na(rxn_df$RXN_IDX_0) && idx_rxn == rxn_df$RXN_IDX_0) {
    fastq2 <- rxn_df$FASTQ2_0
  } else if (!is.na(rxn_df$RXN_IDX_1) && idx_rxn == rxn_df$RXN_IDX_1) {
    fastq2 <- rxn_df$FASTQ2_1
  } else if (!is.na(rxn_df$RXN_IDX_2) && idx_rxn == rxn_df$RXN_IDX_2) {
    fastq2 <- rxn_df$FASTQ2_2
  } else {
    fastq2 <- 'unknown'
  }

  glutamine_repeats_df <-
    glutamine_repeats_df %>%
      mutate(
        sample_name = sample_name,
        sequencer = sequencer,
        idx_rxn = idx_rxn,
        fastq2 = fastq2,
        q2_bin = as.integer((Q2_START %/% 100) * 100),
        match =
          ifelse(
            PREFIX_OFFSET +
              PREFIX_MISMATCHES +
              SUFFIX_MISMATCHES == 0,
            'exact',
            'non-exact'
          )
      )

  glutamine_repeats_df <-
    glutamine_repeats_df %>%
      mutate(
        q2_tail_start =
          ifelse(
            q2_bin < 400,
            paste0(q2_bin, ' - ', as.integer(q2_bin + 99)),
            paste0(q2_bin, '+'))
      )

  return(glutamine_repeats_df)
}

generate_plot <- function(glutamine_repeats_df, sequencer_name, filter_description, color_map) {
  glutamine_repeats_df_sequencer <-
    glutamine_repeats_df %>%
      filter(sequencer == sequencer_name)

  total_reads <- format(sum(glutamine_repeats_df_sequencer$COUNT), nsmall = 0, big.mark = ",")

  hist_plot <-
    glutamine_repeats_df_sequencer %>%
      ggplot(aes(x = GLUTAMINE_COUNT)) +
      geom_histogram(aes(weight = COUNT, fill = q2_tail_start), binwidth = 1) +
      scale_y_continuous(label = scales::comma_format(accuracy = 1)) +
      scale_fill_manual(values = color_map) +
      labs(
        title = paste0(sequencer_name, ' 450bp reads with glutamine repeats (', filter_description, ')'),
        subtitle = paste0(total_reads, ' reads; CAG or CAA repeats between CCTTCGAGTCCCTCAAGTCCTTC and CCGCCACCG'),
        x = 'num glutamine repeats',
        y = 'num reads'
      ) +
      theme_bw() +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        strip.text.y = element_text(angle = 0)
      ) +
      facet_grid(rows = vars(match), scales = 'free', space = 'free')

  ggsave(
    plot = hist_plot,
    filename = paste0(plots_dir, 'glutamine_repeats_', tolower(sequencer_name), '_', filter_description, '.png'),
    dpi = 72,
    width = 1000,
    height = 1000,
    units = 'px'
  )
  # show(hist_plot)
  return(hist_plot)
}

rxns_df <-
  rxns_path %>%
    read_delim(
      delim = '\t',
      col_types = cols(
        ORDER = col_integer()
      )
    )

rxn_idx_names <-
  rxns_df %>%
    select(starts_with('RXN_IDX_')) %>%
    pivot_longer(cols = everything(), values_drop_na = TRUE) %>%
    pull(value)

if (!exists('glutamine_repeats_df')) {
  glutamine_repeats_df <-
    bind_rows(
      lapply(
        rxn_idx_names,
        read_glutamine_repeats_df,
        rxns_df = rxns_df
      )
    )

  glutamine_repeats_df$q2_tail_start <- as.factor(glutamine_repeats_df$q2_tail_start)
}

color_labels <- levels(glutamine_repeats_df$q2_tail_start)
color_map <- hsv(1, 1, seq(0, 1, length.out = length(color_labels)))
names(color_map) <- rev(color_labels)

glutamine_repeats_df_tail <-
  glutamine_repeats_df %>%
    filter(GLUTAMINE_COUNT > 85)

generate_plot(glutamine_repeats_df, 'NovaSeq', 'full', color_map)
generate_plot(glutamine_repeats_df_tail, 'NovaSeq', 'tail', color_map)
generate_plot(glutamine_repeats_df, 'MiSeq', 'full', color_map)
generate_plot(glutamine_repeats_df_tail, 'MiSeq', 'tail', color_map)
