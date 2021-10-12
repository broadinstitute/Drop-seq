library(tidyverse)
library(gridExtra)

rxns_path <- '/Users/kshakir/Downloads/rxns.tsv'
qscore_counts <- '/Users/kshakir/Downloads/qscore_counts/'
plots_dir <- '/Users/kshakir/Downloads/'

read_qscore_counts_df <- function(qscore_counts, rxns_df, idx_rxn) {
  qscore_counts_path <- paste0(qscore_counts, idx_rxn, '.tsv')
  qscore_matrix <- as.matrix(read.table(qscore_counts_path, sep = '\t'))
  qscore_matrix <- qscore_matrix[1:450, 1:40, drop = FALSE]
  rownames(qscore_matrix) <- seq(0, nrow(qscore_matrix) - 1)
  colnames(qscore_matrix) <- seq(0, ncol(qscore_matrix) - 1)

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

  qscore_df <-
    qscore_matrix %>%
      as.data.frame() %>%
      rownames_to_column('pos') %>%
      pivot_longer(cols = !'pos', names_to = 'qscore', values_to = 'count') %>%
      mutate_at(vars(pos, qscore), as.integer) %>%
      mutate(
        sample_name = sample_name,
        sequencer = sequencer,
        idx_rxn = idx_rxn,
        fastq2 = fastq2,
        log_10_count = ifelse(count == 0, 0, log10(count))
      )

  return(qscore_df)
}

generate_plot <- function(qscore_counts_df, sequencer_name) {
  qscore_counts_df_sequencer <-
    qscore_counts_df %>%
      filter(sequencer == sequencer_name)

  total_reads <- format(sum(qscore_counts_df_sequencer$count), nsmall = 0, big.mark = ",")

  maxtrix_plot <-
    qscore_counts_df_sequencer %>%
      ggplot(aes(x = pos, y = qscore)) +
      geom_raster(aes(fill = log_10_count)) +
      # scale_y_continuous(label = scales::comma_format(accuracy = 1)) +
      labs(
        title = paste0(sequencer_name, ' 450bp reads - count of quality scores at each read position'),
        subtitle = paste0(total_reads, ' reads'),
        x = 'base position (zero based)',
        y = 'quality score'
      ) +
      theme_bw() +
      theme(
        # panel.grid.major.y = element_blank(),
        # panel.grid.minor.y = element_blank(),
        # panel.border = element_blank(),
        # strip.text.y = element_text(angle = 0)
      )

  ggsave(
    plot = maxtrix_plot,
    filename = paste0(plots_dir, 'qscore_counts_', tolower(sequencer_name), '.png'),
    dpi = 72,
    width = 600,
    height = 400,
    units = 'px'
  )
  # show(maxtrix_plot)
  return(maxtrix_plot)
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

if (!exists('qscore_counts_df')) {
  qscore_counts_df <-
    bind_rows(
      lapply(
        rxn_idx_names,
        read_qscore_counts_df,
        qscore_counts = qscore_counts,
        rxns_df = rxns_df
      )
    )
}

generate_plot(qscore_counts_df, 'NovaSeq')
generate_plot(qscore_counts_df, 'MiSeq')
