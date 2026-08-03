# Internal output-table builders and writer for run_eqtl_finemapping_files().
#
# All functions in this file are internal (not exported).

#' @keywords internal
#' @noRd
get_scalar <- function(row, col, default = NA) {
  if (col %in% names(row)) row[[col]][1] else default
}

#' @keywords internal
#' @noRd
get_opt_column <- function(df, col, n, default = NA) {
  if (col %in% names(df)) df[[col]] else rep(default, n)
}

#' @keywords internal
#' @noRd
build_gene_verbose_rows <- function(dataset, gene_row, gene_cis_pairs, variant_ids_gene,
                                     analysis_status, snp_res, n_donors_total, n_donors_used,
                                     n_donors_dropped_missing_genotype, n_covariate_columns,
                                     susie_converged, L, coverage) {
  n_var <- length(variant_ids_gene)
  gene_id <- gene_row$gene_id
  gene_name <- gene_row$gene_name
  index_variant_id <- gene_row$variant_id
  qval <- gene_row$qval

  cis_match <- match(variant_ids_gene, gene_cis_pairs$variant_id)

  if (is.null(snp_res)) {
    pip <- rep(NA_real_, n_var)
    posterior_mean <- rep(NA_real_, n_var)
    posterior_sd <- rep(NA_real_, n_var)
    in_credible_set <- rep(NA, n_var)
    credible_set_id <- rep(NA_character_, n_var)
    credible_set_coverage <- rep(NA_real_, n_var)
    cs_min <- rep(NA_real_, n_var)
    cs_mean <- rep(NA_real_, n_var)
    cs_median <- rep(NA_real_, n_var)
    cs_alpha <- rep(NA_real_, n_var)
  } else {
    ord <- match(variant_ids_gene, snp_res$variant_id)
    pip <- snp_res$pip[ord]
    posterior_mean <- snp_res$posterior_mean[ord]
    posterior_sd <- snp_res$posterior_sd[ord]
    in_credible_set <- snp_res$in_credible_set[ord]
    credible_set_id <- snp_res$credible_set_id[ord]
    credible_set_coverage <- snp_res$credible_set_coverage[ord]
    cs_min <- snp_res$credible_set_min_abs_corr[ord]
    cs_mean <- snp_res$credible_set_mean_abs_corr[ord]
    cs_median <- snp_res$credible_set_median_abs_corr[ord]
    cs_alpha <- snp_res$credible_set_alpha[ord]
  }

  data.frame(
    dataset = rep(dataset, n_var),
    gene_id = rep(gene_id, n_var),
    gene_name = rep(gene_name, n_var),
    index_variant_id = rep(index_variant_id, n_var),
    variant_id = variant_ids_gene,
    is_index_variant = variant_ids_gene == index_variant_id,
    qval = rep(qval, n_var),
    index_pval_nominal = rep(get_scalar(gene_row, "pval_nominal"), n_var),
    index_slope = rep(get_scalar(gene_row, "slope"), n_var),
    index_slope_se = rep(get_scalar(gene_row, "slope_se"), n_var),
    cis_start_distance = get_opt_column(gene_cis_pairs, "start_distance", nrow(gene_cis_pairs))[cis_match],
    cis_af = get_opt_column(gene_cis_pairs, "af", nrow(gene_cis_pairs))[cis_match],
    cis_ma_samples = get_opt_column(gene_cis_pairs, "ma_samples", nrow(gene_cis_pairs))[cis_match],
    cis_ma_count = get_opt_column(gene_cis_pairs, "ma_count", nrow(gene_cis_pairs))[cis_match],
    cis_pval_nominal = get_opt_column(gene_cis_pairs, "pval_nominal", nrow(gene_cis_pairs))[cis_match],
    cis_slope = get_opt_column(gene_cis_pairs, "slope", nrow(gene_cis_pairs))[cis_match],
    cis_slope_se = get_opt_column(gene_cis_pairs, "slope_se", nrow(gene_cis_pairs))[cis_match],
    pip = pip,
    posterior_mean = posterior_mean,
    posterior_sd = posterior_sd,
    in_credible_set = in_credible_set,
    credible_set_id = credible_set_id,
    credible_set_requested_coverage = rep(coverage, n_var),
    credible_set_coverage = credible_set_coverage,
    credible_set_min_abs_corr = cs_min,
    credible_set_mean_abs_corr = cs_mean,
    credible_set_median_abs_corr = cs_median,
    credible_set_alpha = cs_alpha,
    susie_converged = rep(susie_converged, n_var),
    analysis_status = rep(analysis_status, n_var),
    n_cis_variants = rep(n_var, n_var),
    n_donors_total = rep(n_donors_total, n_var),
    n_donors_used = rep(n_donors_used, n_var),
    n_donors_dropped_missing_genotype = rep(n_donors_dropped_missing_genotype, n_var),
    n_covariate_columns = rep(n_covariate_columns, n_var),
    L = rep(L, n_var),
    coverage = rep(coverage, n_var),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
build_gene_credible_set_rows <- function(dataset, gene_row, verbose_rows, analysis_status) {
  index_variant_id <- gene_row$variant_id

  slim_cols <- c(
    "dataset", "gene_id", "gene_name", "index_variant_id", "variant_id",
    "is_index_variant", "in_credible_set", "credible_set_id", "pip",
    "posterior_mean", "posterior_sd", "credible_set_requested_coverage",
    "credible_set_coverage", "credible_set_min_abs_corr",
    "credible_set_mean_abs_corr", "credible_set_median_abs_corr",
    "credible_set_alpha", "qval",
    "susie_converged", "analysis_status", "n_cis_variants", "n_donors_total",
    "n_donors_used", "n_donors_dropped_missing_genotype"
  )

  if (identical(analysis_status, "success")) {
    sel <- verbose_rows[!is.na(verbose_rows$in_credible_set) & verbose_rows$in_credible_set, ]
    out <- sel[, slim_cols]
  } else {
    idx_row <- verbose_rows[verbose_rows$variant_id == index_variant_id, ]
    idx_row$in_credible_set <- FALSE
    idx_row$credible_set_id <- NA_character_
    idx_row$credible_set_coverage <- NA_real_
    idx_row$credible_set_min_abs_corr <- NA_real_
    idx_row$credible_set_mean_abs_corr <- NA_real_
    idx_row$credible_set_median_abs_corr <- NA_real_
    idx_row$credible_set_alpha <- NA_real_
    idx_row$analysis_status <- analysis_status
    out <- idx_row[, slim_cols]
  }
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
build_gene_report_row <- function(dataset, gene_row, analysis_status, n_cis_variants,
                                   n_donors_total, n_donors_used, n_donors_dropped_missing_genotype,
                                   n_covariate_columns, susie_converged, n_credible_sets,
                                   n_credible_set_snps, index_snp_in_credible_set,
                                   warning_message, error_message, L, coverage) {
  data.frame(
    dataset = dataset,
    gene_id = gene_row$gene_id,
    gene_name = gene_row$gene_name,
    index_variant_id = gene_row$variant_id,
    qval = gene_row$qval,
    status = analysis_status,
    n_cis_variants = n_cis_variants,
    n_donors_total = n_donors_total,
    n_donors_used = n_donors_used,
    n_donors_dropped_missing_genotype = n_donors_dropped_missing_genotype,
    n_covariate_columns = n_covariate_columns,
    susie_converged = susie_converged,
    n_credible_sets = n_credible_sets,
    n_credible_set_snps = n_credible_set_snps,
    index_snp_in_credible_set = index_snp_in_credible_set,
    warning_message = warning_message,
    error_message = error_message,
    L = L,
    coverage = coverage,
    stringsAsFactors = FALSE
  )
}

# Merges each gene's num_var warning (if any) into its report row's
# warning_message, emits a warning() per gene that has one (in gene order, so
# output is deterministic regardless of n_cores), and row-binds every gene's
# verbose/credible-set/report rows into the three final output data.frames.
#' @keywords internal
#' @noRd
combine_gene_results <- function(gene_results, gene_names, num_var_warnings) {
  n_genes <- length(gene_results)
  verbose_list <- vector("list", n_genes)
  slim_list <- vector("list", n_genes)
  report_list <- vector("list", n_genes)

  for (i in seq_len(n_genes)) {
    gene_name <- gene_names[i]
    result_i <- gene_results[[i]]
    report_row <- result_i$report

    warning_message <- report_row$warning_message
    if (gene_name %in% names(num_var_warnings)) {
      nv_msg <- num_var_warnings[[gene_name]]
      warning_message <- if (is.na(warning_message)) nv_msg else paste(warning_message, nv_msg, sep = "; ")
      report_row$warning_message <- warning_message
    }
    if (!is.na(warning_message)) {
      warning(warning_message, call. = FALSE)
    }

    verbose_list[[i]] <- result_i$verbose
    slim_list[[i]] <- result_i$credible_sets
    report_list[[i]] <- report_row
  }

  list(
    verbose = as.data.frame(data.table::rbindlist(verbose_list, fill = TRUE), stringsAsFactors = FALSE),
    credible_sets = as.data.frame(data.table::rbindlist(slim_list, fill = TRUE), stringsAsFactors = FALSE),
    report = as.data.frame(data.table::rbindlist(report_list, fill = TRUE), stringsAsFactors = FALSE)
  )
}

#' @keywords internal
#' @noRd
validate_output_keys <- function(df, gene_label = NULL) {
  key <- paste(df$gene_id, df$variant_id, sep = "")
  if (any(duplicated(key))) {
    stop(sprintf(
      "Unexpected duplicate gene_id/variant_id row(s) in output%s.",
      if (!is.null(gene_label)) sprintf(" for gene '%s'", gene_label) else ""
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
write_finemapping_outputs <- function(verbose_output, credible_set_output, report_output,
                                       verbose_outfile, credible_set_outfile, report_outfile) {
  data.table::fwrite(verbose_output, verbose_outfile, sep = "\t")
  data.table::fwrite(credible_set_output, credible_set_outfile, sep = "\t")
  data.table::fwrite(report_output, report_outfile, sep = "\t")
  invisible(TRUE)
}
