# Internal per-gene analysis helpers: donor/covariate assembly, QR
# residualization, and susieR fitting / credible-set extraction.
#
# All functions in this file are internal (not exported).

#' @keywords internal
#' @noRd
build_gene_data <- function(gene_id, variant_ids, expression_matrix, genotype_matrix,
                             covariate_matrix, donor_order) {
  expression_full <- expression_matrix[gene_id, donor_order]
  genotype_full <- t(genotype_matrix[variant_ids, donor_order, drop = FALSE])
  missing_mask <- apply(genotype_full, 1, function(row) any(is.na(row)))
  n_donors_total <- length(donor_order)
  n_donors_dropped_missing_genotype <- sum(missing_mask)
  n_donors_used <- n_donors_total - n_donors_dropped_missing_genotype

  keep <- donor_order[!missing_mask]
  genotype <- genotype_full[keep, , drop = FALSE]
  expression <- expression_full[keep]
  if (nrow(covariate_matrix) > 0) {
    covariates <- t(covariate_matrix[, keep, drop = FALSE])
  } else {
    covariates <- matrix(numeric(0), nrow = length(keep), ncol = 0)
  }

  list(
    genotype = genotype, expression = expression, covariates = covariates,
    n_donors_total = n_donors_total,
    n_donors_used = n_donors_used,
    n_donors_dropped_missing_genotype = n_donors_dropped_missing_genotype
  )
}

#' @keywords internal
#' @noRd
validate_gene_inputs <- function(gene_data, gene_id, index_variant_id) {
  if (gene_data$n_donors_used < 2) {
    stop(sprintf(
      "Gene '%s' (index SNP '%s'): only %d donor(s) remain after removing donors with missing genotype dosage; at least 2 are required.",
      gene_id, index_variant_id, gene_data$n_donors_used
    ), call. = FALSE)
  }
  if (any(is.na(gene_data$genotype))) {
    stop(sprintf("Gene '%s': unexpected missing genotype value(s) remain after donor filtering.",
                 gene_id), call. = FALSE)
  }
  variances <- apply(gene_data$genotype, 2, stats::var)
  zero_var <- colnames(gene_data$genotype)[!is.finite(variances) | variances <= 0]
  if (length(zero_var) > 0) {
    stop(sprintf("Gene '%s': monomorphic or zero-variance cis variant(s) detected: %s.",
                 gene_id, paste(zero_var, collapse = ", ")), call. = FALSE)
  }
  if (any(!is.finite(gene_data$expression))) {
    stop(sprintf("Gene '%s': missing or non-finite expression value(s) remain after donor filtering.",
                 gene_id), call. = FALSE)
  }
  if (ncol(gene_data$covariates) > 0 && any(!is.finite(gene_data$covariates))) {
    stop(sprintf("Gene '%s': missing or non-finite covariate value(s) remain after donor filtering.",
                 gene_id), call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
residualize_against_covariates <- function(X, y, Z, gene_id) {
  n <- nrow(X)
  Z_design <- cbind(Intercept = rep(1, n), Z)
  qr_z <- qr(Z_design, LAPACK = FALSE)
  y_resid <- as.numeric(qr.resid(qr_z, y))
  X_resid <- qr.resid(qr_z, X)
  names(y_resid) <- rownames(X)
  rownames(X_resid) <- rownames(X)
  colnames(X_resid) <- colnames(X)

  # A vector that lies entirely in the column space of Z_design (e.g. a constant
  # response regressed against an intercept) residualizes to values that are zero
  # only up to floating-point error, not exactly 0. Compare each residual's
  # variance against a tolerance scaled to the corresponding pre-residualization
  # variable's own magnitude, rather than testing for exact non-positivity.
  zero_tol <- sqrt(.Machine$double.eps)
  y_scale <- max(mean(y^2), 1)
  y_var <- stats::var(y_resid)
  if (!is.finite(y_var) || y_var <= zero_tol * y_scale) {
    stop(sprintf("Gene '%s': residualized expression has zero or non-finite variance.", gene_id),
         call. = FALSE)
  }
  col_var <- apply(X_resid, 2, stats::var)
  col_scale <- pmax(apply(X, 2, function(col) mean(col^2)), 1)
  bad <- colnames(X_resid)[!is.finite(col_var) | col_var <= zero_tol * col_scale]
  if (length(bad) > 0) {
    stop(sprintf(
      "Gene '%s': residualized genotype has zero or non-finite variance for variant(s): %s.",
      gene_id, paste(bad, collapse = ", ")
    ), call. = FALSE)
  }
  list(X = X_resid, y = y_resid)
}

#' @keywords internal
#' @noRd
extract_susie_snp_results <- function(fit, cs, variant_ids, gene_id) {
  n_var <- length(variant_ids)
  pip <- as.numeric(fit$pip)
  posterior_mean <- as.numeric(susieR::susie_get_posterior_mean(fit))
  posterior_sd <- as.numeric(susieR::susie_get_posterior_sd(fit))

  in_credible_set <- rep(FALSE, n_var)
  credible_set_id <- rep(NA_character_, n_var)
  cs_coverage <- rep(NA_real_, n_var)
  cs_min <- rep(NA_real_, n_var)
  cs_mean <- rep(NA_real_, n_var)
  cs_median <- rep(NA_real_, n_var)
  cs_alpha <- rep(NA_real_, n_var)

  cs_list <- cs$cs
  if (!is.null(cs_list) && length(cs_list) > 0) {
    cs_names <- names(cs_list)
    membership_count <- integer(n_var)
    for (i in seq_along(cs_list)) {
      idx <- cs_list[[i]]
      membership_count[idx] <- membership_count[idx] + 1L
      in_credible_set[idx] <- TRUE
      credible_set_id[idx] <- cs_names[i]
      cs_coverage[idx] <- cs$coverage[i]
      cs_min[idx] <- cs$purity[cs_names[i], "min.abs.corr"]
      cs_mean[idx] <- cs$purity[cs_names[i], "mean.abs.corr"]
      cs_median[idx] <- cs$purity[cs_names[i], "median.abs.corr"]
      cs_alpha[idx] <- fit$alpha[cs$cs_index[i], idx]
    }
    if (any(membership_count > 1)) {
      bad_variants <- variant_ids[membership_count > 1]
      stop(sprintf(
        "Gene '%s': variant(s) unexpectedly assigned to multiple credible sets: %s.",
        gene_id, paste(bad_variants, collapse = ", ")
      ), call. = FALSE)
    }
    # fit$alpha[l, ] is single-effect component l's own posterior inclusion distribution
    # over ALL cis variants for this gene, and sums to exactly 1 across them (unlike the
    # marginal `pip`, which unions evidence across every active component and need not
    # sum to 1, or to `coverage`, across a credible set's members). When the gene has
    # exactly one credible set, that component is this gene's only source of PIP mass, so
    # its alpha is reported for every cis variant (not just the set's members), making
    # sum(credible_set_alpha) == 1 across the gene. With multiple credible sets, a
    # non-member variant's mass is split across several components with no single
    # attributable value, so alpha is left NA outside each set's own members, and
    # sum(credible_set_alpha) within one set's members approximates its `coverage`
    # instead.
    if (length(cs_list) == 1) {
      cs_alpha <- fit$alpha[cs$cs_index[1], ]
    }
  }

  data.frame(
    variant_id = variant_ids,
    pip = pip,
    posterior_mean = posterior_mean,
    posterior_sd = posterior_sd,
    in_credible_set = in_credible_set,
    credible_set_id = credible_set_id,
    credible_set_coverage = cs_coverage,
    credible_set_min_abs_corr = cs_min,
    credible_set_mean_abs_corr = cs_mean,
    credible_set_median_abs_corr = cs_median,
    credible_set_alpha = cs_alpha,
    stringsAsFactors = FALSE
  )
}

# Runs the full per-gene analysis (data assembly, residualization, SuSiE fit,
# credible-set extraction, and all three output-row builders) for one gene and
# returns its results as data. This function does not call warning() itself:
# a gene-level SuSiE/credible-set failure is recorded only in the returned
# report row's warning_message, so that the caller can emit warnings in a
# fixed, gene-index order after collecting every gene's result. That is what
# keeps warning output deterministic when genes are processed in parallel
# (worker completion order is not guaranteed to match gene order).
#' @keywords internal
#' @noRd
process_one_gene <- function(gene_row, gene_cis_pairs, expression_matrix, genotype_matrix,
                              covariate_matrix, donor_order, dataset, L, coverage) {
  gene_name <- gene_row$gene_name
  index_variant_id <- gene_row$variant_id
  variant_ids_gene <- gene_cis_pairs$variant_id

  gene_data <- build_gene_data(gene_name, variant_ids_gene, expression_matrix,
                                genotype_matrix, covariate_matrix, donor_order)
  validate_gene_inputs(gene_data, gene_name, index_variant_id)
  resid <- residualize_against_covariates(gene_data$genotype, gene_data$expression,
                                           gene_data$covariates, gene_name)

  fit <- NULL
  cs <- NULL
  snp_res <- NULL
  analysis_status <- NA_character_
  error_message <- NA_character_

  fit_or_error <- tryCatch(
    susieR::susie(resid$X, resid$y, L = L, coverage = coverage, verbose = FALSE),
    error = function(e) e
  )

  if (inherits(fit_or_error, "error")) {
    analysis_status <- "susie_error"
    error_message <- conditionMessage(fit_or_error)
  } else {
    fit <- fit_or_error
    if (!is.null(fit$converged) && !isTRUE(fit$converged)) {
      analysis_status <- "non_converged"
    } else {
      cs_or_error <- tryCatch(
        susieR::susie_get_cs(fit, X = resid$X, coverage = coverage),
        error = function(e) e
      )
      if (inherits(cs_or_error, "error")) {
        analysis_status <- "credible_set_error"
        error_message <- conditionMessage(cs_or_error)
      } else {
        cs <- cs_or_error
        snp_res_or_error <- tryCatch(
          extract_susie_snp_results(fit, cs, variant_ids_gene, gene_name),
          error = function(e) e
        )
        if (inherits(snp_res_or_error, "error")) {
          analysis_status <- "credible_set_error"
          error_message <- conditionMessage(snp_res_or_error)
        } else {
          snp_res <- snp_res_or_error
          n_cs <- if (is.null(cs$cs)) 0L else length(cs$cs)
          analysis_status <- if (n_cs > 0) "success" else "no_credible_set"
        }
      }
    }
  }

  fit_warning_message <- NA_character_
  if (analysis_status %in% c("susie_error", "non_converged", "credible_set_error")) {
    fit_warning_message <- sprintf(
      "Gene '%s' (index SNP '%s'): %s%s", gene_name, index_variant_id, analysis_status,
      if (!is.na(error_message)) paste0(" - ", error_message) else ""
    )
  }

  susie_converged_value <- if (!is.null(fit) && !is.null(fit$converged)) isTRUE(fit$converged) else NA

  verbose_rows <- build_gene_verbose_rows(
    dataset = dataset, gene_row = gene_row, gene_cis_pairs = gene_cis_pairs,
    variant_ids_gene = variant_ids_gene, analysis_status = analysis_status, snp_res = snp_res,
    n_donors_total = gene_data$n_donors_total, n_donors_used = gene_data$n_donors_used,
    n_donors_dropped_missing_genotype = gene_data$n_donors_dropped_missing_genotype,
    n_covariate_columns = ncol(gene_data$covariates), susie_converged = susie_converged_value,
    L = L, coverage = coverage
  )
  validate_output_keys(verbose_rows, gene_name)

  slim_rows <- build_gene_credible_set_rows(dataset, gene_row, verbose_rows, analysis_status)

  n_credible_sets <- NA_integer_
  n_credible_set_snps <- NA_integer_
  index_snp_in_credible_set <- NA

  if (identical(analysis_status, "success")) {
    n_credible_sets <- length(cs$cs)
    n_credible_set_snps <- sum(verbose_rows$in_credible_set, na.rm = TRUE)
    index_snp_in_credible_set <- isTRUE(verbose_rows$in_credible_set[verbose_rows$is_index_variant])
  } else if (identical(analysis_status, "no_credible_set")) {
    n_credible_sets <- 0L
    n_credible_set_snps <- 0L
    index_snp_in_credible_set <- FALSE
  }

  report_row <- build_gene_report_row(
    dataset = dataset, gene_row = gene_row, analysis_status = analysis_status,
    n_cis_variants = length(variant_ids_gene), n_donors_total = gene_data$n_donors_total,
    n_donors_used = gene_data$n_donors_used,
    n_donors_dropped_missing_genotype = gene_data$n_donors_dropped_missing_genotype,
    n_covariate_columns = ncol(gene_data$covariates), susie_converged = susie_converged_value,
    n_credible_sets = n_credible_sets, n_credible_set_snps = n_credible_set_snps,
    index_snp_in_credible_set = index_snp_in_credible_set,
    warning_message = fit_warning_message, error_message = error_message, L = L, coverage = coverage
  )

  list(verbose = verbose_rows, credible_sets = slim_rows, report = report_row)
}
