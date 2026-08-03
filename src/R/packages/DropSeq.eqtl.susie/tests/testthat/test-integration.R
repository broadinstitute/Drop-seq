test_that("end-to-end run over the synthetic fixture produces correct schemas, sort order, and files", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())

  res <- run_eqtl_finemapping_files(
    index_eqtl_file = fx$index_eqtl_file, cis_pairs_file = fx$cis_pairs_file,
    expression_file = fx$expression_file, genotype_file = fx$genotype_file,
    covariate_file = fx$covariate_file,
    verbose_outfile = outs$verbose, credible_set_outfile = outs$credible_sets,
    report_outfile = outs$report
  )

  # (1) GENE1 has a reported credible set containing the strong shared-effect variant.
  # Filtering is by gene_name: it is the cross-file join key (matching phenotype_id/pid),
  # while gene_id ("ENSG_GENE1", ...) is a separate, preserved Ensembl-style annotation.
  gene1_report <- res$report[res$report$gene_name == "GENE1", ]
  expect_equal(gene1_report$status, "success")
  expect_equal(gene1_report$gene_id, "ENSG_GENE1")
  expect_gt(gene1_report$n_credible_sets, 0)
  expect_true(fx$shared_variant %in% res$credible_sets$variant_id[res$credible_sets$gene_name == "GENE1"])

  # (2) GENE2 has no signal and therefore no credible set (fallback index-SNP row).
  gene2_report <- res$report[res$report$gene_name == "GENE2", ]
  expect_equal(gene2_report$status, "no_credible_set")
  gene2_slim <- res$credible_sets[res$credible_sets$gene_name == "GENE2", ]
  expect_equal(nrow(gene2_slim), 1)
  expect_true(gene2_slim$is_index_variant)
  expect_false(gene2_slim$in_credible_set)

  # (3) GENE3 has missing dosage values that remove donors.
  gene3_report <- res$report[res$report$gene_name == "GENE3", ]
  expect_equal(gene3_report$n_donors_dropped_missing_genotype, length(fx$missing_donor_idx))
  expect_lt(gene3_report$n_donors_used, gene3_report$n_donors_total)

  # (4) The same variant (shared_variant) occurs for two different genes.
  shared_rows <- res$verbose[res$verbose$variant_id == fx$shared_variant, ]
  expect_setequal(shared_rows$gene_name, c("GENE1", "GENE3"))
  expect_setequal(shared_rows$gene_id, c("ENSG_GENE1", "ENSG_GENE3"))

  # (6) Exact output schemas.
  expect_equal(names(res$verbose), c(
    "dataset", "gene_id", "gene_name", "index_variant_id", "variant_id", "is_index_variant",
    "qval", "index_pval_nominal", "index_slope", "index_slope_se", "cis_start_distance",
    "cis_af", "cis_ma_samples", "cis_ma_count", "cis_pval_nominal", "cis_slope", "cis_slope_se",
    "pip", "posterior_mean", "posterior_sd", "in_credible_set", "credible_set_id",
    "credible_set_requested_coverage", "credible_set_coverage", "credible_set_min_abs_corr",
    "credible_set_mean_abs_corr", "credible_set_median_abs_corr", "credible_set_alpha",
    "susie_converged", "analysis_status", "n_cis_variants", "n_donors_total", "n_donors_used",
    "n_donors_dropped_missing_genotype", "n_covariate_columns", "L", "coverage"
  ))
  expect_equal(names(res$credible_sets), c(
    "dataset", "gene_id", "gene_name", "index_variant_id", "variant_id", "is_index_variant",
    "in_credible_set", "credible_set_id", "pip", "posterior_mean", "posterior_sd",
    "credible_set_requested_coverage", "credible_set_coverage", "credible_set_min_abs_corr",
    "credible_set_mean_abs_corr", "credible_set_median_abs_corr", "credible_set_alpha",
    "qval", "susie_converged",
    "analysis_status", "n_cis_variants", "n_donors_total", "n_donors_used",
    "n_donors_dropped_missing_genotype"
  ))

  # credible_set_alpha is the credible set's own single-effect alpha distribution:
  # it sums to exactly 1 across ALL cis variants tested for the gene (not just the
  # variants inside the credible set itself), unlike the marginal `pip` column,
  # which is not constrained to sum to 1 or to `coverage`.
  gene1_verbose <- res$verbose[res$verbose$gene_name == "GENE1", ]
  expect_equal(sum(gene1_verbose$credible_set_alpha, na.rm = TRUE), 1, tolerance = 1e-6)
  # GENE2 reports no credible set at all, so it has no alpha distribution to report.
  gene2_verbose <- res$verbose[res$verbose$gene_name == "GENE2", ]
  expect_true(all(is.na(gene2_verbose$credible_set_alpha)))
  expect_equal(names(res$report), c(
    "dataset", "gene_id", "gene_name", "index_variant_id", "qval", "status", "n_cis_variants",
    "n_donors_total", "n_donors_used", "n_donors_dropped_missing_genotype", "n_covariate_columns",
    "susie_converged", "n_credible_sets", "n_credible_set_snps", "index_snp_in_credible_set",
    "warning_message", "error_message", "L", "coverage"
  ))

  # Sort order: gene_id, in_credible_set desc, pip desc, variant_id.
  expect_false(is.unsorted(res$verbose$gene_id))
  for (g in unique(res$verbose$gene_id)) {
    sub <- res$verbose[res$verbose$gene_id == g, ]
    # TRUE rows (rank 1) precede FALSE rows (rank 2) precede NA rows (rank 3).
    rank_val <- ifelse(is.na(sub$in_credible_set), 3, ifelse(sub$in_credible_set, 1, 2))
    expect_false(is.unsorted(rank_val))
    for (r in unique(rank_val)) {
      pip_grp <- sub$pip[rank_val == r]
      pip_grp <- pip_grp[!is.na(pip_grp)]
      if (length(pip_grp) > 1) expect_false(is.unsorted(rev(pip_grp)))
    }
  }

  # (7) invisible return and all three files written.
  expect_true(file.exists(outs$verbose))
  expect_true(file.exists(outs$credible_sets))
  expect_true(file.exists(outs$report))
  reread_verbose <- data.table::fread(outs$verbose)
  expect_equal(nrow(reread_verbose), nrow(res$verbose))
})

test_that("one row per gene-SNP pair is enforced in both the verbose and slim outputs", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  res <- run_eqtl_finemapping_files(
    index_eqtl_file = fx$index_eqtl_file, cis_pairs_file = fx$cis_pairs_file,
    expression_file = fx$expression_file, genotype_file = fx$genotype_file,
    covariate_file = fx$covariate_file,
    verbose_outfile = outs$verbose, credible_set_outfile = outs$credible_sets,
    report_outfile = outs$report
  )
  key_v <- paste(res$verbose$gene_id, res$verbose$variant_id)
  expect_equal(anyDuplicated(key_v), 0)
  key_s <- paste(res$credible_sets$gene_id, res$credible_sets$variant_id)
  expect_equal(anyDuplicated(key_s), 0)
})

test_that("an exact covariate subset is honored end to end", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  res <- run_eqtl_finemapping_files(
    index_eqtl_file = fx$index_eqtl_file, cis_pairs_file = fx$cis_pairs_file,
    expression_file = fx$expression_file, genotype_file = fx$genotype_file,
    covariate_file = fx$covariate_file, covariates = c("PC1"),
    verbose_outfile = outs$verbose, credible_set_outfile = outs$credible_sets,
    report_outfile = outs$report
  )
  expect_true(all(res$report$n_covariate_columns == 1))
})

test_that("n_cores > 1 gives identical results to sequential execution, including fatal errors and warnings", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)

  outs_seq <- out_paths(local_temp_dir())
  res_seq <- run_eqtl_finemapping_files(
    index_eqtl_file = fx$index_eqtl_file, cis_pairs_file = fx$cis_pairs_file,
    expression_file = fx$expression_file, genotype_file = fx$genotype_file,
    covariate_file = fx$covariate_file, n_cores = 1,
    verbose_outfile = outs_seq$verbose, credible_set_outfile = outs_seq$credible_sets,
    report_outfile = outs_seq$report
  )

  outs_par <- out_paths(local_temp_dir())
  res_par <- suppressMessages(run_eqtl_finemapping_files(
    index_eqtl_file = fx$index_eqtl_file, cis_pairs_file = fx$cis_pairs_file,
    expression_file = fx$expression_file, genotype_file = fx$genotype_file,
    covariate_file = fx$covariate_file, n_cores = 2,
    verbose_outfile = outs_par$verbose, credible_set_outfile = outs_par$credible_sets,
    report_outfile = outs_par$report
  ))

  expect_equal(res_seq$verbose, res_par$verbose)
  expect_equal(res_seq$credible_sets, res_par$credible_sets)
  expect_equal(res_seq$report, res_par$report)

  # A fatal, run-aborting error (monomorphic variant) must still stop the whole run,
  # in the same deterministic way, regardless of n_cores.
  fx2 <- make_synthetic_dataset(local_temp_dir())
  rewrite_gz_tsv(fx2$genotype_file, function(dt) {
    dt[dt$pid == fx2$gene1_variants[1], fx2$donors] <- 1
    dt
  })
  outs2 <- out_paths(local_temp_dir())
  expect_error(
    suppressMessages(run_eqtl_finemapping_files(
      index_eqtl_file = fx2$index_eqtl_file, cis_pairs_file = fx2$cis_pairs_file,
      expression_file = fx2$expression_file, genotype_file = fx2$genotype_file,
      covariate_file = fx2$covariate_file, n_cores = 2,
      verbose_outfile = outs2$verbose, credible_set_outfile = outs2$credible_sets,
      report_outfile = outs2$report
    )),
    "monomorphic or zero-variance"
  )
})
