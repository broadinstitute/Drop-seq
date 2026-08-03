test_that("a gene-level susieR::susie() error is caught, warned, and does not stop other genes", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())

  original_susie <- susieR::susie
  local_mocked_bindings(
    susie = function(X, y, ...) {
      if (identical(colnames(X)[1], fx$gene1_variants[1])) stop("forced susie failure")
      original_susie(X, y, ...)
    },
    .package = "susieR"
  )

  res <- expect_warning(
    run_eqtl_finemapping_files(
      index_eqtl_file = fx$index_eqtl_file, cis_pairs_file = fx$cis_pairs_file,
      expression_file = fx$expression_file, genotype_file = fx$genotype_file,
      covariate_file = fx$covariate_file,
      verbose_outfile = outs$verbose, credible_set_outfile = outs$credible_sets,
      report_outfile = outs$report
    ),
    "susie_error"
  )

  gene1_report <- res$report[res$report$gene_name == "GENE1", ]
  expect_equal(gene1_report$status, "susie_error")
  expect_match(gene1_report$error_message, "forced susie failure")

  gene1_verbose <- res$verbose[res$verbose$gene_name == "GENE1", ]
  expect_true(all(is.na(gene1_verbose$pip)))
  expect_true(all(is.na(gene1_verbose$in_credible_set)))

  gene1_slim <- res$credible_sets[res$credible_sets$gene_name == "GENE1", ]
  expect_equal(nrow(gene1_slim), 1)
  expect_equal(gene1_slim$variant_id, fx$shared_variant) # the tensorQTL index SNP
  expect_false(gene1_slim$in_credible_set)

  # Other genes are unaffected.
  expect_true("GENE2" %in% res$report$gene_name)
  expect_true("GENE3" %in% res$report$gene_name)
  expect_equal(res$report$status[res$report$gene_name == "GENE3"], "success")
})

test_that("a non-converged fit is treated as a recoverable, warned fallback", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())

  original_susie <- susieR::susie
  local_mocked_bindings(
    susie = function(X, y, ...) {
      fit <- original_susie(X, y, ...)
      fit$converged <- FALSE
      fit
    },
    .package = "susieR"
  )

  res <- expect_warning(
    run_eqtl_finemapping_files(
      index_eqtl_file = fx$index_eqtl_file, cis_pairs_file = fx$cis_pairs_file,
      expression_file = fx$expression_file, genotype_file = fx$genotype_file,
      covariate_file = fx$covariate_file,
      verbose_outfile = outs$verbose, credible_set_outfile = outs$credible_sets,
      report_outfile = outs$report
    ),
    "non_converged"
  )
  expect_true(all(res$report$status == "non_converged"))
  expect_true(all(res$report$susie_converged == FALSE))
})

test_that("a susie_get_cs() error is a recoverable credible_set_error fallback", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())

  # susieR::susie() itself calls susie_get_cs() internally (to populate fit$sets)
  # before our own explicit post-fit call. Fail only every other call so the
  # internal call succeeds (fit is returned) and our explicit call is the one
  # that fails, producing "credible_set_error" rather than "susie_error".
  call_count <- 0
  original_get_cs <- susieR::susie_get_cs
  local_mocked_bindings(
    susie_get_cs = function(...) {
      call_count <<- call_count + 1
      if (call_count %% 2 == 0) stop("forced cs extraction failure")
      original_get_cs(...)
    },
    .package = "susieR"
  )

  res <- expect_warning(
    run_eqtl_finemapping_files(
      index_eqtl_file = fx$index_eqtl_file, cis_pairs_file = fx$cis_pairs_file,
      expression_file = fx$expression_file, genotype_file = fx$genotype_file,
      covariate_file = fx$covariate_file,
      verbose_outfile = outs$verbose, credible_set_outfile = outs$credible_sets,
      report_outfile = outs$report
    ),
    "credible_set_error"
  )
  expect_true(all(res$report$status == "credible_set_error"))
  expect_true(all(grepl("forced cs extraction failure", res$report$error_message)))
  expect_true(all(!res$credible_sets$in_credible_set))
})

test_that("a SNP unexpectedly assigned to multiple credible sets is a recoverable credible_set_error fallback, not a fatal stop", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())

  # susie_get_cs() is called once internally by susie() (to populate fit$sets) and
  # once explicitly by process_one_gene(); target only GENE1's explicit call (the
  # second call for that gene's X) and inject an overlapping credible set so one
  # SNP is assigned to two sets, without disturbing any other gene.
  call_count <- 0
  original_get_cs <- susieR::susie_get_cs
  local_mocked_bindings(
    susie_get_cs = function(res, X = NULL, ...) {
      cs <- original_get_cs(res, X = X, ...)
      if (!is.null(X) && identical(colnames(X)[1], fx$gene1_variants[1])) {
        call_count <<- call_count + 1
        if (call_count %% 2 == 0 && length(cs$cs) > 0) {
          cs$cs$L_dup <- cs$cs[[1]][1]
          cs$coverage <- c(cs$coverage, cs$coverage[1])
          cs$cs_index <- c(cs$cs_index, cs$cs_index[1])
          cs$purity <- rbind(cs$purity, cs$purity[1, ])
          rownames(cs$purity)[nrow(cs$purity)] <- "L_dup"
        }
      }
      cs
    },
    .package = "susieR"
  )

  res <- expect_warning(
    run_eqtl_finemapping_files(
      index_eqtl_file = fx$index_eqtl_file, cis_pairs_file = fx$cis_pairs_file,
      expression_file = fx$expression_file, genotype_file = fx$genotype_file,
      covariate_file = fx$covariate_file,
      verbose_outfile = outs$verbose, credible_set_outfile = outs$credible_sets,
      report_outfile = outs$report
    ),
    "credible_set_error"
  )

  gene1_report <- res$report[res$report$gene_name == "GENE1", ]
  expect_equal(gene1_report$status, "credible_set_error")
  expect_match(gene1_report$error_message, "unexpectedly assigned to multiple credible sets")

  gene1_slim <- res$credible_sets[res$credible_sets$gene_name == "GENE1", ]
  expect_equal(nrow(gene1_slim), 1)
  expect_false(gene1_slim$in_credible_set)

  # Other genes are unaffected and the run completes rather than aborting.
  expect_equal(res$report$status[res$report$gene_name == "GENE3"], "success")
})

test_that("a successfully fit gene with no credible set is not a warning condition", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  res <- expect_silent(
    run_eqtl_finemapping_files(
      index_eqtl_file = fx$index_eqtl_file, cis_pairs_file = fx$cis_pairs_file,
      expression_file = fx$expression_file, genotype_file = fx$genotype_file,
      covariate_file = fx$covariate_file,
      verbose_outfile = outs$verbose, credible_set_outfile = outs$credible_sets,
      report_outfile = outs$report
    )
  )
  gene2_report <- res$report[res$report$gene_name == "GENE2", ]
  expect_equal(gene2_report$status, "no_credible_set")
  expect_true(is.na(gene2_report$warning_message))
})
