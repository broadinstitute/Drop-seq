test_that("run_eqtl_finemapping derives dataset from the directory basename and delegates to run_eqtl_finemapping_files", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs_a <- out_paths(local_temp_dir())
  outs_b <- out_paths(local_temp_dir())

  set.seed(99)
  res_a <- run_eqtl_finemapping(
    dataset_dir = fx$dataset_dir,
    verbose_outfile = outs_a$verbose, credible_set_outfile = outs_a$credible_sets,
    report_outfile = outs_a$report
  )
  set.seed(99)
  res_b <- run_eqtl_finemapping_files(
    index_eqtl_file = fx$index_eqtl_file, cis_pairs_file = fx$cis_pairs_file,
    expression_file = fx$expression_file, genotype_file = fx$genotype_file,
    covariate_file = fx$covariate_file,
    verbose_outfile = outs_b$verbose, credible_set_outfile = outs_b$credible_sets,
    report_outfile = outs_b$report, dataset = fx$dataset
  )

  expect_equal(res_a$verbose, res_b$verbose)
  expect_equal(res_a$credible_sets, res_b$credible_sets)
  expect_equal(res_a$report, res_b$report)
  expect_true(all(res_a$report$dataset == fx$dataset))
})

test_that("run_eqtl_finemapping stops early (before delegating) if a standard file is missing", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  file.remove(fx$covariate_file)
  outs <- out_paths(local_temp_dir())
  expect_error(
    run_eqtl_finemapping(
      dataset_dir = fx$dataset_dir,
      verbose_outfile = outs$verbose, credible_set_outfile = outs$credible_sets,
      report_outfile = outs$report
    ),
    "covariate_file"
  )
  expect_false(file.exists(outs$verbose))
})

test_that("both exported functions return an invisible list of three base data.frames", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  res <- withVisible(run_eqtl_finemapping(
    dataset_dir = fx$dataset_dir,
    verbose_outfile = outs$verbose, credible_set_outfile = outs$credible_sets,
    report_outfile = outs$report
  ))
  expect_false(res$visible)
  expect_named(res$value, c("verbose", "credible_sets", "report"))
  expect_true(is.data.frame(res$value$verbose) && !inherits(res$value$verbose, "data.table"))
  expect_true(is.data.frame(res$value$credible_sets) && !inherits(res$value$credible_sets, "data.table"))
  expect_true(is.data.frame(res$value$report) && !inherits(res$value$report, "data.table"))
})
