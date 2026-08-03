run_fx <- function(fx, outs) {
  run_eqtl_finemapping_files(
    index_eqtl_file = fx$index_eqtl_file, cis_pairs_file = fx$cis_pairs_file,
    expression_file = fx$expression_file, genotype_file = fx$genotype_file,
    covariate_file = fx$covariate_file,
    verbose_outfile = outs$verbose, credible_set_outfile = outs$credible_sets,
    report_outfile = outs$report
  )
}

test_that("a full run over the synthetic fixture succeeds with no fatal errors (sanity baseline)", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  expect_no_error(run_fx(fx, outs))
})

test_that("a donor-set mismatch between expression and genotype is fatal", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  rewrite_gz_tsv(fx$genotype_file, function(dt) {
    dt[[fx$donors[1]]] <- NULL
    dt
  })
  expect_error(run_fx(fx, outs), "Donor sets do not match")
  expect_false(file.exists(outs$verbose))
})

test_that("a missing (NA) expression value for a significant gene is fatal", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  rewrite_gz_tsv(fx$expression_file, function(dt) {
    dt[dt$pid == "GENE1", fx$donors[1]] <- NA
    dt
  })
  expect_error(run_fx(fx, outs), "missing or non-finite expression")
})

test_that("a missing (NA) selected covariate value is fatal", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  rewrite_txt_tsv(fx$covariate_file, function(dt) {
    dt[dt$ID == "PC1", fx$donors[1]] <- NA
    dt
  })
  expect_error(run_fx(fx, outs), "missing or non-finite selected covariate")
})

test_that("an index SNP absent from its gene's cis-pair rows is fatal", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  rewrite_gz_tsv(fx$cis_pairs_file, function(dt) {
    dt[!(dt$phenotype_id == "GENE1" & dt$variant_id == fx$shared_variant), ]
  })
  expect_error(run_fx(fx, outs), "not present among")
})

test_that("a cis variant missing from the genotype BED is fatal", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  rewrite_gz_tsv(fx$genotype_file, function(dt) dt[dt$pid != fx$gene2_variants[1], ])
  expect_error(run_fx(fx, outs), "no row for required cis variant")
})

test_that("a duplicate retained gene_id row in the index eQTL file is fatal", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  rewrite_gz_tsv(fx$index_eqtl_file, function(dt) rbind(dt, dt[dt$gene_id == "ENSG_GENE1", ]))
  expect_error(run_fx(fx, outs), "duplicate retained gene_id")
})

test_that("a duplicate retained gene_name row in the index eQTL file is fatal (gene_name is the join key)", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  rewrite_gz_tsv(fx$index_eqtl_file, function(dt) {
    dup <- dt[dt$gene_name == "GENE1", ]
    dup$gene_id <- "ENSG_GENE1_DUP" # distinct gene_id, same gene_name: still ambiguous for joins
    rbind(dt, dup)
  })
  expect_error(run_fx(fx, outs), "duplicate retained gene_name")
})

test_that("a duplicate variant_id row within one gene's cis-pairs is fatal", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  rewrite_gz_tsv(fx$cis_pairs_file, function(dt) {
    rbind(dt, dt[dt$phenotype_id == "GENE1" & dt$variant_id == fx$gene1_variants[1], ])
  })
  expect_error(run_fx(fx, outs), "duplicate variant_id")
})

test_that("a monomorphic (constant-dosage) cis variant is fatal", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  rewrite_gz_tsv(fx$genotype_file, function(dt) {
    dt[dt$pid == fx$gene1_variants[1], fx$donors] <- 1
    dt
  })
  expect_error(run_fx(fx, outs), "monomorphic or zero-variance")
})

test_that("insufficient retained donors after genotype-missingness filtering is fatal", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir, n_donors = 20)
  outs <- out_paths(local_temp_dir())
  rewrite_gz_tsv(fx$genotype_file, function(dt) {
    # Blank out one cis variant's dosage for all but one donor of GENE1.
    row <- which(dt$pid == fx$gene1_variants[1])
    dt[row, fx$donors[-1]] <- NA
    dt
  })
  expect_error(run_fx(fx, outs), "at least 2 are required")
})

test_that("no output files are written when a fatal validation error occurs before gene-level fitting", {
  dir <- local_temp_dir()
  fx <- make_synthetic_dataset(dir)
  outs <- out_paths(local_temp_dir())
  rewrite_gz_tsv(fx$genotype_file, function(dt) {
    dt[[fx$donors[1]]] <- NULL
    dt
  })
  expect_error(run_fx(fx, outs))
  expect_false(file.exists(outs$verbose))
  expect_false(file.exists(outs$credible_sets))
  expect_false(file.exists(outs$report))
})
