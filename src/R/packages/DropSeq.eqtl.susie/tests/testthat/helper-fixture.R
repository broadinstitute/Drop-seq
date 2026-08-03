# Shared synthetic-dataset fixture used across tests. Not part of the
# installed package; sourced automatically by testthat.

local_temp_dir <- function() {
  d <- tempfile("dseqtlsusie-test-")
  dir.create(d, recursive = TRUE)
  d
}

make_synthetic_dataset <- function(dir, dataset = "TESTCT__TESTREGION", n_donors = 20,
                                    seed = 1, include_num_var = TRUE) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  if (!is.null(old_seed)) {
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
  }

  ds_dir <- file.path(dir, dataset)
  dir.create(ds_dir, showWarnings = FALSE, recursive = TRUE)
  donors <- paste0("D", seq_len(n_donors))

  pc1 <- stats::rnorm(n_donors)
  batch <- rep(c(0, 1), length.out = n_donors)
  cov_df <- data.frame(ID = c("PC1", "Batch1", "Batch2"), rbind(pc1, batch, 1 - batch))
  names(cov_df) <- c("ID", donors)
  utils::write.table(cov_df, file.path(ds_dir, paste0(dataset, ".covariates_peer.txt")),
                      sep = "\t", quote = FALSE, row.names = FALSE)

  make_dosage <- function(n, p = 0.3) stats::rbinom(n, 2, p) + stats::rnorm(n, 0, 0.05)

  gene1_variants <- paste0("V1_", seq_len(6))
  gene2_variants <- paste0("V2_", seq_len(5))
  gene3_variants <- paste0("V3_", seq_len(5))
  shared_variant <- "SNP_SHARED"
  gene1_variants[length(gene1_variants)] <- shared_variant
  gene3_variants[length(gene3_variants)] <- shared_variant

  all_variants <- unique(c(gene1_variants, gene2_variants, gene3_variants))
  geno_mat <- matrix(NA_real_, nrow = length(all_variants), ncol = n_donors,
                      dimnames = list(all_variants, donors))
  for (v in all_variants) geno_mat[v, ] <- make_dosage(n_donors)

  missing_donor_idx <- if (n_donors >= 6) c(2, 4, 6) else integer(0)
  geno_mat["V3_2", missing_donor_idx] <- NA_real_

  geno_bed <- data.frame(
    `#chr` = "chr1", start = seq_along(all_variants) * 100,
    end = seq_along(all_variants) * 100 + 1, pid = all_variants,
    geno_mat, check.names = FALSE
  )
  names(geno_bed) <- c("#chr", "start", "end", "pid", donors)
  gz <- gzfile(file.path(ds_dir, paste0(dataset, ".genotype_matrix.bed.gz")), "w")
  utils::write.table(geno_bed, gz, sep = "\t", quote = FALSE, row.names = FALSE)
  close(gz)

  y1 <- 3 * geno_mat[shared_variant, ] + 0.5 * pc1 + stats::rnorm(n_donors, 0, 0.3)
  y2 <- stats::rnorm(n_donors)
  y3 <- 2 * geno_mat["V3_1", ] + stats::rnorm(n_donors, 0, 0.3)
  expr_mat <- rbind(GENE1 = y1, GENE2 = y2, GENE3 = y3)
  expr_bed <- data.frame(
    `#chr` = "chr1", start = c(1, 2, 3), end = c(2, 3, 4),
    pid = c("GENE1", "GENE2", "GENE3"), expr_mat, check.names = FALSE
  )
  names(expr_bed) <- c("#chr", "start", "end", "pid", donors)
  gz2 <- gzfile(file.path(ds_dir, paste0(dataset, ".gene_expression_normalized.bed.gz")), "w")
  utils::write.table(expr_bed, gz2, sep = "\t", quote = FALSE, row.names = FALSE)
  close(gz2)

  build_pairs <- function(gene, variants) {
    data.frame(
      phenotype_id = gene, variant_id = variants,
      start_distance = seq_along(variants) * 10,
      af = stats::runif(length(variants), 0.05, 0.5),
      ma_samples = sample(5:20, length(variants), replace = TRUE),
      ma_count = sample(5:30, length(variants), replace = TRUE),
      pval_nominal = stats::runif(length(variants), 1e-8, 0.5),
      slope = stats::rnorm(length(variants)),
      slope_se = stats::runif(length(variants), 0.1, 0.5)
    )
  }
  pairs <- rbind(
    build_pairs("GENE1", gene1_variants),
    build_pairs("GENE2", gene2_variants),
    build_pairs("GENE3", gene3_variants),
    build_pairs("GENE4_NOT_SIG", c("V4_1", "V4_2"))
  )
  gz3 <- gzfile(file.path(ds_dir, paste0(dataset, ".cis_qtl_pairs.txt.gz")), "w")
  utils::write.table(pairs, gz3, sep = "\t", quote = FALSE, row.names = FALSE)
  close(gz3)

  # gene_name is the cross-file join key (it matches phenotype_id in the cis-pairs
  # file and pid in the expression BED); gene_id is a separate Ensembl-style
  # annotation that must NOT be used for joining, matching real tensorQTL output.
  ann <- data.frame(
    gene_id = c("ENSG_GENE1", "ENSG_GENE2", "ENSG_GENE3", "ENSG_GENE4"),
    gene_name = c("GENE1", "GENE2", "GENE3", "GENE4_NOT_SIG"),
    variant_id = c(shared_variant, gene2_variants[1], "V3_2", "V4_1"),
    qval = c(0.01, 0.02, 0.03, 0.5),
    pval_nominal = c(1e-10, 0.2, 1e-6, 0.6),
    slope = c(3, 0.01, 2, 0.001),
    slope_se = c(0.2, 0.3, 0.25, 0.4)
  )
  if (include_num_var) {
    ann$num_var <- c(length(gene1_variants), length(gene2_variants), length(gene3_variants), 2)
  }
  ann_path <- file.path(ds_dir, paste0(dataset, ".cis_qtl_ann.txt.gz"))
  con <- gzfile(ann_path, "w")
  writeLines(c("## metadata line 1", "## metadata line 2"), con)
  utils::write.table(ann, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)

  list(
    dataset_dir = ds_dir, dataset = dataset, donors = donors,
    index_eqtl_file = ann_path,
    cis_pairs_file = file.path(ds_dir, paste0(dataset, ".cis_qtl_pairs.txt.gz")),
    expression_file = file.path(ds_dir, paste0(dataset, ".gene_expression_normalized.bed.gz")),
    genotype_file = file.path(ds_dir, paste0(dataset, ".genotype_matrix.bed.gz")),
    covariate_file = file.path(ds_dir, paste0(dataset, ".covariates_peer.txt")),
    gene1_variants = gene1_variants, gene2_variants = gene2_variants,
    gene3_variants = gene3_variants, shared_variant = shared_variant,
    missing_donor_idx = missing_donor_idx
  )
}

out_paths <- function(dir) {
  list(
    verbose = file.path(dir, "verbose.tsv.gz"),
    credible_sets = file.path(dir, "cs.tsv.gz"),
    report = file.path(dir, "report.tsv.gz")
  )
}

rewrite_gz_tsv <- function(path, transform) {
  dt <- data.table::fread(path, header = TRUE)
  dt <- transform(dt)
  gz <- gzfile(path, "w")
  utils::write.table(dt, gz, sep = "\t", quote = FALSE, row.names = FALSE)
  close(gz)
  invisible(NULL)
}

rewrite_txt_tsv <- function(path, transform) {
  dt <- data.table::fread(path, header = TRUE)
  dt <- transform(dt)
  utils::write.table(dt, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(NULL)
}
