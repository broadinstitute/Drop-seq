# Internal file-reading and validation helpers for run_eqtl_finemapping_files().
#
# All functions in this file are internal (not exported). They read and
# validate the five tensorQTL-derived input files described in
# susie_eqtl_finemapping_package_plan.md.

#' @keywords internal
#' @noRd
construct_dataset_paths <- function(dataset_dir) {
  if (!is.character(dataset_dir) || length(dataset_dir) != 1 ||
      is.na(dataset_dir) || !nzchar(dataset_dir)) {
    stop("'dataset_dir' must be a single non-empty character string.", call. = FALSE)
  }
  # Existence of dataset_dir and its five standard files is checked by
  # validate_explicit_input_paths() on the constructed paths below; this function
  # is pure path construction and does not itself require dataset_dir to exist.
  dataset <- basename(sub("/+$", "", dataset_dir))
  list(
    dataset = dataset,
    index_eqtl_file = file.path(dataset_dir, paste0(dataset, ".cis_qtl_ann.txt.gz")),
    cis_pairs_file = file.path(dataset_dir, paste0(dataset, ".cis_qtl_pairs.txt.gz")),
    expression_file = file.path(dataset_dir, paste0(dataset, ".gene_expression_normalized.bed.gz")),
    genotype_file = file.path(dataset_dir, paste0(dataset, ".genotype_matrix.bed.gz")),
    covariate_file = file.path(dataset_dir, paste0(dataset, ".covariates_peer.txt"))
  )
}

#' @keywords internal
#' @noRd
validate_explicit_input_paths <- function(index_eqtl_file, cis_pairs_file, expression_file,
                                           genotype_file, covariate_file) {
  paths <- c(
    index_eqtl_file = index_eqtl_file, cis_pairs_file = cis_pairs_file,
    expression_file = expression_file, genotype_file = genotype_file,
    covariate_file = covariate_file
  )
  problems <- character(0)
  for (nm in names(paths)) {
    p <- paths[[nm]]
    if (!is.character(p) || length(p) != 1 || is.na(p) || !nzchar(p)) {
      problems <- c(problems, sprintf("%s: must be a single non-empty character string.", nm))
      next
    }
    if (!file.exists(p)) {
      problems <- c(problems, sprintf("%s: file does not exist ('%s').", nm, p))
      next
    }
    if (dir.exists(p)) {
      problems <- c(problems, sprintf("%s: path is a directory, not a file ('%s').", nm, p))
      next
    }
    if (file.access(p, 4) != 0) {
      problems <- c(problems, sprintf("%s: file is not readable ('%s').", nm, p))
    }
  }
  if (length(problems) > 0) {
    stop(paste(c("Invalid input file path(s):", problems), collapse = "\n"), call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
validate_scalar_arguments <- function(qvalue_threshold, L, coverage, n_cores,
                                       verbose_outfile, credible_set_outfile, report_outfile) {
  if (!is.numeric(qvalue_threshold) || length(qvalue_threshold) != 1 ||
      !is.finite(qvalue_threshold) || qvalue_threshold < 0 || qvalue_threshold > 1) {
    stop("'qvalue_threshold' must be a single finite numeric value between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(L) || length(L) != 1 || !is.finite(L) || L <= 0 || L != as.integer(L)) {
    stop("'L' must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(coverage) || length(coverage) != 1 || !is.finite(coverage) ||
      coverage < 0 || coverage > 1) {
    stop("'coverage' must be a single finite numeric value between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1 || !is.finite(n_cores) ||
      n_cores <= 0 || n_cores != as.integer(n_cores)) {
    stop("'n_cores' must be a single positive integer.", call. = FALSE)
  }
  outfiles <- c(verbose = verbose_outfile, credible_sets = credible_set_outfile, report = report_outfile)
  bad <- vapply(outfiles, function(x) !is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x), logical(1))
  if (any(bad)) {
    stop("'verbose_outfile', 'credible_set_outfile', and 'report_outfile' must each be a single non-empty character string.",
         call. = FALSE)
  }
  if (length(unique(outfiles)) != length(outfiles)) {
    stop("'verbose_outfile', 'credible_set_outfile', and 'report_outfile' must be distinct paths.", call. = FALSE)
  }
  for (f in outfiles) {
    parent <- dirname(f)
    if (!dir.exists(parent)) {
      stop(sprintf("Output directory does not exist: '%s'.", parent), call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
derive_dataset_label <- function(index_eqtl_file, cis_pairs_file, expression_file,
                                  genotype_file, covariate_file) {
  suffixes <- c(
    "\\.cis_qtl_ann\\.txt\\.gz$",
    "\\.cis_qtl_pairs\\.txt\\.gz$",
    "\\.gene_expression_normalized\\.bed\\.gz$",
    "\\.genotype_matrix\\.bed\\.gz$",
    "\\.covariates_peer\\.txt$"
  )
  files <- c(index_eqtl_file, cis_pairs_file, expression_file, genotype_file, covariate_file)
  bases <- basename(files)
  matched <- vapply(seq_along(files), function(i) grepl(suffixes[i], bases[i]), logical(1))
  prefixes <- vapply(seq_along(files), function(i) sub(suffixes[i], "", bases[i]), character(1))
  if (all(matched) && length(unique(prefixes)) == 1 && nzchar(prefixes[1])) {
    prefixes[1]
  } else {
    NA_character_
  }
}

#' @keywords internal
#' @noRd
read_bed_donors <- function(path, file_role) {
  header <- data.table::fread(path, header = TRUE, nrows = 0)
  required_prefix <- c("#chr", "start", "end", "pid")
  actual_prefix <- names(header)[seq_len(min(4, ncol(header)))]
  if (!identical(actual_prefix, required_prefix)) {
    stop(sprintf(
      "%s file '%s' must begin with columns %s; found %s.",
      file_role, path, paste(required_prefix, collapse = ", "), paste(actual_prefix, collapse = ", ")
    ), call. = FALSE)
  }
  donors <- names(header)[-seq_len(4)]
  if (length(donors) == 0) {
    stop(sprintf("%s file '%s' has no donor columns.", file_role, path), call. = FALSE)
  }
  dup <- donors[duplicated(donors)]
  if (length(dup) > 0) {
    stop(sprintf("%s file '%s' has duplicated donor column name(s): %s.",
                 file_role, path, paste(unique(dup), collapse = ", ")), call. = FALSE)
  }
  donors
}

#' @keywords internal
#' @noRd
read_covariate_donors <- function(path) {
  header <- data.table::fread(path, header = TRUE, nrows = 0)
  if (ncol(header) < 2) {
    stop(sprintf("Covariate file '%s' must have an ID column plus at least one donor column.", path),
         call. = FALSE)
  }
  donors <- names(header)[-1]
  dup <- donors[duplicated(donors)]
  if (length(dup) > 0) {
    stop(sprintf("Covariate file '%s' has duplicated donor column name(s): %s.",
                 path, paste(unique(dup), collapse = ", ")), call. = FALSE)
  }
  donors
}

#' @keywords internal
#' @noRd
validate_identical_donors <- function(expression_donors, genotype_donors, covariate_donors) {
  problems <- character(0)
  check_one <- function(label, donors) {
    missing_d <- setdiff(expression_donors, donors)
    extra_d <- setdiff(donors, expression_donors)
    msgs <- character(0)
    if (length(missing_d) > 0) {
      msgs <- c(msgs, sprintf("%s is missing donor(s) present in the expression file: %s.",
                               label, paste(missing_d, collapse = ", ")))
    }
    if (length(extra_d) > 0) {
      msgs <- c(msgs, sprintf("%s has extra donor(s) not present in the expression file: %s.",
                               label, paste(extra_d, collapse = ", ")))
    }
    msgs
  }
  problems <- c(problems, check_one("Genotype file", genotype_donors))
  problems <- c(problems, check_one("Covariate file", covariate_donors))
  if (length(problems) > 0) {
    stop(paste(c("Donor sets do not match across input files:", problems), collapse = "\n"), call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
count_leading_comment_lines <- function(path) {
  con <- if (R.utils::isGzipped(path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  n <- 0L
  repeat {
    line <- readLines(con, n = 1L)
    if (length(line) == 0) break
    if (!startsWith(line, "#")) break
    n <- n + 1L
  }
  n
}

#' @keywords internal
#' @noRd
read_index_eqtls <- function(index_eqtl_file, qvalue_threshold) {
  skip_n <- count_leading_comment_lines(index_eqtl_file)
  dt <- data.table::fread(index_eqtl_file, header = TRUE, skip = skip_n)
  required_cols <- c("gene_id", "gene_name", "variant_id", "qval")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop(sprintf("Index eQTL file '%s' is missing required column(s): %s.",
                 index_eqtl_file, paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (any(is.na(dt$gene_id) | !nzchar(as.character(dt$gene_id)))) {
    stop(sprintf("Index eQTL file '%s' has missing 'gene_id' values.", index_eqtl_file), call. = FALSE)
  }
  if (any(is.na(dt$gene_name) | !nzchar(as.character(dt$gene_name)))) {
    stop(sprintf("Index eQTL file '%s' has missing 'gene_name' values.", index_eqtl_file), call. = FALSE)
  }
  if (any(is.na(dt$variant_id) | !nzchar(as.character(dt$variant_id)))) {
    stop(sprintf("Index eQTL file '%s' has missing 'variant_id' values.", index_eqtl_file), call. = FALSE)
  }
  if (!is.numeric(dt$qval)) {
    stop(sprintf("Index eQTL file '%s' column 'qval' must be numeric.", index_eqtl_file), call. = FALSE)
  }
  if (any(is.na(dt$qval))) {
    stop(sprintf("Index eQTL file '%s' has missing 'qval' values.", index_eqtl_file), call. = FALSE)
  }
  retained <- dt[dt$qval <= qvalue_threshold, ]
  if (nrow(retained) == 0) {
    stop(sprintf("No significant eQTLs found in '%s' at qvalue_threshold = %g.",
                 index_eqtl_file, qvalue_threshold), call. = FALSE)
  }
  dup_genes <- retained$gene_id[duplicated(retained$gene_id)]
  if (length(dup_genes) > 0) {
    stop(sprintf("Index eQTL file '%s' has duplicate retained gene_id value(s): %s.",
                 index_eqtl_file, paste(unique(dup_genes), collapse = ", ")), call. = FALSE)
  }
  # gene_name is the cross-file join key (it matches phenotype_id in the cis-pairs
  # file and pid in the expression BED, not the Ensembl-style gene_id column), so
  # a duplicate here would make those joins ambiguous.
  dup_gene_names <- retained$gene_name[duplicated(retained$gene_name)]
  if (length(dup_gene_names) > 0) {
    stop(sprintf("Index eQTL file '%s' has duplicate retained gene_name value(s): %s.",
                 index_eqtl_file, paste(unique(dup_gene_names), collapse = ", ")), call. = FALSE)
  }
  key <- paste(retained$gene_id, retained$variant_id, sep = "")
  if (any(duplicated(key))) {
    stop(sprintf("Index eQTL file '%s' has duplicate gene_id/variant_id row(s).", index_eqtl_file),
         call. = FALSE)
  }
  retained
}

#' @keywords internal
#' @noRd
read_cis_pairs <- function(cis_pairs_file, index_eqtls) {
  dt <- data.table::fread(cis_pairs_file, header = TRUE)
  required_cols <- c("phenotype_id", "variant_id")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop(sprintf("Cis-pairs file '%s' is missing required column(s): %s.",
                 cis_pairs_file, paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  # The cross-file join key is gene_name (it matches phenotype_id here and pid in
  # the expression BED); gene_id (Ensembl) is a preserved annotation, not a key.
  gene_names <- index_eqtls$gene_name
  filtered <- dt[dt$phenotype_id %in% gene_names, ]
  gene_groups <- split(filtered, filtered$phenotype_id, drop = TRUE)

  num_var_warnings <- character(0)
  has_num_var <- "num_var" %in% names(index_eqtls)

  for (i in seq_len(nrow(index_eqtls))) {
    gene_id <- index_eqtls$gene_id[i]
    gene_name <- index_eqtls$gene_name[i]
    index_variant_id <- index_eqtls$variant_id[i]
    gene_rows <- gene_groups[[gene_name]]
    if (is.null(gene_rows) || nrow(gene_rows) == 0) {
      stop(sprintf("Gene '%s' (gene_id '%s') has no rows in cis-pairs file '%s'.",
                   gene_name, gene_id, cis_pairs_file), call. = FALSE)
    }
    dup_variants <- gene_rows$variant_id[duplicated(gene_rows$variant_id)]
    if (length(dup_variants) > 0) {
      stop(sprintf("Gene '%s' has duplicate variant_id value(s) in cis-pairs file '%s': %s.",
                   gene_name, cis_pairs_file, paste(unique(dup_variants), collapse = ", ")), call. = FALSE)
    }
    if (!(index_variant_id %in% gene_rows$variant_id)) {
      stop(sprintf(
        "Index SNP '%s' for gene '%s' is not present among that gene's cis-pair rows in '%s'.",
        index_variant_id, gene_name, cis_pairs_file
      ), call. = FALSE)
    }
    if (has_num_var) {
      expected <- index_eqtls$num_var[i]
      if (!is.na(expected) && expected != nrow(gene_rows)) {
        msg <- sprintf(
          "Gene '%s' (index SNP '%s'): num_var (%s) does not match the number of cis-pair rows (%d).",
          gene_name, index_variant_id, expected, nrow(gene_rows)
        )
        warning(msg, call. = FALSE)
        num_var_warnings[gene_name] <- msg
      }
    }
  }

  list(cis_pairs = filtered, num_var_warnings = num_var_warnings)
}

#' @keywords internal
#' @noRd
resolve_covariates <- function(covariate_file, covariates, donor_order) {
  dt <- data.table::fread(covariate_file, header = TRUE)
  if (ncol(dt) < 2) {
    stop(sprintf("Covariate file '%s' must have an ID column plus at least one donor column.",
                 covariate_file), call. = FALSE)
  }
  id_col <- names(dt)[1]
  row_names <- as.character(dt[[id_col]])
  dup_rows <- row_names[duplicated(row_names)]
  if (length(dup_rows) > 0) {
    stop(sprintf("Covariate file '%s' has duplicated covariate row name(s): %s.",
                 covariate_file, paste(unique(dup_rows), collapse = ", ")), call. = FALSE)
  }
  donor_cols <- names(dt)[-1]
  missing_in_file <- setdiff(donor_order, donor_cols)
  if (length(missing_in_file) > 0) {
    stop(sprintf("Covariate file '%s' is missing donor column(s): %s.",
                 covariate_file, paste(missing_in_file, collapse = ", ")), call. = FALSE)
  }

  if (is.null(covariates)) {
    selected_idx <- seq_len(nrow(dt))
    selected_names <- row_names
  } else {
    if (!is.character(covariates) || length(covariates) == 0 || any(is.na(covariates))) {
      stop("'covariates' must be a non-empty character vector of covariate row names.", call. = FALSE)
    }
    dup_requested <- covariates[duplicated(covariates)]
    if (length(dup_requested) > 0) {
      stop(sprintf("'covariates' contains duplicated name(s): %s.",
                   paste(unique(dup_requested), collapse = ", ")), call. = FALSE)
    }
    missing_requested <- setdiff(covariates, row_names)
    if (length(missing_requested) > 0) {
      stop(sprintf("Requested covariate row name(s) not found in '%s': %s.",
                   covariate_file, paste(missing_requested, collapse = ", ")), call. = FALSE)
    }
    selected_idx <- match(covariates, row_names)
    selected_names <- covariates
  }

  if (length(selected_idx) == 0 && !is.null(covariates)) {
    stop("'covariates' selected zero rows; at least one covariate row is required.", call. = FALSE)
  }

  sub_dt <- dt[selected_idx, donor_order, with = FALSE]
  non_numeric_cols <- names(sub_dt)[!vapply(sub_dt, is.numeric, logical(1))]
  if (length(non_numeric_cols) > 0) {
    stop(sprintf("Covariate file '%s' has non-numeric selected covariate value(s) in donor column(s): %s.",
                 covariate_file, paste(non_numeric_cols, collapse = ", ")), call. = FALSE)
  }
  mat <- as.matrix(sub_dt)
  storage.mode(mat) <- "double"
  rownames(mat) <- selected_names
  if (any(!is.finite(mat))) {
    stop(sprintf("Covariate file '%s' has missing or non-finite selected covariate value(s).",
                 covariate_file), call. = FALSE)
  }
  mat
}

#' @keywords internal
#' @noRd
read_expression_rows <- function(expression_file, pids, donor_order) {
  dt <- data.table::fread(expression_file, header = TRUE)
  required_prefix <- c("#chr", "start", "end", "pid")
  actual_prefix <- names(dt)[seq_len(min(4, ncol(dt)))]
  if (!identical(actual_prefix, required_prefix)) {
    stop(sprintf(
      "Expression file '%s' must begin with columns %s; found %s.",
      expression_file, paste(required_prefix, collapse = ", "), paste(actual_prefix, collapse = ", ")
    ), call. = FALSE)
  }
  missing_donors <- setdiff(donor_order, names(dt))
  if (length(missing_donors) > 0) {
    stop(sprintf("Expression file '%s' is missing donor column(s): %s.",
                 expression_file, paste(missing_donors, collapse = ", ")), call. = FALSE)
  }
  sub <- dt[dt$pid %in% pids, ]
  pid_counts <- table(sub$pid)
  missing_genes <- setdiff(pids, names(pid_counts))
  if (length(missing_genes) > 0) {
    stop(sprintf("Expression file '%s' has no row for phenotype ID(s): %s.",
                 expression_file, paste(missing_genes, collapse = ", ")), call. = FALSE)
  }
  dup_genes <- names(pid_counts)[pid_counts > 1]
  if (length(dup_genes) > 0) {
    stop(sprintf("Expression file '%s' has duplicate row(s) for phenotype ID(s): %s.",
                 expression_file, paste(dup_genes, collapse = ", ")), call. = FALSE)
  }
  donor_sub <- sub[, donor_order, with = FALSE]
  non_numeric <- names(donor_sub)[!vapply(donor_sub, is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(sprintf("Expression file '%s' has non-numeric donor value(s) in column(s): %s.",
                 expression_file, paste(non_numeric, collapse = ", ")), call. = FALSE)
  }
  mat <- as.matrix(donor_sub)
  storage.mode(mat) <- "double"
  rownames(mat) <- sub$pid
  mat <- mat[pids, , drop = FALSE]
  if (any(!is.finite(mat))) {
    stop(sprintf(
      "Expression file '%s' has missing or non-finite expression value(s) for one or more significant genes.",
      expression_file
    ), call. = FALSE)
  }
  mat
}

#' @keywords internal
#' @noRd
read_genotype_rows <- function(genotype_file, variant_ids, donor_order) {
  dt <- data.table::fread(genotype_file, header = TRUE)
  required_prefix <- c("#chr", "start", "end", "pid")
  actual_prefix <- names(dt)[seq_len(min(4, ncol(dt)))]
  if (!identical(actual_prefix, required_prefix)) {
    stop(sprintf(
      "Genotype file '%s' must begin with columns %s; found %s.",
      genotype_file, paste(required_prefix, collapse = ", "), paste(actual_prefix, collapse = ", ")
    ), call. = FALSE)
  }
  missing_donors <- setdiff(donor_order, names(dt))
  if (length(missing_donors) > 0) {
    stop(sprintf("Genotype file '%s' is missing donor column(s): %s.",
                 genotype_file, paste(missing_donors, collapse = ", ")), call. = FALSE)
  }
  sub <- dt[dt$pid %in% variant_ids, ]
  pid_counts <- table(sub$pid)
  missing_variants <- setdiff(variant_ids, names(pid_counts))
  if (length(missing_variants) > 0) {
    stop(sprintf("Genotype file '%s' has no row for required cis variant(s): %s.",
                 genotype_file, paste(missing_variants, collapse = ", ")), call. = FALSE)
  }
  dup_variants <- names(pid_counts)[pid_counts > 1]
  if (length(dup_variants) > 0) {
    stop(sprintf("Genotype file '%s' has duplicate row(s) for variant(s): %s.",
                 genotype_file, paste(dup_variants, collapse = ", ")), call. = FALSE)
  }
  donor_sub <- sub[, donor_order, with = FALSE]
  non_numeric <- names(donor_sub)[!vapply(donor_sub, is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(sprintf("Genotype file '%s' has non-numeric dosage value(s) in column(s): %s.",
                 genotype_file, paste(non_numeric, collapse = ", ")), call. = FALSE)
  }
  mat <- as.matrix(donor_sub)
  storage.mode(mat) <- "double"
  rownames(mat) <- sub$pid
  mat <- mat[variant_ids, , drop = FALSE]
  bad <- is.infinite(mat) | is.nan(mat)
  if (any(bad)) {
    stop(sprintf(
      "Genotype file '%s' has non-finite (Inf/NaN) dosage value(s); missing dosage must be encoded as NA.",
      genotype_file
    ), call. = FALSE)
  }
  mat
}
