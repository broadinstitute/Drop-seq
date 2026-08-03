#' Run individual-level SuSiE fine-mapping from explicit tensorQTL-derived files
#'
#' For every gene whose tensorQTL index association has \code{qval <= qvalue_threshold}
#' in \code{index_eqtl_file}, this function identifies all cis variants tested for that
#' gene (from \code{cis_pairs_file}), loads the gene's normalized expression vector (from
#' \code{expression_file}) and donor-level dosages for all cis variants (from
#' \code{genotype_file}), optionally residualizes expression and genotypes against the
#' covariate matrix (from \code{covariate_file}) using QR decomposition, and fits
#' \code{susieR::susie()} using individual-level data. Credible sets are extracted with
#' \code{susieR::susie_get_cs()} for every successful, converged fit.
#'
#' \code{susieR::susie_get_cs()} computes credible-set purity faster for genes with many
#' cis variants when the \pkg{Rfast} package is installed; if it is absent,
#' \code{susieR} prints an informational hint and falls back to a slower calculation.
#' This has no effect on the returned results, only on speed. \pkg{Rfast} is an optional
#' (\code{Suggests}) dependency of this package and is not required for correctness.
#'
#' The output \code{pip} column is the marginal posterior inclusion probability
#' (\code{fit$pip} from \code{susieR::susie()}), which unions evidence across every
#' active single-effect component for the gene; it need not sum to 1 (or to
#' \code{coverage}) across the members of one credible set. The output
#' \code{credible_set_alpha} column is instead a single-effect component's own
#' \code{alpha} distribution (\code{fit$alpha[l, ]}), which by construction sums to
#' exactly 1 across all cis variants tested for that component. When a gene has exactly
#' one credible set, its \code{credible_set_alpha} is reported for every cis variant in
#' the verbose output (not only the credible set's members), so
#' \code{sum(credible_set_alpha)} across the whole gene equals 1. When a gene has more
#' than one credible set, each set's underlying variant is a mixture of several
#' components' mass with no single attributable value, so \code{credible_set_alpha} is
#' \code{NA} outside each set's own members, and \code{sum(credible_set_alpha)} across
#' just one set's members instead approximates that set's \code{coverage}.
#'
#' Donors are matched by exact donor ID across the expression BED, genotype BED, and
#' covariate file; the three donor sets must be identical (a mismatch is a fatal error).
#' Donors with a missing genotype dosage for any cis variant of a given gene are dropped
#' for that gene's analysis only. No external LD reference is used and no SuSiE fit
#' objects are saved to disk.
#'
#' Genes are joined across \code{index_eqtl_file}, \code{cis_pairs_file}, and
#' \code{expression_file} by \code{gene_name}, not by \code{gene_id}. In tensorQTL output,
#' \code{phenotype_id} in the cis-pairs file and \code{pid} in the expression BED carry the
#' original phenotype identifier used throughout the tensorQTL run (typically a gene symbol,
#' falling back to an Ensembl ID when no symbol is available); this is exactly what ends up
#' in the annotated index-eQTL file's \code{gene_name} column. \code{gene_id} in
#' \code{index_eqtl_file} is a separately annotated Ensembl identifier and is carried
#' through to the output tables but is not used as a join key.
#'
#' A gene-level SuSiE fitting error, non-convergence, or credible-set extraction error is
#' caught, reported via \code{\link{warning}} (identifying \code{gene_name} and the index
#' variant) and recorded in the processing report; it does not stop the run for other
#' genes. In that case the gene's tensorQTL index SNP is emitted as a fallback row in the
#' slim credible-set output with \code{in_credible_set = FALSE}, and all SuSiE-derived
#' fields are \code{NA} in the verbose output. A successfully fit gene that reports no
#' credible set is not a warning condition.
#'
#' @param index_eqtl_file Path to the tensorQTL annotated index-eQTL file
#'   (\code{*.cis_qtl_ann.txt.gz}). May begin with metadata lines starting with \code{#}.
#'   Must contain columns \code{gene_id}, \code{gene_name}, \code{variant_id}, \code{qval}.
#' @param cis_pairs_file Path to the tensorQTL all-cis-pairs file
#'   (\code{*.cis_qtl_pairs.txt.gz}). Must contain columns \code{phenotype_id},
#'   \code{variant_id}.
#' @param expression_file Path to the normalized expression BED
#'   (\code{*.gene_expression_normalized.bed.gz}) with columns \code{#chr}, \code{start},
#'   \code{end}, \code{pid}, followed by one column per donor.
#' @param genotype_file Path to the genotype dosage BED
#'   (\code{*.genotype_matrix.bed.gz}) with columns \code{#chr}, \code{start}, \code{end},
#'   \code{pid}, followed by one column per donor. Missing dosage values are permitted.
#' @param covariate_file Path to the tensorQTL covariate file (\code{*.covariates_peer.txt})
#'   with one covariate per row, an identifier column, and one column per donor.
#' @param covariates Optional character vector of exact covariate row names to use from
#'   \code{covariate_file}. \code{NULL} (the default) uses every covariate row in the file.
#'   No prefix matching or one-hot expansion is performed; requested names must match row
#'   names exactly.
#' @param qvalue_threshold Maximum tensorQTL \code{qval} for an index eQTL to be selected
#'   for fine-mapping. Default \code{0.05}.
#' @param L Maximum number of single effects passed to \code{susieR::susie()}. Default
#'   \code{10}.
#' @param coverage Requested credible-set coverage passed to \code{susieR::susie()} and
#'   \code{susieR::susie_get_cs()}. Default \code{0.95}.
#' @param n_cores Number of CPU cores to use for the per-gene SuSiE fitting step, via
#'   \code{pbmcapply::pbmclapply()} (fork-based; shows a progress bar). Default \code{1}
#'   (sequential; no forking). Fork-based parallelism is only available on Unix-alikes
#'   (Linux, macOS); on Windows \code{pbmclapply()} runs sequentially regardless of
#'   \code{n_cores}. \code{susieR::susie()} fitting is deterministic, so the three
#'   returned tables are numerically identical regardless of \code{n_cores}; warnings for
#'   gene-level failures are always emitted in the same gene order regardless of how many
#'   cores are used or the order in which workers finish.
#' @param verbose_outfile Path to write the verbose, one-row-per-analyzed-gene-SNP-pair
#'   result table (recommended suffix \code{.tsv.gz}).
#' @param credible_set_outfile Path to write the slim, one-row-per-credible-set-gene-SNP-pair
#'   result table (recommended suffix \code{.tsv.gz}).
#' @param report_outfile Path to write the one-row-per-gene processing report
#'   (recommended suffix \code{.tsv.gz} or \code{.tsv}).
#' @param dataset Optional dataset label recorded in the \code{dataset} output column.
#'   When \code{NULL} (the default), a label is derived from the five input file names
#'   when they share an unambiguous common prefix; otherwise \code{NA_character_} is used.
#'
#' @return Invisibly, a \code{list} with three base \code{data.frame} elements:
#'   \code{verbose}, \code{credible_sets}, and \code{report}. The same three tables are
#'   written to \code{verbose_outfile}, \code{credible_set_outfile}, and
#'   \code{report_outfile} respectively, as tab-separated files (gzip-compressed when the
#'   supplied path ends in \code{.gz}) via \code{data.table::fwrite()}.
#'
#' @section Errors: This function calls \code{\link{stop}} for missing or unreadable
#'   input files, invalid arguments, no eQTLs passing \code{qvalue_threshold}, a donor-set
#'   mismatch across the expression/genotype/covariate files, missing or duplicated
#'   required rows, an index SNP absent from its gene's cis-pair rows, missing or
#'   non-finite expression or covariate values, monomorphic cis variants, and other
#'   structural input inconsistencies described in the package vignette-level
#'   documentation.
#'
#' @section Warnings: This function calls \code{\link{warning}} (and continues
#'   processing other genes) for a gene-level SuSiE fitting error, non-convergence, a
#'   credible-set extraction error, or a mismatch between the annotated \code{num_var}
#'   and the observed number of cis-pair rows for a gene.
#'
#' @examples
#' \dontrun{
#' run_eqtl_finemapping_files(
#'   index_eqtl_file = "SPN_D1__CaH.cis_qtl_ann.txt.gz",
#'   cis_pairs_file = "SPN_D1__CaH.cis_qtl_pairs.txt.gz",
#'   expression_file = "SPN_D1__CaH.gene_expression_normalized.bed.gz",
#'   genotype_file = "SPN_D1__CaH.genotype_matrix.bed.gz",
#'   covariate_file = "SPN_D1__CaH.covariates_peer.txt",
#'   verbose_outfile = "SPN_D1__CaH.susie_verbose.tsv.gz",
#'   credible_set_outfile = "SPN_D1__CaH.susie_credible_sets.tsv.gz",
#'   report_outfile = "SPN_D1__CaH.susie_report.tsv.gz"
#' )
#' }
#'
#' @export
run_eqtl_finemapping_files <- function(index_eqtl_file, cis_pairs_file, expression_file,
                                        genotype_file, covariate_file, covariates = NULL,
                                        qvalue_threshold = 0.05, L = 10, coverage = 0.95,
                                        n_cores = 1,
                                        verbose_outfile, credible_set_outfile, report_outfile,
                                        dataset = NULL) {
  validate_explicit_input_paths(index_eqtl_file, cis_pairs_file, expression_file,
                                 genotype_file, covariate_file)
  validate_scalar_arguments(qvalue_threshold, L, coverage, n_cores, verbose_outfile,
                             credible_set_outfile, report_outfile)

  if (is.null(dataset)) {
    dataset <- derive_dataset_label(index_eqtl_file, cis_pairs_file, expression_file,
                                     genotype_file, covariate_file)
  } else if (!is.character(dataset) || length(dataset) != 1) {
    stop("'dataset' must be a single character string or NULL.", call. = FALSE)
  }

  expression_donors <- read_bed_donors(expression_file, "Expression")
  genotype_donors <- read_bed_donors(genotype_file, "Genotype")
  covariate_donors <- read_covariate_donors(covariate_file)
  validate_identical_donors(expression_donors, genotype_donors, covariate_donors)
  donor_order <- expression_donors

  index_eqtls <- read_index_eqtls(index_eqtl_file, qvalue_threshold)
  cis_pairs_result <- read_cis_pairs(cis_pairs_file, index_eqtls)
  cis_pairs <- cis_pairs_result$cis_pairs
  num_var_warnings <- cis_pairs_result$num_var_warnings

  covariate_matrix <- resolve_covariates(covariate_file, covariates, donor_order)

  # gene_name (not the Ensembl-style gene_id) is the cross-file join key: it is
  # what matches phenotype_id in the cis-pairs file and pid in the expression BED.
  gene_names <- index_eqtls$gene_name
  expression_matrix <- read_expression_rows(expression_file, gene_names, donor_order)

  cis_pairs_by_gene <- split(cis_pairs, cis_pairs$phenotype_id, drop = TRUE)
  all_variant_ids <- unique(cis_pairs$variant_id)
  genotype_matrix <- read_genotype_rows(genotype_file, all_variant_ids, donor_order)

  n_genes <- nrow(index_eqtls)
  gene_rows <- lapply(seq_len(n_genes), function(i) index_eqtls[i, , drop = FALSE])

  run_one_gene <- function(gene_row) {
    gene_cis_pairs <- cis_pairs_by_gene[[gene_row$gene_name]]
    tryCatch(
      process_one_gene(gene_row, gene_cis_pairs, expression_matrix, genotype_matrix,
                        covariate_matrix, donor_order, dataset, L, coverage),
      error = function(e) e
    )
  }

  # A gene-level SuSiE/credible-set failure is caught inside process_one_gene() and
  # never raises an R condition here; only a genuinely fatal, run-aborting error (e.g.
  # a monomorphic variant) reaches this tryCatch. Collecting results into a plain list
  # first - rather than calling stop()/warning() from inside the (possibly forked)
  # worker - keeps fatal-error detection and warning order deterministic and identical
  # whether n_cores is 1 or many, since both are resolved below in fixed gene order.
  if (n_cores > 1) {
    gene_results <- pbmcapply::pbmclapply(gene_rows, run_one_gene, mc.cores = n_cores)
  } else {
    gene_results <- lapply(gene_rows, run_one_gene)
  }

  for (i in seq_len(n_genes)) {
    if (inherits(gene_results[[i]], "error")) {
      stop(conditionMessage(gene_results[[i]]), call. = FALSE)
    }
  }

  combined <- combine_gene_results(gene_results, index_eqtls$gene_name, num_var_warnings)
  verbose_output <- combined$verbose
  credible_set_output <- combined$credible_sets
  report_output <- combined$report

  ord <- order(verbose_output$gene_id, verbose_output$in_credible_set, verbose_output$pip,
               verbose_output$variant_id, method = "radix", decreasing = c(FALSE, TRUE, TRUE, FALSE))
  verbose_output <- verbose_output[ord, ]
  rownames(verbose_output) <- NULL

  validate_output_keys(verbose_output)
  validate_output_keys(credible_set_output)

  write_finemapping_outputs(verbose_output, credible_set_output, report_output,
                             verbose_outfile, credible_set_outfile, report_outfile)

  invisible(list(verbose = verbose_output, credible_sets = credible_set_output, report = report_output))
}

#' Run individual-level SuSiE fine-mapping from a standard tensorQTL dataset directory
#'
#' Directory-layout adapter for \code{\link{run_eqtl_finemapping_files}}. Given a dataset
#' directory named \code{<cell_type>__<region>}, this function derives \code{<dataset>}
#' from \code{basename(dataset_dir)}, constructs the five standard tensorQTL-derived file
#' paths expected directly inside \code{dataset_dir}:
#' \itemize{
#'   \item \code{<dataset>.cis_qtl_ann.txt.gz}
#'   \item \code{<dataset>.cis_qtl_pairs.txt.gz}
#'   \item \code{<dataset>.gene_expression_normalized.bed.gz}
#'   \item \code{<dataset>.genotype_matrix.bed.gz}
#'   \item \code{<dataset>.covariates_peer.txt}
#' }
#' validates that all five files exist and are readable, and delegates all analysis to
#' \code{\link{run_eqtl_finemapping_files}} with \code{dataset} set to the directory
#' basename. This function performs no analysis of its own.
#'
#' @param dataset_dir Path to one dataset directory named \code{<cell_type>__<region>}
#'   containing the five standard tensorQTL-derived files described above.
#' @param covariates Optional character vector of exact covariate row names to use.
#'   \code{NULL} (the default) uses every covariate row in the covariate file. See
#'   \code{\link{run_eqtl_finemapping_files}}.
#' @param qvalue_threshold Maximum tensorQTL \code{qval} for an index eQTL to be selected
#'   for fine-mapping. Default \code{0.05}.
#' @param L Maximum number of single effects passed to \code{susieR::susie()}. Default
#'   \code{10}.
#' @param coverage Requested credible-set coverage. Default \code{0.95}.
#' @param n_cores Number of CPU cores to use for the per-gene SuSiE fitting step. Default
#'   \code{1} (sequential). See \code{\link{run_eqtl_finemapping_files}} for details.
#' @param verbose_outfile Path to write the verbose, one-row-per-analyzed-gene-SNP-pair
#'   result table (recommended suffix \code{.tsv.gz}).
#' @param credible_set_outfile Path to write the slim, one-row-per-credible-set-gene-SNP-pair
#'   result table (recommended suffix \code{.tsv.gz}).
#' @param report_outfile Path to write the one-row-per-gene processing report
#'   (recommended suffix \code{.tsv.gz} or \code{.tsv}).
#'
#' @return Invisibly, the same \code{list(verbose, credible_sets, report)} of base
#'   \code{data.frame}s returned by \code{\link{run_eqtl_finemapping_files}}. The same
#'   three tables are written to the supplied output paths.
#'
#' @seealso \code{\link{run_eqtl_finemapping_files}} for the underlying analysis, its
#'   error conditions, and its warning conditions.
#'
#' @examples
#' \dontrun{
#' run_eqtl_finemapping(
#'   dataset_dir = "SPN_D1__CaH",
#'   verbose_outfile = "SPN_D1__CaH.susie_verbose.tsv.gz",
#'   credible_set_outfile = "SPN_D1__CaH.susie_credible_sets.tsv.gz",
#'   report_outfile = "SPN_D1__CaH.susie_report.tsv.gz"
#' )
#' }
#'
#' @export
run_eqtl_finemapping <- function(dataset_dir, covariates = NULL, qvalue_threshold = 0.05,
                                  L = 10, coverage = 0.95, n_cores = 1,
                                  verbose_outfile, credible_set_outfile, report_outfile) {
  paths <- construct_dataset_paths(dataset_dir)
  validate_explicit_input_paths(paths$index_eqtl_file, paths$cis_pairs_file,
                                 paths$expression_file, paths$genotype_file, paths$covariate_file)
  run_eqtl_finemapping_files(
    index_eqtl_file = paths$index_eqtl_file,
    cis_pairs_file = paths$cis_pairs_file,
    expression_file = paths$expression_file,
    genotype_file = paths$genotype_file,
    covariate_file = paths$covariate_file,
    covariates = covariates,
    qvalue_threshold = qvalue_threshold,
    L = L,
    coverage = coverage,
    n_cores = n_cores,
    verbose_outfile = verbose_outfile,
    credible_set_outfile = credible_set_outfile,
    report_outfile = report_outfile,
    dataset = paths$dataset
  )
}
