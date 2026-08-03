# Project Plan: Run SuSiE Fine-Mapping for Significant eQTLs from Explicit Files or a Dataset Directory

## Objective

Implement a small R package that runs individual-level SuSiE fine-mapping for every significant eQTL gene in one tensorQTL-derived dataset. Inputs may be supplied as five explicit files or through a standard `<cell_type>__<region>` directory adapter.

Each core-function call processes one matched set of index-eQTL, cis-pair, expression, genotype, and covariate files. For each gene whose tensorQTL index association has `qval <= 0.05`, the package will:

1. identify all cis variants tested for that gene;
2. load the gene's normalized expression vector;
3. load donor-level dosage values for all cis variants;
4. optionally adjust expression and genotypes for the same covariate matrix used by tensorQTL;
5. run `susieR::susie()` using individual-level data;
6. emit a verbose SNP-gene result table;
7. emit a slim credible-set table;
8. emit a one-row-per-gene processing summary.

The package will not use an external LD reference and will not save SuSiE fit objects.

## Exported Functions

The package should expose two main functions. The first is the core implementation and receives every input file explicitly. The second is a directory-layout adapter that discovers the standard tensorQTL-derived files and delegates all analysis to the core function.

### Core file-based function

```r
run_eqtl_finemapping_files(
    index_eqtl_file,
    cis_pairs_file,
    expression_file,
    genotype_file,
    covariate_file,
    covariates = NULL,
    qvalue_threshold = 0.05,
    L = 10,
    coverage = 0.95,
    verbose_outfile,
    credible_set_outfile,
    report_outfile,
    dataset = NULL
)
```

This function performs all reading, validation, residualization, SuSiE fitting, output construction, writing, and invisible return. It must not infer related file paths from any one input path.

### Directory adapter

```r
run_eqtl_finemapping(
    dataset_dir,
    covariates = NULL,
    qvalue_threshold = 0.05,
    L = 10,
    coverage = 0.95,
    verbose_outfile,
    credible_set_outfile,
    report_outfile
)
```

The adapter derives `<dataset>` from `basename(dataset_dir)`, constructs the five standard file paths, validates that they exist, and calls `run_eqtl_finemapping_files()` with `dataset = <dataset>`. It must not duplicate the core analysis logic.

### Shared analysis arguments

- `dataset_dir`
  - Adapter-only path to one dataset directory named `<cell_type>__<region>`.
  - The directory is processed independently of all other datasets.

- `index_eqtl_file`, `cis_pairs_file`, `expression_file`, `genotype_file`, `covariate_file`
  - Core-function paths to the five required inputs.
  - Each path is supplied independently and explicitly.

- `dataset`
  - Optional dataset label used in output columns by the core function.
  - When `NULL`, derive a conservative label from the common filename prefix when this can be done unambiguously; otherwise use `NA_character_`.
  - The directory adapter always supplies the directory basename explicitly.

- `covariates`
  - Optional character vector containing exact covariate row names from the tensorQTL covariate file.
  - `NULL` means use every covariate row in the file.
  - No prefix matching or one-hot expansion is performed.
  - A supplied vector must match exact row names.

- `qvalue_threshold`
  - Maximum `qval` for selecting index eQTLs.
  - Default: `0.05`.

- `L`
  - Maximum number of single effects passed to `susieR::susie()`.
  - Default: `10`.

- `coverage`
  - Requested credible-set coverage.
  - Default: `0.95`.

- `verbose_outfile`
  - Required path for the verbose SNP-gene result table.
  - Recommended suffix: `.tsv.gz`.

- `credible_set_outfile`
  - Required path for the slim credible-set result table.
  - Recommended suffix: `.tsv.gz`.

- `report_outfile`
  - Required path for the one-row-per-gene processing report.
  - Recommended suffix: `.tsv.gz` or `.tsv`.

## Return Value

Write all three output files and invisibly return:

```r
invisible(list(
    verbose = verbose_output,
    credible_sets = credible_set_output,
    report = report_output
))
```

Each returned table must be a base `data.frame`. Both exported functions return the same structure because the adapter returns the result of the core function invisibly.

## Standard Dataset Directory Layout Used by the Adapter

Let `<dataset>` be the basename of `dataset_dir`, for example:

```text
SPN_D1__CaH
```

Require these files directly inside `dataset_dir`:

```text
<dataset>.cis_qtl_ann.txt.gz
<dataset>.cis_qtl_pairs.txt.gz
<dataset>.gene_expression_normalized.bed.gz
<dataset>.genotype_matrix.bed.gz
<dataset>.covariates_peer.txt
```

Example:

```text
SPN_D1__CaH/
├── SPN_D1__CaH.cis_qtl_ann.txt.gz
├── SPN_D1__CaH.cis_qtl_pairs.txt.gz
├── SPN_D1__CaH.gene_expression_normalized.bed.gz
├── SPN_D1__CaH.genotype_matrix.bed.gz
└── SPN_D1__CaH.covariates_peer.txt
```

Stop before reading large data matrices if any required file is absent.

## Input File Expectations

### Annotated index-eQTL file

File suffix:

```text
cis_qtl_ann.txt.gz
```

The file may begin with comment lines starting with `#`.

Require at least:

```text
gene_id
gene_name
variant_id
qval
```

Useful fields to preserve when present include:

```text
num_var
pval_nominal
slope
slope_se
pval_perm
pval_beta
chr
variant_pos
ref
alt
rsid
gene_chr
gene_start
gene_end
strand
```

Each retained row defines one gene-level fine-mapping analysis:

```r
qval <= qvalue_threshold
```

Use `gene_id` as the canonical gene identifier and `variant_id` as the original tensorQTL index SNP.

Require exactly one retained index row per `gene_id`.

### All cis-pair file

File suffix:

```text
cis_qtl_pairs.txt.gz
```

Require at least:

```text
phenotype_id
variant_id
```

Useful fields to preserve when present include:

```text
start_distance
af
ma_samples
ma_count
pval_nominal
slope
slope_se
```

For each retained `gene_id`, select all rows where:

```r
phenotype_id == gene_id
```

The selected `variant_id` values define the full cis-variant set supplied to SuSiE for that gene.

Require that the index SNP from `cis_qtl_ann.txt.gz` occurs among the gene's rows in `cis_qtl_pairs.txt.gz`.

### Expression BED

File suffix:

```text
gene_expression_normalized.bed.gz
```

Expect four BED annotation columns followed by donor columns:

```text
#chr
start
end
pid
<donor columns...>
```

Use `pid` to match `gene_id`.

Expression values are already normalized by the eQTL workflow. Preserve them as supplied before optional covariate residualization.

### Genotype dosage BED

File suffix:

```text
genotype_matrix.bed.gz
```

Expect four BED annotation columns followed by donor columns:

```text
#chr
start
end
pid
<donor columns...>
```

Use `pid` to match `variant_id`.

Genotypes are dosage values. Missing dosage values are permitted and handled per gene.

### Covariate file

File suffix:

```text
covariates_peer.txt
```

Expected orientation:

```text
ID    donor1    donor2    donor3 ...
```

Each row is one tensorQTL covariate and donors are columns.

The selected rows constitute the exact covariate design matrix used for adjustment.

## Processing Workflow

### 1. Resolve and validate paths

For `run_eqtl_finemapping_files()`, validate the five explicitly supplied input paths and stop if any are missing or not regular readable files. Do not infer one input path from another.

For `run_eqtl_finemapping()`, derive `<dataset>` from `basename(dataset_dir)`, construct the five standard paths, validate them, and immediately delegate to `run_eqtl_finemapping_files()`.

Validate scalar arguments:

- `qvalue_threshold` must be numeric, finite, and between 0 and 1;
- `L` must be a positive integer;
- `coverage` must be numeric, finite, and between 0 and 1;
- output paths must be distinct and non-empty.

### 2. Read headers and assert exact donor consistency

Read donor IDs from:

- the expression BED header;
- the genotype BED header;
- the covariate-file header.

Require the donor sets to be identical across all three files.

A mismatch is a fatal error. The error message must report, separately for each file:

- donors missing relative to the expression file;
- extra donors relative to the expression file.

After validation, reorder genotype and covariate columns to the expression donor order.

Do not use a donor intersection.

### 3. Read and select significant index eQTLs

Read `cis_qtl_ann.txt.gz`, ignoring leading `#` metadata lines.

Require the needed columns and validate:

- no missing `gene_id`;
- no missing `variant_id`;
- no missing or non-numeric `qval`;
- no duplicate retained `gene_id` rows;
- no duplicate retained gene-index-variant rows.

Retain rows satisfying:

```r
qval <= qvalue_threshold
```

Stop if no rows remain.

Preserve the annotated tensorQTL statistics needed for output.

### 4. Read and validate cis-pair membership

Read `cis_qtl_pairs.txt.gz`.

Retain rows for the significant genes only.

For every significant gene:

- require at least one cis-pair row;
- require unique `variant_id` values within the gene;
- require the index SNP to be present among the cis variants;
- optionally compare the number of selected cis variants with `num_var` when `num_var` is available;
- treat a `num_var` mismatch as a fatal consistency error.

The same variant may legitimately occur for multiple genes.

The expected unique key in this table is:

```text
phenotype_id + variant_id
```

### 5. Read and validate covariates

Read the complete covariate matrix before expensive gene-level processing.

If `covariates = NULL`:

- use every row except the identifier/header row.

If `covariates` is supplied:

- require a non-empty character vector;
- reject duplicated requested names;
- require exact row-name matches;
- stop and list any missing requested names.

After donor reordering, require:

- all selected covariate values to be numeric;
- no missing or non-finite selected covariate values;
- no duplicated covariate row names;
- at least one selected covariate row when a non-`NULL` selection is requested.

Linear dependence among selected covariates, including dependence introduced by adding an intercept to a complete one-hot encoding, is permitted. The package uses the column space of the design matrix for residualization and does not interpret individual covariate coefficients.

Any missing expression or covariate value is a fatal error.

### 6. Read required expression rows

Read only significant-gene rows from the normalized-expression BED where practical.

For every significant gene:

- require exactly one matching `pid` row;
- require all donor values to be numeric and finite;
- stop if any expression value is missing;
- stop on duplicate expression rows.

Expression is the response vector `y` for that gene.

### 7. Read required genotype rows

Build the union of all required cis `variant_id` values across significant genes.

Read only those genotype rows where practical.

For every required variant:

- require exactly one matching genotype BED row;
- stop on missing variants;
- stop on duplicate genotype rows;
- require non-missing values to be numeric and finite.

The same loaded genotype row may be reused for multiple genes.

### 8. Build one analysis data set per gene

For each significant gene:

1. retrieve its expression vector;
2. retrieve all cis-variant dosage rows in the order recorded in `cis_qtl_pairs.txt.gz`;
3. transpose dosage values to construct donor-by-variant matrix `X`;
4. identify donors with a missing dosage in any cis variant;
5. remove those donors from `X`, `y`, and the selected covariate matrix for this gene only.

Record:

```text
n_donors_total
n_donors_used
n_donors_dropped_missing_genotype
```

Do not remove or impute individual SNP values independently. A donor with any missing dosage for the gene is excluded from that entire gene-level fit.

After donor removal, require:

- at least two retained donors;
- no remaining missing genotype values;
- no monomorphic or zero-variance cis variants;
- no missing or non-finite expression values;
- no missing or non-finite covariate values;
- sufficient retained donors for covariate residualization and SuSiE fitting.

A monomorphic or zero-variance cis variant is a fatal run-level error because it indicates an unexpected input inconsistency.

### 9. Remove covariate effects

Follow the individual-level SuSiE vignette approach by regressing the same covariate matrix out of both `X` and `y`.

Let `Z` be the donor-by-covariate matrix for the retained donors. Add an intercept column explicitly.

Use QR decomposition rather than explicitly solving normal equations:

```r
qr_z <- qr(Z, LAPACK = FALSE)
y_resid <- qr.resid(qr_z, y)
X_resid <- qr.resid(qr_z, X)
```

`qr.resid()` accepts a matrix response, so the genotype columns can be residualized together. Redundant one-hot encoded covariate columns are acceptable because the residual projection remains well-defined even when individual regression coefficients are not unique.

Require the residualized response and each residualized genotype column to have non-zero finite variance.

The outputs are interpreted on the residualized scale.

### 10. Run SuSiE

Fit one individual-level model per gene:

```r
fit <- susieR::susie(
    X_resid,
    y_resid,
    L = L,
    coverage = coverage,
    verbose = FALSE
)
```

Use explicit namespace qualification for all third-party functions.

Do not use `susie_rss()` and do not compute or import an external LD matrix.

After every successful SuSiE fit, call `susieR::susie_get_cs()` explicitly using the fitted model, residualized genotype matrix, and caller-specified `coverage`:

```r
cs <- susieR::susie_get_cs(
    fit,
    X = X_resid,
    coverage = coverage
)
```

Use the returned credible sets to construct the slim credible-set output, assign SNP-level credible-set membership in the verbose output, and record credible-set counts and purity summaries in the per-gene report. Credible-set extraction is a required step, not an optional fallback.

Record:

- `fit$converged` when available;
- SNP-level PIP from `fit$pip`;
- posterior mean effect from the SuSiE posterior summaries;
- posterior SD when available or derivable from the fit;
- credible-set membership;
- credible-set coverage and purity summaries where available.

### 11. Handle SuSiE errors and non-convergence

A SuSiE error or non-converged fit must not stop the directory-level run.

Catch fitting and credible-set-extraction errors per gene.

For a failed or non-converged gene:

- issue an R warning identifying `gene_id` and `index_variant_id`;
- record the complete warning/error text in the report;
- emit every cis gene-SNP pair in the verbose output;
- set all SuSiE-derived SNP fields to `NA`;
- emit the original tensorQTL index SNP in the slim output;
- set `in_credible_set = FALSE` for the fallback index-SNP row;
- mark the gene status clearly.

This behavior is intentionally conservative because the frequency and causes of failures must be measurable before choosing a stricter policy.

### 12. Construct the verbose SNP-gene output

Emit one row for every cis gene-SNP pair analyzed.

The same variant may appear for multiple genes. Within one gene, each variant must appear once.

Proposed columns:

```text
dataset
gene_id
gene_name
index_variant_id
variant_id
is_index_variant
qval
index_pval_nominal
index_slope
index_slope_se
cis_start_distance
cis_af
cis_ma_samples
cis_ma_count
cis_pval_nominal
cis_slope
cis_slope_se
pip
posterior_mean
posterior_sd
in_credible_set
credible_set_id
credible_set_requested_coverage
credible_set_coverage
credible_set_min_abs_corr
credible_set_mean_abs_corr
credible_set_median_abs_corr
susie_converged
analysis_status
n_cis_variants
n_donors_total
n_donors_used
n_donors_dropped_missing_genotype
n_covariate_columns
L
coverage
```

Notes:

- `index_variant_id` is repeated for every SNP row for that gene.
- `is_index_variant` is logical.
- `in_credible_set` is logical when the fit succeeds and `NA` when credible-set status cannot be evaluated because fitting failed.
- `credible_set_id` may be `NA` for SNPs outside credible sets.
- If a SNP is unexpectedly assigned to multiple credible sets for the same gene, treat this as a fatal internal consistency error rather than duplicating or collapsing the row.
- Preserve one row per `gene_id + variant_id`.

Sort by:

```text
gene_id
in_credible_set descending
pip descending
variant_id
```

### 13. Construct the slim credible-set output

The slim table must also use one row per gene-SNP pair.

For a successfully fitted gene with one or more credible sets:

- emit every unique SNP that belongs to a reported credible set;
- do not automatically emit the tensorQTL index SNP if it is outside all credible sets.

For a gene with no credible set:

- emit the tensorQTL index SNP as one fallback row;
- set `in_credible_set = FALSE`;
- retain the index SNP's PIP and posterior summaries when SuSiE otherwise succeeded;
- set `analysis_status = "no_credible_set"`.

For a failed or non-converged gene:

- emit the tensorQTL index SNP as one fallback row;
- set `in_credible_set = FALSE`;
- set SuSiE-derived fields to `NA`;
- record failure status.

Proposed columns:

```text
dataset
gene_id
gene_name
index_variant_id
variant_id
is_index_variant
in_credible_set
credible_set_id
pip
posterior_mean
posterior_sd
credible_set_requested_coverage
credible_set_coverage
credible_set_min_abs_corr
credible_set_mean_abs_corr
credible_set_median_abs_corr
qval
susie_converged
analysis_status
n_cis_variants
n_donors_total
n_donors_used
n_donors_dropped_missing_genotype
```

The table key is:

```text
gene_id + variant_id
```

The same `variant_id` may occur in multiple rows when associated with different genes.

### 14. Construct the processing report

Write one row per significant gene.

Required columns:

```text
dataset
gene_id
gene_name
index_variant_id
qval
status
n_cis_variants
n_donors_total
n_donors_used
n_donors_dropped_missing_genotype
n_covariate_columns
susie_converged
n_credible_sets
n_credible_set_snps
index_snp_in_credible_set
warning_message
error_message
L
coverage
```

Suggested `status` values:

```text
success
no_credible_set
non_converged
susie_error
credible_set_error
```

The report should make it possible to determine:

- how often donors are removed for missing dosages;
- how often SuSiE fails;
- how often it does not converge;
- how often no credible set is reported;
- how often the tensorQTL index SNP is retained in a credible set.

### 15. Write outputs

Write the verbose and slim outputs as tab-separated gzip-compressed files using:

```r
data.table::fwrite(output, outfile, sep = "\t")
```

Use the supplied output paths directly; the caller controls suffixes.

Write the report similarly. Compression should follow the supplied `.gz` suffix when supported by `data.table::fwrite()`.

Create parent directories only if that behavior is documented explicitly; otherwise require them to exist and fail early.

Do not write partial final output files if a fatal validation error occurs before gene-level fitting.

For runtime failures after fitting has begun, build results in memory and write all outputs once processing finishes.

## Error and Warning Policy

### Fatal errors

Use `stop()` for:

- missing required files;
- missing required columns;
- invalid function arguments;
- no significant eQTLs at the requested q-value threshold;
- donor-set mismatch across expression, genotype, and covariate files;
- missing or duplicate significant expression rows;
- missing or duplicate required genotype rows;
- missing or duplicate cis-pair rows;
- index SNP absent from its gene's cis-pair rows;
- disagreement between `num_var` and cis-pair count;
- missing or non-finite expression values;
- missing or non-finite selected covariate values;
- missing requested covariates;
- duplicated selected covariate names;
- monomorphic or zero-variance cis variants;
- insufficient retained donors after genotype-missingness filtering;
- unexpected duplicate gene-SNP rows in either output;
- unexpected membership of one gene-SNP pair in multiple credible sets.

### Recoverable gene-level conditions

Use `warning()` and continue for:

- SuSiE fitting errors;
- SuSiE non-convergence;
- credible-set extraction errors.

Warnings must identify the gene and index SNP.

The same details must be recorded in the report and represented by `NA` SuSiE fields in the verbose output.

A successful fit with no credible set is not an R warning. Record it as `status = "no_credible_set"`.

## Suggested Internal Helpers

```r
construct_dataset_paths()
validate_explicit_input_paths()
read_bed_donors()
read_covariate_donors()
validate_identical_donors()
read_index_eqtls()
read_cis_pairs()
resolve_covariates()
read_expression_rows()
read_genotype_rows()
validate_gene_inputs()
subset_complete_genotype_donors()
residualize_against_covariates()
run_susie_for_gene()
extract_susie_snp_results()
build_verbose_output()
build_credible_set_output()
build_gene_report()
validate_output_keys()
write_finemapping_outputs()
```

Keep `run_eqtl_finemapping_files()` focused on orchestration. Keep `run_eqtl_finemapping()` limited to standard path construction, validation, and delegation.

## Implementation Notes

- Use `data.table::fread()` and `data.table::fwrite()` for large tabular files.
- Avoid `dplyr`.
- Use explicit namespace qualification for third-party functions, including `susieR::susie()` and `data.table::fread()`.
- Use base R or `data.table` for joins, filtering, and reshaping.
- Preserve source tensorQTL columns when they provide useful interpretation.
- Avoid repeatedly reading the genotype BED once per gene; read the union of required variants once.
- Avoid repeatedly reading expression once per gene; read the selected expression rows once.
- Process SuSiE fits sequentially in the first implementation unless parallel execution is added deliberately with deterministic warning and output handling.
- Do not save `.rds` fit files.
- Do not use external LD or proxy-population LD.
- Do not transform or rank-normalize expression beyond the optional linear covariate residualization specified here.
- Record package versions in package metadata; optionally include `susieR` version in each report row or as a report attribute.

## Testing Plan

### Unit tests

Test helpers using small synthetic files for:

- exact path construction in the directory adapter;
- explicit-path validation in the core function;
- confirmation that the adapter delegates to the core function without duplicated analysis logic;
- comment-line handling in `cis_qtl_ann.txt.gz`;
- q-value filtering;
- exact covariate selection;
- `covariates = NULL` selecting all rows;
- donor-set equality independent of donor order;
- donor-column reordering;
- cis-pair selection by gene;
- index SNP membership validation;
- missing-genotype donor removal;
- QR residualization of both `X` and `y`;
- residualization with redundant one-hot encoded covariate columns;
- credible-set row extraction;
- no-credible-set fallback behavior;
- failed-fit fallback behavior;
- one-row-per-gene-SNP output enforcement.

### Integration tests

Create a compact synthetic dataset directory containing all five expected files and verify:

1. one significant gene with a reported credible set;
2. one significant gene with no credible set;
3. one gene with missing dosage values that removes donors;
4. the same variant occurring for two different genes;
5. a forced SuSiE error or mocked failure;
6. exact output schemas and sort order;
7. invisible return value and all three written files.

### Fatal validation tests

Confirm that the run stops for:

- mismatched donor sets;
- missing expression;
- missing covariate values;
- absent index SNP;
- missing cis variant in the genotype BED;
- duplicate gene or variant rows;
- monomorphic genotype rows;
- insufficient donors after genotype filtering.

## Package Quality and Check Requirements

The completed package must pass:

```r
devtools::check()
```

with no errors, warnings, or notes attributable to the package. This requirement is equivalent in intent to passing `R CMD check` on the supported development platform.

Implementation requirements for check compliance:

- Add complete roxygen2 documentation for both exported functions and use `@export`.
- Document all exported arguments, return values, side effects, output files, errors, and warning behavior.
- Declare package dependencies correctly in `DESCRIPTION`.
- Use explicit namespace qualification for third-party functions, for example `data.table::fread()` and `susieR::susie()`, rather than relying on broad package imports.
- Do not use `library()` or `require()` in package code.
- Ensure all referenced functions are either base R functions, explicitly namespace-qualified, imported narrowly, or defined internally.
- Avoid `R CMD check` notes for undefined global variables. For `data.table` code, use forms that are check-safe without creating ambiguous `..name` bindings.
- Do not modify the global environment, working directory, user options, random-number state, or other process-wide state without restoring it with `on.exit()`.
- Use `tempfile()` or test fixtures for files created during tests; tests must not write into the package source tree.
- Include `testthat` tests for the core function, directory adapter, validation failures, fallback outputs, and file writing.
- Ensure examples are either fast and self-contained or marked appropriately when they require external input files.
- Do not include generated output files, local paths, credentials, or large test data in the package source.
- Run `devtools::document()` before the final check so `NAMESPACE` and generated help files are synchronized with the R source.

The implementation is not complete until the final `devtools::check()` succeeds cleanly.

## Acceptance Criteria

The package is complete when it can:

1. process one explicit set of five input files per core-function call;
2. provide a directory adapter that discovers the five standard tensorQTL-derived files from the directory basename and delegates to the core function;
3. select index eQTLs using `qval <= 0.05` by default;
4. identify all cis variants for each selected gene from `cis_qtl_pairs.txt.gz`;
5. assert identical donor sets across expression, genotype, and covariates before fitting;
6. use all tensorQTL covariates by default or an exact caller-supplied subset;
7. stop on missing expression, missing covariates, monomorphic variants, or structural input mismatches;
8. drop donors per gene when any required genotype dosage is missing;
9. residualize both expression and genotype against the same covariate design;
10. run individual-level `susieR::susie()` with configurable `L` and `coverage`;
11. write one verbose row for every analyzed gene-SNP pair;
12. write one slim row for every credible-set gene-SNP pair;
13. emit the index SNP as a fallback slim row when no credible set is reported or fitting fails;
14. permit the same SNP to appear for multiple genes while enforcing one row per gene-SNP pair;
15. continue after gene-level SuSiE errors or non-convergence and record `NA` result fields;
16. write a one-row-per-gene processing report;
17. invisibly return all three output data frames;
18. avoid external LD matrices and saved SuSiE fit objects;
19. pass `devtools::check()` with no package-attributable errors, warnings, or notes.

## Methodological Basis

The implementation follows the individual-level SuSiE fine-mapping workflow:

- fit expression `y` against the donor-by-variant genotype matrix `X`;
- use `L = 10` as the default maximum number of effects;
- use configurable credible-set coverage, defaulting to 95%;
- report SNP posterior inclusion probabilities;
- regress covariates out of both `X` and `y` before fitting;
- use the observed genotype matrix rather than an external LD panel.


## Clarification on Redundant One-Hot Covariates

A full-rank covariate matrix is not required for this workflow. The tensorQTL covariate file may contain a complete set of one-hot encoded levels, making the design matrix intentionally redundant when an intercept is included.

The implementation should retain the selected covariate rows as supplied and use QR-based residualization for both expression and genotype values. The package does not need to calculate or report the covariate rank, rank-deficiency status, or residual degrees of freedom. It should instead validate the resulting residualized expression vector and genotype columns directly, requiring finite values and non-zero variance before calling SuSiE.
