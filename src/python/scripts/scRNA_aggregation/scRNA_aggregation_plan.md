# Plan: Serial shell-script port of ScRnaAggregationWorkflow for Nextflow outputs

## Context
`ScRnaAggregationWorkflow.scala` (Zamboni) aggregates the per-library outputs of the old UGER pipeline: scPred `summary.txt` cell metadata, donor DGEs, and `buildSummaryTable` library metrics. The new Nextflow pipeline writes different artifacts to GCS. scPred is gone and MapMyCells (MMC) replaces it, the cell metadata is `<uei>.cmd.tsv`, and the QC summary is now `summarize_scRNA_experiments.py`. We want a simple serial bash script that reproduces the Zamboni DAG using the new locator and summary tools, so each step can be tested on its own. It will later be the spec for a Nextflow workflow.

Reference run (legacy): `/broad/bican_um1_mccarroll/.../2026-07-20_v192_10X-GEMX-3P_EC/unfiltered/` (the commands come from `logs/workflow.log`).
Test data (new): `gs://mccarroll_scrnaseq_central_temp/nextflow/pipeline-testing/jn_SMRI_V11_rxn1_output/...`

## Facts that shape the design
- `cmd.tsv` uses **lowercase** `cell_barcode` and (village only) `donor`. The legacy files used `CELL_BARCODE`/`DONOR`/`predClass`.
- The locator's `cell_metadata` field already prefers the village `cmd.tsv`, which has donors, over the standard_analysis one. `mmc` is `{model: {mmc_annotations, ...}}` or `NA`.
- The MMC `<uei>.csv` is **comma-separated** with `#` header lines. One of those lines is `# taxonomy hierarchy = ["neighborhood", "class", ...]`. Its key column is `cell_id`, and each level has a `<level>_name` column.
- `join_and_filter_tsv` (`src/dropseq/aggregation/join_and_filter_tsv.py`) reads join files as TSV only. `--set` adds columns at the **end**. All of the Python tools take local files, so GCS inputs are staged with `gcloud storage cp` (the script exports `CLOUDSDK_CORE_PROJECT=bican-um1-mccarroll` unless it is already set, because the gcloud `core/project` is a project number on the Mac).
- `donorCellTypeReport` (`transcriptome/R/packages/DropSeq.aggregation/R/donorCellTypeReport.R`) hardcodes `DONOR` and `predClass`. It is ported to Python (`donor_cell_type_report`), and the R file stays unchanged.
- The CLIs for `MakeTripletDge` (`barcode_column` in its manifest), `cat_tsvs` (`--index-col`) and `MakeMetacellsFromTripletDge` (`GROUP_COLUMNS`) already let us pass the new column names.

## Paired comparison data (same library, both pipelines)
- Zamboni: `/broad/mccarroll/nemesh/zamboni_vs_nextflow/zamboni_out/libraries/2024-03-27_SMRI_V11_rxn1/GRCh38_ensembl_v43.isa.exonic+intronic/std_analysis/svm_nuclei_noCbrbInit_10X/`. It has the donors DGE and `MapMyCells/HMBA_Human_WB_v0.5/*.{csv,tsv}`. There is no scPred `summary.txt` and no aggregation run.
- Nextflow: `gs://mccarroll_scrnaseq_central_temp/nextflow/pipeline-testing/jn_SMRI_V11_rxn1_output/`, with manifests in `/broad/mccarroll/nemesh/zamboni_vs_nextflow/manifest/`.
- In Zamboni, `DropSeq.cellclassification::convertPrefixMMCCsvClp` (`transcriptome/R/packages/DropSeq.cellclassification/R/MapMyCells.R:10`) turns the MMC CSV into a TSV. It reads the CSV skipping `#` lines, prefixes every column except `cell_id` with `<model>_`, and writes it tab-separated. Nextflow publishes only the CSV.

## Input manifest (same schema as Zamboni's `LaunchScRnaAggregation`)
```yaml
dgeDefaults:               # projected onto each dge that does not set the key
  filters: {doublet: {exclude: doublet}}
dges:                      # a single dictionary is also accepted
  - dge: gs://.../<uei>.donors.digital_expression.txt.gz    # required: a donor or non-donor DGE
    donor: N1              # optional, for a library that is not a village
    filters:               # keyed by cell metadata column: include, exclude, includeFile, excludeFile, min, max
      pct_mt: {max: 0.1}
    joins:                 # secondary files joined on to the cell metadata, in order
      - {joinFile: extra.tsv, leftColumn: cell_barcode, joinColumn: cell}
```
- `libraryId` is **not supported**: the library ID is the `uei` of the DGE. Giving it is an error.
- `donor` is handled by `locate_scRNA_artifacts`, which records it in the library's artifact YAML as `donor` (`NA` if not set).
- `filters` and `joins` are handled by `resolve_aggregation_manifest`, a stand-in for the Nextflow filter and join step that goes away with it. Column names must be the Nextflow metadata names (`cell_barcode`, `doublet`, `<model>_<level>_name`, ...), not the scPred ones.

## Donor rule
Every Nextflow `cmd.tsv` already has the right `donor` column (`JOIN_CELL_METADATA` in `mccarroll_nextflow`: from `donor_cell_map` for a village, or a constant for a library run with the `donor` parameter). A library with neither has no `donor` column. So, per library:
1. A `donor` in the manifest or `dgeDefaults` overrides it: `join_and_filter_tsv --set donor X`.
2. Otherwise the `donor` column of `cmd.tsv` is used as it is.
3. Otherwise the library has no donor. `cat_tsvs` leaves the donor empty for its cells, and `MakeMetacellsFromTripletDge` skips cells with an empty group value.

## Decisions (confirmed)
- **Python only.** R is needed only for plotting, if that is ever needed. `donorCellTypeReport` is ported.
- **The script has no inline Python and no Zamboni artifacts** (no `workflow.json`, no `finished.txt`). It only calls tools and reads one TSV samplesheet, so it ports directly to Nextflow (`splitCsv`).
- **Every tool is run from `--dropseq-tools`**, Java and Python, so a deployment pins the whole workflow to one version. The Python tools need a wrapper there: `transcriptome/python/src/dropseq/command_line_programs.tsv` lists them (`summarize_scRNA_experiments`, `locate_scRNA_artifacts`, `resolve_aggregation_manifest`, `donor_cell_type_report`, `convert_mmc_csv`, besides `cat_tsvs` and `join_and_filter_tsv`).
- Keep the new column names: `cell_barcode`, `donor`, and MMC columns prefixed with the model as in Zamboni (e.g. `HMBA_Human_WB_v0.5_class_name`). Only the R report gets a temporary copy with renamed columns.
- The cell-type column is `<model>_<level>_name`. `MMC_LEVEL` is a parameter; if it is not set, use the first level in the CSV's `taxonomy hierarchy` line.
- `MMC_MODEL` is a parameter; if it is not set, each library uses its first model. The script stops if that gives different cell-type columns for different libraries.
- The DGE aggregated for a library is `dge_donors` if it exists, else `dge_selected_cells`. (`user_dge` is the DGE the user gave, and is used only to match manifest entries to artifacts.)
- One YAML per library, from `locate_scRNA_artifacts --output-dir`. The files are flat, so the script and the Nextflow workflow never parse the manifest.

## Tool changes (in `public/src/python`, all with tests; `uv run python -m unittest ...`)
1. **`convert_mmc_csv`** (`src/dropseq/mmc/convert_mmc_csv.py`, registered in `pyproject.toml`): a Python port of `convertPrefixMMCCsvClp`. `--input csv --output tsv --column-prefix <model>_` skips `cell_id`. `read_taxonomy_hierarchy(path)` parses the `# taxonomy hierarchy = [...]` line. The fixture test matches the Zamboni `.tsv`.
2. **`join_and_filter_tsv --set-first`**: `--set` columns come first, in the order given, so `PREFIX` is column 1. The default behaviour does not change.
3. **`locate_scRNA_artifacts`**:
   - new `user_dge` field (the DGE given) and `donor` field (`NA` if not set), both after `uei`;
   - new `reference` field (the `reference` FASTA of the alignment `properties.yaml`) and `reduced_gtf` field (the one `*.reduced.gtf` in the reference directory; `NA` if there is none, an error if there is more than one). The script and the user never name the reduced GTF;
   - reads `dgeDefaults`, and the `donor` of each entry (`read_dge_entries`);
   - `--donor` for a single `--dge`;
   - the command line run fails if two DGEs have the same `uei` (`check_unique_ueis`). `locate_datasets` stays permissive, so several DGE types of one experiment can still be located.
4. **`resolve_aggregation_manifest`** (new): applies `dgeDefaults`, matches each manifest DGE to its artifact file by `user_dge`, and writes the `join_and_filter_tsv` arguments for its joins and filters (one argument per line), plus the samplesheet `libraries.tsv`, with a header and one row per library: `library_id artifact_file dge cell_metadata mmc_model mmc_annotations reduced_gtf args_file` (`NA` for missing values). `dge` is `dge_donors`, else `dge_selected_cells`. `mmc_model` is `--mmc-model` or the library's first model. A donor in the artifacts (from the manifest) becomes `--set donor X` at the start of the arguments file. It checks that no library ID equals the analysis ID, that all libraries have the same `reference` (or alignment directory name if that is `NA`; `--reference` skips the check) and the same `reduced_gtf`, and that each has a DGE and cell metadata, and rejects unknown manifest keys, including `libraryId`.
5. **`convert_mmc_csv --level` and `--cell-type-column-output`**: `--level` defaults to the top level of the taxonomy hierarchy. The output file has the cell type column name `<prefix><level>_name`, checked against the converted columns.
6. **`donor_cell_type_report`** (new, port of `donorCellTypeReport`): `--input cell_metadata --output report --cell-type-column COL`, with columns `donor cell_type num_nuclei fraction_nuclei median_umis_per_nucleus`. Cells with no donor are ignored. Cells with a donor but no cell type count in the donor's total but get no row; R counted them as one extra NA cell type.

## Deliverable layout
The plan and the script go in one directory, because together they are the spec for the later Nextflow pipeline:
- `public/src/python/scripts/scRNA_aggregation/scRNA_aggregation_plan.md` is this plan. Copy it in as the first implementation step, and update it whenever the script changes so the two stay in sync.
- `public/src/python/scripts/scRNA_aggregation/run_scRNA_aggregation.sh` is the script.

## Script: `run_scRNA_aggregation.sh` (bash, `set -euo pipefail`, serial)
**Parameters:** `--manifest`, `--analysis-id`, `--out-dir`, `--dropseq-tools` (directory with all of the tools, Java and Python); optional `--tmp-dir`, `--mmc-model`, `--mmc-level`, `--reference`, `--java-memory`.
**Helpers (inline in the script):** `run` (echoes the command), `fetch` (`gcloud storage cp` for `gs://`, `cp` otherwise) and `has_column`. There is no inline Python.

| # | Step | Command / output |
|---|------|------------------|
| 0 | Set up | Create `OUT_DIR` and `TMP_DIR`. |
| 1 | Locate artifacts, resolve the manifest | `locate_scRNA_artifacts --manifest $MANIFEST --output-dir $OUT_DIR/artifacts` writes one flat `<uei>.yaml` per library. Then `resolve_aggregation_manifest ... --output-dir $OUT_DIR/resolved_manifest` writes the samplesheet `libraries.tsv` and the per-library donor, filter and join arguments. |
| 2 | Library metrics (replaces `buildSummaryTable`) | For each library, run `summarize_scRNA_experiments --artifacts <uei>.yaml --summary artifacts/<uei>.library_metrics.txt --tearsheet-dir tearsheets/`. Then `cat_tsvs` → `<ID>.library_metrics.txt`. |
| 3 | Per-library cell metadata (replaces the scPred step) | For each library, from its samplesheet row: stage `cell_metadata`, `dge` and, if present, `mmc_annotations` into `TMP_DIR/<uei>/`. Run `convert_mmc_csv --column-prefix <model>_ [--level L] --cell-type-column-output cell_type_column.txt` → `mmc.tsv`. Then run `join_and_filter_tsv --input cmd.tsv --set-first --set PREFIX <uei> [--join mmc.tsv cell_barcode cell_id] <args file: [--set donor X] and the manifest joins and filters>` → `<uei>.joined_filtered_cell_metadata.txt`. The MMC columns are appended at the end; cells with no MMC result get empty values. Record the cell-type column from `cell_type_column.txt` (`<model>_<level>_name`). Fail if it differs between libraries. Warn and skip the cell-type steps if no library has MMC. |
| 4 | Write the MakeTripletDge manifest | A heredoc loop that writes, for each library, `dge` (local copy), `prefix`, `barcode_list` (the step 3 output) and `barcode_column: cell_barcode` → `<ID>.make_triplet_dge.yaml`. |
| 5 | Combine expression (unchanged) | Stage the `reduced_gtf` of the samplesheet (the same file for all libraries, checked by the resolver) into `TMP_DIR`, unless it is already there, then `MakeTripletDge MANIFEST=… REDUCED_GTF=… OUTPUT=matrix.mtx.gz OUTPUT_CELLS=barcodes.tsv.gz OUTPUT_FEATURES=features.tsv.gz TMP_DIR=…` |
| 6 | Combine metadata (unchanged) | `cat_tsvs --index-col PREFIX --index-col cell_barcode …` → `<ID>.cell_metadata.txt` |
| 7 | Donor / cell-type report | Only if the combined metadata has a `donor` column and an MMC cell-type column. `donor_cell_type_report --input <ID>.cell_metadata.txt --output <ID>.donor_cell_type.txt --cell-type-column <celltype_col>` |
| 8 | Donor metacells | Only if there is a `donor` column. `MakeMetacellsFromTripletDge … GROUP_COLUMNS=donor` → `<ID>.donors.metacells.txt.gz` |
| 9 | Donor × cell-type metacells | Only if there is a `donor` column and an MMC cell-type column. `GROUP_COLUMNS=donor GROUP_COLUMNS=<level>_name` → `<ID>.donor_cell_type.metacells.txt.gz` |

Each step echoes its command line, so the log reads like the Zamboni `workflow.log`. Cell filters come only from the manifest; with none, the result is "unfiltered", as in the legacy run.

## Out of scope (for now)
- The legacy PDFs (cell_type_counts, donor_qc, etc.) were not produced by this workflow and are not ported. The Nextflow implementation comes later.

## Status
Tool changes 1 to 6 are implemented and unit tested. The script has no inline Python. The script ran end to end, steps 0 to 9, on one village library (SMRI_V11_rxn1, 238 cells, one MMC model) in about 40 seconds, with the Java tools from `/Users/nemesh/jn_branch`. Checked: 17 donor metacell columns, 55 donor × neighborhood columns (equal to the rows of `donor_cell_type.txt`), 238 barcodes, no `workflow.json`, and a second run did not download the reduced GTF again. A second run with `dgeDefaults` filters (`doublet`, `pct_mt`), a manifest `donor` (NX) and `--mmc-level class` also worked: 193 of the 238 cells remained, all with donor NX, and the cell type was the class. **Not yet run:** a library without donors, manifest joins, several libraries, a library without MMC, `--mmc-model`, and a run with all of the tools from a deployment (`/Users/nemesh/jn_branch` needs a rebuild to get the wrappers listed in `command_line_programs.tsv`).

## Example run
```
cd public/src/python
D=/Users/nemesh/nextflow_aggregation_shell_script
bash scripts/scRNA_aggregation/run_scRNA_aggregation.sh \
  --manifest $D/manifest.yaml --analysis-id SMRI_V11_test --out-dir $D/out --tmp-dir $D/tmp \
  --dropseq-tools /Users/nemesh/jn_branch 2>&1 | tee $D/run.log
```
`--dropseq-tools` must contain every tool the script calls. Until the deployment is rebuilt, use a directory of links to the `.venv/bin` Python tools and small `exec` wrappers for the Java tools (a Java wrapper cannot be symlinked, because it finds `loadDotKits.sh` next to itself).

## Running and developing
- **Run the script from a deployment** (`--dropseq-tools /Users/nemesh/jn_branch`), with plain `bash`, no `uv`: the script calls every tool from that directory, so this tests the code Nextflow and other users will run. The deployment needs a rebuild to get the wrappers listed in `command_line_programs.tsv`.
- **Use `uv run` in `public/src/python` only for development:** unit tests and trying a tool on the working tree (an editable install in `.venv`).
- The script exports `CLOUDSDK_CORE_PROJECT=bican-um1-mccarroll` unless it is already set.

## Next steps
1. The user rebuilds `jn_branch`, clears the old outputs and reruns the example above, to confirm that every tool comes from the deployment.
2. The user runs the Zamboni aggregation workflow on the Zamboni results of the same libraries, and compares: cells per donor, donor metacell correlation, and cell type fractions at the same MMC level. Expected differences: column names, scPred class against MMC level, the report's NA cell type counts, slightly different selected cells, and the library metrics columns. Zamboni is a sanity check, not the target.
3. Run the untested cases (listed in Status).
4. Commit on `jn_aggregation_workflow` (and the `command_line_programs.tsv` change in the `transcriptome` repo) when asked.
5. Build the Nextflow version, replacing the resolver and the manifest parsing with Nextflow-native code.

## Verification
1. Unit tests, from `public/src/python`: `uv run python -m unittest discover -s tests/agregation -t .` and `uv run python -m unittest tests.mmc.test_convert_mmc_csv` (`pytest` is not installed in the uv environment).
2. Converter sanity check: the Zamboni `convertPrefixMMCCsvClp` is only a reference for the *converter's* column-prefixing semantics, and the fixture test already matches it. Beyond that, Zamboni output is not the target. The target is the Nextflow implementation, whose experiment and cell metadata have changed, missing and added columns compared with Zamboni (e.g. lowercase `cell_barcode`/`donor`, no scPred columns, scDblFinder `doublet`, new library-metrics columns).
3. Run the script on a manifest of the SMRI_V11 village donor DGEs (rxn1, and any other reactions available) with a local `OUT_DIR`. Then check that:
   - each `<uei>.joined_filtered_cell_metadata.txt` has `PREFIX` as column 1, the same row count as `cmd.tsv` (with no manifest filters), and the MMC columns (from `cell_id`) at the end;
   - the number of rows in `barcodes.tsv.gz` equals the number of rows in the combined `cell_metadata.txt`;
   - the metacell files have one column per donor and per donor × `<model>_<level>_name`;
   - `library_metrics.txt` has one row per library.
4. Also run it with a library that is not a village (no donor, and with a manifest `donor`), and with a manifest `filters` and `joins` entry.
5. Judge the outputs against the Nextflow artifacts, not against Zamboni. Check that every column in a library's `cmd.tsv` survives into `cell_metadata.txt`, that the added columns are only `PREFIX` and the MMC columns, and that the barcodes and donors agree with the Nextflow `donor_cell_map.txt`. Differences from the legacy v192 directory are expected and are not failures; just note which columns changed.
