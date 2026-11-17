# Hot spots, dropped ideas, open ideas

## Single-thread hot spots that were worth fixing (apply only when the profile shows them)
- `htsjdk Histogram<String>.increment` is TreeMap-backed (~20 string compares per increment over millions of keys). Count in a `HashMap<String,long[]>`/array and build the `Histogram` once at the end (output order is TreeMap order anyway; double sums of whole counts are exact).
- `rec.getReadString()` + substring: use `rec.getReadBases()` with `BaseRange.getSequenceForBaseRange(List<BaseRange>, byte[])`; build Strings lazily. Small gain, only when not inflate/deflate-bound.
- Exact/ED1 barcode lookup: `beadsynthesis/PackedBarcodeSet` (2-bit packed longs, open addressing; falls back to Strings for N/lowercase/odd lists). Keep candidate order identical.
- Interval annotation: range-fill `getLocusFunctionsByInterval` and `int[]` allocation rework in `TagReadWithGeneFunction`.
- `PrintStream`/`ErrorCheckingPrintStream` over a raw stream sends every `println` straight to the underlying stream: put a 1 MiB `BufferedOutputStream` between them (largest single gain for a text-merge tool). Writing to `/dev/null` through an unbuffered `PrintStream` is a misleading "no gzip" control.
- UMI edit-distance collapse (`UMICollection.getDigitalExpression`, used by every `UMIIterator` consumer): `HammingDistance.withinHammingDistance(a, b, max)` compares in place and stops past `max`; a group with one UMI skips the collapse. About 40% of (cell, gene) groups have one UMI, 17% two, 43% three or more.
- Pairwise UMI collapse is O(n^2) per (cell, gene) group and group sizes are very skewed (10x: groups over 100 UMIs are 0.3% of groups but ~80% of compares). `collapseBarcodesGeneric` also paid `ArrayList.removeAll` per core even with nothing to remove (skip when empty, use a Set). For Hamming ED 1, `PackedEd1Collapser` (utils.editdistance, on `PackedBarcodeSet`) looks up the 3 x length single-base changes instead; used from `MapBarcodesByEditDistance.collapseBarcodes` when ED is 1, no indels and #barcodes >= 3 x barcode length (cost model, no tuned constant). N, mixed length or >31 bases fall back to the pairwise search, so correctness never depends on the data. Slice 76.6 -> 48 s; DigitalExpression 11%.
- Do not count reads in a separate pass when the main pass already reads everything (`SelectCellsByNumTranscripts` counts reads per cell inside the UMI pass; saved a whole prefilter pass, 4 of 32 min).
- Compute shared lists once per SNP, not per cell; avoid re-filtering records pass 1 already filtered (verify the premise first: a cheap position pre-filter failed because pass 1 keeps far more SNPs than expected).
- Equal records must not be ordered by `HashMap` iteration: tie-break explicitly (e.g. by barcode, last).
- Text tools reusing htsjdk `ProgressLogger` print a meaningless `Last read position: :0`: subclass it, override `log(String...)` and strip that suffix (`UMICollectionByCellParser.CellProgressLogger`). Count and label what the tool actually iterates, ideally as a fraction of the cell list.
- A text-merge tool with output compression off the calling thread is capped by that thread's parse + `MergingIterator` + join; extra compression threads add nothing.

## Tried and dropped (do not redo without new evidence)
- Parallel sort of each in-memory chunk in `CombiningSortingCollection` (2 and 4 threads): wall -0.3% / -1.5% for +30% CPU.
- Rebuilding `CombiningSortingCollection` on htsjdk `SortingCollection` with a `TreeMap` pre-combine: +4% wall. A plain `SortingCollection` merging duplicates at iteration would sort/spill ~2.2x the records (not run).
- A `GeneFunctionProcessor.processReads` shortcut for one-record groups (~85% of groups): ~0.6 min. The profile blamed `TreeMap.put` on the resolver; it was htsjdk's `SortingCollection$MergingIterator`.
- Name sort / key sort chunk sizes (500k / 2M): not tried, user not interested.
- Early row filtering to `CELL_BC_FILE` in the text parser as a speedup: selected cells hold ~85% of rows (excluded cells have far fewer UMIs), ~2% gain. Kept only for the log.
- Replacing `progress.record(rec)` with position-only; byte-based barcode extraction as a standalone speedup (no effect when I/O-bound).
- Likelihood math in ACTS is ~0% of time; do not thread it. A cheap position pre-filter in ACTS pass 2 was dropped.
- SelectCellsByNumTranscripts options D, C, the htsjdk rebuild, N and M were tried or dropped by decision.

## Still unaddressed (ideas, not built)
- Parallel BGZF inflate within one file; user ruled out relying on BGZF internals (public htsjdk APIs only). Parallel BGZF BAM writer: see bam-io.
- Sorted text-merge tools (`MergeMolecularBarcodeDistributionByGene`): one reader/parser thread per input feeding the same `MergingIterator` in the same input order.
- Sort-bound tools (`BamTagOfTagCounts`, `CompareDropSeqAlignments`, `SpermSeqMarkDuplicates`): sorting stays serial. Sorting writers (`presorted=false`) keep a `SortingCollection` in front of any parallel writer.
- `SplitBamByCell` calling-thread cost above 4 threads: queryname merge comparator (shared code, ask first), per-record `XC` decode (mates adjacent, could skip the second), HashMap lookup.
- Writer-bound tools needing the batched writer: `GenerateSyntheticDoublets`, `DetectBeadSynthesisErrors` (cleanBAM), `CollapseBarcodesInPlace`, `FilterBam`, `DetectBeadSubstitutionErrors` (repair), `FilterBamByGeneFunction`, `DiscoverMetaGenes`.
- `SingleCellRnaSeqMetricsCollector` after the sort removal: GTF load (~16 s on the slice run, 1.7 GB file) and the extra read-count pass when `--NUM_CORE_BARCODES` is used.

## Nextflow follow-ups (user's repo `mccarroll_nextflow`, not edited by agents)
A threaded program needs more than `process_single`: raise CPUs, pass `--NUM_THREADS ${task.cpus}` and the two codec flags (one shared `params.codecArgs`). Note `-resume` cache invalidation and byte-different BAMs when codecs/levels change. Suggested CPUs so far: `countBarcodeSequences` ~4, `correctScrnaReadPairs` 2-3, `splitBamByCell` ~4, `mergeMolecularBarcodeDistributionByGene` `--NUM_THREADS 2` (~3 CPUs).
