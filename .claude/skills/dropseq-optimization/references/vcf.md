# VCF/BCF path: `ParallelVCFReader` (`vcftools/`)

```java
try (ParallelVCFReader r = new ParallelVCFReader(path, requireIndex, numThreads)) {   // also (path,numThreads), (path,requireIndex,numThreads,batchSize,queueDepth)
    VCFHeader header = r.getFileHeader();
    Iterator<VariantContext> it = r.iterator();     // once only
}
```
- One reader thread cuts lines (or BCF records) in file order, N workers each with their own codec decode batches and **force full genotype decode**, and the consumer takes `Future`s from a bounded queue in submission order. Same records and order as `VCFFileReader`; everything downstream stays single-threaded.
- `numThreads == 1` delegates to plain `VCFFileReader` (lazy genotypes). There is deliberately no parallel+lazy mode (the lazy parser would race on the worker's codec buffers).
- Defaults: batch 100, queue depth 2 x threads (1x is slower; 1000-record batches add head-of-line blocking). Memory ~ (depth+2) x batch records (a 345-sample record is ~32 KB).
- Formats: `.vcf`, `.vcf.gz`, uncompressed BCF 2.1 (`.bcf`). BGZF-compressed BCF is rejected until htsjdk 6. BCF decodes ~8x faster than VCF.gz but is ~13x larger (watch disk throughput). INFO floats may re-format from BCF (numerically equal).
- Scales ~linearly to ~8 threads on pure decode; the ceiling after that is the caller's thread (genotype-walking filters, `SNPInfoCollection` building) or the reader thread (~15-20x).
- `requireIndex=true` opens/closes a `VCFFileReader(vcf,true)` so a missing index fails as before.

## Wiring recipe
Done for `GatherDigitalAlleleCounts` and `AssignCellsToSamples` (both passes) via `SampleAssignmentVCFUtils.getSNPInfoCollection(..., numThreads, log)` and `getVCFIterator(vcfReader, records, ...)`: add `NUM_THREADS`, keep a `VCFFileReader` only for header/sample checks, iterate records from `ParallelVCFReader`, **close it in `finally`** (walks end early or throw), add an end-to-end test at 1 and 4 threads, then one full-size timing with byte-identical output comparison.

Candidates: whole file in order AND touches genotypes (not `query()`/random access). Left: `DetectDoublets` (strong), `CensusSeqUtils.getVCFIterator` (CensusSeq, RollCall), `FindMonomorphicSitesInDonorPool`/`CsiAnalysis`, `PrepareEqtlGenotypeData` (full-file path), `CreateSnpIntervalFromVcf` (if GQ/het filters on).

If the caller's thread is the limit after threading the reader: do per-SNP work once per SNP, not once per cell, and move optional output (`VCF_OUTPUT`) to a background writer (done for ACTS at `NUM_THREADS > 1`, htsjdk `Options.USE_ASYNC_IO`).

Parked, not built: region sharding with one `VCFFileReader` per thread (needs an indexed VCF, per-shard filter chains, ordered merge of partial `SNPInfoCollection`s). Only if one-thread decode becomes the ceiling.

Tools that read both a BAM and a VCF: one visible `--NUM_THREADS`; hidden `BAM_DECODE_THREADS` (default 1 = stock reader); extra BAM threads come out of N (VCF reader gets `max(1, N - (B-1))`). Give the BAM reader a thread only when a test shows it pays; the VCF usually benefits more.
