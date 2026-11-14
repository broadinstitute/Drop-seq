# Verification, tests and benchmark recipes

## Commands
```
cd /Users/nemesh/dropseqrna/transcriptome_java/public
./gradlew --offline :dropseq:installDist -q                 # builds dropseq/build/install/dropseq
./gradlew --offline :dropseq:test --tests '*Name'           # targeted; full suite ~70 s
JAVA_OPTS=-Xmx8g dropseq/build/install/dropseq/bin/dropseq <Program> --USE_JDK_DEFLATER --USE_JDK_INFLATER ...   # JAVA_OPTS and DROPSEQ_OPTS reach the JVM
```
Results: `dropseq/build/test-results/test/TEST-*.xml` (`tests= skipped= failures= errors=` on the `<testsuite>` line). Inspect library APIs with `javap -p -cp dropseq/build/install/dropseq/lib/htsjdk-5.0.1.jar <class>` (no sources jar). Edit with Edit, not `sed -i`.

## Tests (few, critical, end to end)
- Favor integration tests that run the real code: real files, the tool's `doWork`, read the output back as a consumer would. No boundary/size permutations or tests of unimportant parts.
- Per threaded tool: the same end-to-end case at 1 and 4 threads (compare to the same expected output), a multi-input case when relevant, and `NUM_THREADS < 1` rejected via `customCommandLineValidation` (protected: the test sits in the tool's package). `@Argument(minValue = 1)` is not enforced by the parser.
- Mutation-check each new test (break the code, e.g. deliver batches out of order or remove backpressure, and see it fail). Use the build exit status and delete the old result XML first: a failed compile leaves old XML that looks like a pass. If mutating a value fails no test, add the output assertion first.
- Compare data, not bytes: decompressed text, records, header minus `@PG`. Out-of-order hook for deterministic completion order: a package-private `beforeBatchDecode(int)`-style callback (see `ParallelVCFReaderTest`, `ParallelBamReaderTest`).

## Slices and baselines
- Iterate on a slice, finish on the full data. Slice by cell barcode, not by reads from the top of the file, so all alignments of a read stay together: `samtools view -D XC:<barcodes> -b`, taking every 4th barcode of the `BamTagHistogram` output (`count<TAB>barcode`). Use a slice with ALL read counts for anything touching cell selection. Slices of the benchmark BAM: `/downloads/SELECT_CELLS_BY_NUM_TRANSCRIPTS/slice25` and `slice25all` (`barcodes25.txt` beside the BAM).
- Slice a BAM via SAM text (`samtools view -h | awk`, or `head -n` on SAM lines), never `head` the BAM stream. `samtools view -l` is a library filter and silently drops reads. `samtools view -u` makes level-0 BAMs.
- Baseline the OLD code with a git worktree of master built separately (`git worktree add <scratch>/base HEAD`, build there), so old and new run from the same inputs and flags. Remove the worktree afterwards.

## Benchmarks
- Compare one variable at a time; read `@Argument` defaults and the log echo before comparing runs. Run-to-run noise is 5-15%; repeat and rotate order, state whether the page cache is warm. Slice runs: differences under ~15 s are noise.
- Launch each benchmark as its own background task: run the script itself with `run_in_background` (not `nohup &` in a shell), put outputs in the user's monitoring directory (`/downloads/SELECT_CELLS_BY_NUM_TRANSCRIPTS/...`) when asked, and check the output files when reporting status. Don't poll.
- Record wall, user CPU, wall x N, peak RSS (`/usr/bin/time -l`). Production reads the central `workflow.log` (per-step logs can be truncated); `/tmp` is noexec on UGER.

## Output identity
- Text/`.gz`: `cmp <(gzip -dc a.gz) <(gzip -dc b.gz) && echo IDENTICAL`, or `gzip -dc f | md5` on each in parallel. Use `gzip -dc`, not `zcat`, on macOS. Exclude `#` command/header lines (`grep -v '^#'`).
- BAM: `samtools view` records identical, header minus `@PG`, metrics/manifest minus `#` lines. At full size do NOT `cmp` two `samtools view` pipes (far too slow): run `samtools view -@3 f | md5` on each BAM in parallel and compare.
