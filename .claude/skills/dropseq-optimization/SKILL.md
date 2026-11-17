---
name: dropseq-optimization
description: Speed up a Drop-seq Java tool (BAM/VCF readers, multithreading, codecs, writers, hot spots). Use when asked to optimize, profile, add --NUM_THREADS to, or benchmark a command in this repo. Holds the reusable code paths, thread-balancing rules, what was already tried and dropped, and the verification and workflow rules.
---

# Drop-seq optimization

Repo: `dropseq/` module, Java 21, htsjdk 5.0.1, Picard 3.5.0, TestNG, Barclay args. Plans, results and data paths of finished work live in `~/*.md` plan files and git history, not here.

## 1. Triage: what bounds the tool?

Profile before coding (JFR `settings=profile` plus per-thread CPU, e.g. `jdk.ThreadCPULoad#period=200ms`). A top frame on the calling thread often shows only where a record's memory is first touched; confirm a suspect with a timing experiment. When a frame surprises you, read the stack below it (`jfr print --stack-depth`) before ranking options on it. `jdk.ExecutionSample` does not see native deflate/inflate: count `jdk.NativeMethodSample` frames (`Deflater.deflate`, `Inflater.inflate`) to size the codec share before choosing a writer/reader fix. Measure distributions the idea depends on (group sizes, rows per cell) on the slice before quantifying a gain; rows are not proportional to cells.

| Bound by | Signature | Lever |
|---|---|---|
| BAM read (inflate/decode) | little per-record work | `SamFileMergeUtil` threads: [bam-io](references/bam-io.md) |
| BAM write (encode + deflate) | tool rewrites BAM, writer thread ~100% | codecs, batched writer: [bam-io](references/bam-io.md) |
| VCF decode | scans whole VCF, touches genotypes | `ParallelVCFReader`: [vcf](references/vcf.md) |
| Calling thread (per-record Java work) | one core at 100%, workers idle | move work to workers: [threading](references/threading-patterns.md) |
| Sort / SortingCollection | sorting dominates | not solved by readers. First ask whether the sort is needed at all: [hotspots](references/hotspots-and-dropped.md) |

Order that has paid off every time: (1) codec flags (no code), (2) wire the reader, (3) fix the writer, (4) move per-record work off the calling thread, (5) micro-optimize only what the profile shows. Threading before hot spots.

## 2. Rules (all tools)

- `--NUM_THREADS`: visible, default 1, validated >= 1 in `customCommandLineValidation`. 1 = today's code path and output, untouched. Hidden `READER_THREADS`/`WORKER_THREADS`-style args (given together, then `NUM_THREADS` is ignored) allow sweeps; log the resolved split at INFO.
- Read through `SamFileMergeUtil` only. Do not modify `UnsortedMergingSamRecordIterator` or `OrderedConcurrentMapper`, and do not use `BamDecodePool` directly in commands. Ask before touching shared `SamFileMergeUtil`.
- Do not hand-write what a public library API does (gzip/BGZF framing, CRC, codecs, format bytes). Magic numbers, bit shifts and byte offsets in new code are a warning sign. Prefer htsjdk/Picard over a third-party library it wraps. When two library routes exist, build both, measure, let the user choose.
- Public repo: keep public methods (delegating wrappers instead of deleting).
- Code layout: `doWork` directly after the `@Argument` fields; helpers and inner classes at the bottom of the file. Tuning-only args are `@Hidden`. New files get the current-year (2026) MIT header. Tests live in the matching class under `src/tests/java`.
- Avoid try/catch/finally unless required (errors crossing threads that would hang a consumer, checked exceptions the API forces). No "output is identical for any thread count" sentences in `@Argument` docs.
- Local runs always pass `--USE_JDK_DEFLATER true --USE_JDK_INFLATER true` (Mac has no Intel codec; this selects libdeflate). GCP tuning and compression level are the user's/team's call: report, do not change.
- Threads must be cost efficient. Production uses even thread counts. Report wall, user CPU, wall x N, peak RSS (`/usr/bin/time -l`) for every threaded change, and judge the marginal gain per doubling. Production sizes VMs by memory (80% rule), so the RSS number matters.
- Large refactors: discuss, write a plan to `~/<tool>_tuning_plan.md` (self-contained: paths, commands, baselines, decisions, a "resume here" section), get approval before coding. Update it after each step.
- Never run Gradle or a benchmark while the user is timing something, or while another benchmark runs (edit files, build when it finishes). No commit or push unless asked; the user opens PRs. When a change has a separable part (e.g. a data-structure swap next to the threading), keep it a separate commit so it can be dropped.

## 2b. Keeping this skill current

If development turns up a general-purpose optimization rule (one that applies beyond the tool at hand), propose adding it here, but only after asking the user and getting permission. Do not edit the skill on your own. Put tool-specific results, numbers and data paths in the tool's plan file, not here, and fold a new rule into the existing section it belongs to rather than appending a duplicate.

## 3. Verification and benchmarking

See [verification](references/verification.md): slicing, identity checks, test rules, run recipes.

## 4. Reference index

- [bam-io](references/bam-io.md): parallel reader, codecs, BAM writers, ordered mapper, parallel gzip text output.
- [vcf](references/vcf.md): `ParallelVCFReader`.
- [threading-patterns](references/threading-patterns.md): consumer-side counting, shard-by-cell workers, reader/writer balance, memory.
- [hotspots-and-dropped](references/hotspots-and-dropped.md): single-thread fixes that paid off, ideas tried and dropped, unbuilt ideas, Nextflow follow-ups.
- [verification](references/verification.md)
