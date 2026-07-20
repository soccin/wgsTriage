# The background distribution: how it is gathered and computed

Status: work in progress. This describes the importer as it stands on
2026-07-20 and the state of the background files currently committed. It is
written to be read before rebuilding the background, and section 6 is the part
that matters if you are about to do that.

Companion documents: `docs/METHODS.md` (where the thresholds come from),
`docs/DECISIONS.md` (design rules). The producer is
`bin/wgsTriageBackground.R`; the shared parsers and threshold definitions are
in `R/qcLib.R`; the consumer is `bin/wgsTriage.R`.

---

## 1. What the background is for

`wgsTriage.R` applies fixed thresholds to a new cohort. The thresholds alone
answer "is this sample outside the acceptable range". They do not answer "how
far outside, compared to what we normally see", and that second answer is what
makes a FAIL defensible to someone upstream who does not want to hear it. The
background supplies it: a robust central value and spread for each metric,
computed over previously mapped samples that passed every threshold they could
be checked against.

Consumption is narrower than production. `wgsTriage.R` reads only
`backgroundStats.tsv`, and from it only two columns:

- `median`, used by `foldOf()` to express a value as a multiple of normal
  ("90x higher than normal" rather than "12.0%"),
- `n`, quoted in the report so a reader can see how much evidence a
  comparison rests on.

The quantiles, the MAD, the per-class coverage table and the metric-coverage
table are computed and two of them are committed, but nothing reads them
programmatically. They exist for human review of the background itself.

---

## 2. Inputs

The importer takes one argument, an archive directory of previously mapped
projects, and scans it recursively. Three file kinds are recognised by name:

| Kind | Pattern | Reader | Contributes |
|---|---|---|---|
| alignment | `.asm.txt`, or Picard's `.alignment_summary_metrics` | `readAsmMetrics` (qcLib) | `pctChimeras`, `pctSoftclip`, `pctReadUsed`, `pctImproperPairs`, read lengths |
| coverage | `.wgs.txt`, `.wgs_metrics`, `.raw_wgs_metrics` | `readWgsMetrics` (qcLib) | `meanCoverage`, `pctExcOverlap`, `pctExcTotal`, `pctExcDupe` |
| samtools | `multiqc_samtools_stats.txt` | `readMultiqcSamtoolsAny` (local) | `supplementaryRate`, `pctProperlyPaired`, `insertSizeAverage`, `interChromRate`, `pctMapped` |

No directory layout is assumed. The archive genuinely contains the same cohort
under `<proj>/out/metrics`, `<proj>/Map/out/metrics`,
`<proj>/results/r_002/mapping/metrics` and, where an `rsync -R` copy was made,
buried under a replica of the source absolute path.

The archive is read only. Nothing outside the tree is followed: manifest files
listing paths elsewhere are ignored by design, and `dir_ls(recurse = TRUE)`
does not descend through symlinked directories.

---

## 3. The pipeline

Six stages, in execution order.

### 3.1 Discovery

One recursive traversal, then classification by filename regex. Two identifiers
are derived per file, and both are load-bearing because together they form the
join key.

**Sample** — `sampleFromFile()` strips the metric extension, then strips
pipeline stage suffixes (`STAGE_SUFFIX`: `smap`, `md`, `recal`, `sorted`,
`dedup`, `markdup`, `bqsr`, `final`) repeatedly to a fixed point, because they
stack as `.md.recal`. Only known tokens are removed. Sample names in this
archive legitimately contain dots, so truncating at the first dot would be
wrong.

**Project** — `projectFromPath()` walks *up* from the file, discarding path
components that match `SCAFFOLD_DIR` (pipeline scaffolding: `out`, `metrics`,
`mapping`, `multiqc`, `r_002`, `Set1`, and so on) or that exactly equal the
sample name or the file name, and takes the first component that remains. This
replaced a fixed sequence of string edits on the path tail, which assumed the
wrapper directories always appeared at the same depth. They do not. Walking up
is stable under every layout seen; string surgery on the tail was not.

The same label must come out for a sample's alignment, coverage and samtools
files, or those files will not join into one row. Section 6.1 is about
verifying that.

### 3.2 Parsing

`collectMetricFiles()` wraps each reader in `safely(quietly(...))`. Errors and
warnings are both captured per file and neither reaches the console: a
malformed file produces both, and letting them print buries the summary under
vroom noise that does not say which file caused it. The audit says which file.

It returns two things, and the separation matters: `data` (parsed rows) and
`audit` (one row per file *discovered*, including the ones that failed). A
background that quietly shrank is indistinguishable from a background that was
always small, so every file discovered is accounted for in the audit whether or
not it contributed.

Two shapes are reconciled here. A Picard file names one sample, in its
filename. A multiqc table names many samples, internally. The file-derived name
is carried separately and a name from inside the file wins, so both shapes
emerge with one unambiguous `sample` column.

`readMultiqcSamtoolsAny()` is a deliberate local replacement for qcLib's
`readMultiqcSamtools`, which is still used by `wgsTriage.R` for current
cohorts. Two differences:

- qcLib keeps only rows suffixed `.recal`. Across this archive most samtools
  rows carry `.md` with no `.recal` counterpart at all, so that rule discarded
  the large majority of the samtools evidence and left `supplementaryRate` —
  the tightest metric in the filter set — resting on a fraction of the samples
  it could have used. The local reader prefers `.recal`, falls back to `.md`,
  then to anything, and records which in a `samtoolsStage` column.
- qcLib discards the whole file if any one of six wanted columns is absent.
  The local reader materialises a missing column as all-NA, so a table lacking
  insert size still contributes a usable supplementary alignment rate.

### 3.3 Duplicate resolution

Within each kind, collapse to one row per `(project, sample)`. **Newest mtime
wins**, ties broken on path so repeated runs agree. Every collision is written
to the audit with disposition `supersededDuplicate` and a warning is raised.

Note what is *not* claimed here: mtime is a proxy for recency, not provenance.
Section 6.2.

### 3.4 Joining

`joinOnProjectAndSample()` inner-joins on `(project, sample)` and keeps
anything unmatched on either side as its own row. It is applied twice:
alignment against coverage, then the resulting Picard frame against samtools.

**There is deliberately no fallback join on sample name alone.** An earlier
draft had one, on the theory that a differing project label meant one run filed
in two places. Picard headers disproved it: `NK_KHYG1_CL_D`, for example, has
alignment metrics from a `.smap.bam` against `human_g1k_v37_decoy` in November
2025 and coverage metrics from an `.md.cram` against b37 in February 2026.
Those are two separate mappings of one biological sample, and merging them
invents a sample whose chimera rate and whose coverage were never measured on
the same alignment.

Nothing is lost by leaving them apart, because reference ranges are built per
metric from whatever is not NA — both rows contribute exactly what they
measured. Names appearing under more than one project label are reported as a
warning for a human to adjudicate (section 6.3).

Samples that only ever appeared in a multiqc table are kept. They carry real
samtools metrics and contribute to the ranges that are already thinnest.

After the joins the key is checked for duplication rather than assumed, since a
silent duplicate would double-weight a sample in every range.

### 3.5 Threshold evaluation

Every threshold metric is materialised, as an all-NA column where the archive
has no such data. Without this, a metric absent everywhere is never evaluated
and therefore never reported as MISSING, which reads exactly like a metric that
was checked and passed.

Evaluation then uses qcLib's `evaluateThresholds()` and `sampleVerdict()` —
the same functions the report applies to a current cohort, so the background
carries verdicts on the same terms. Those functions key on a single column, so
they are passed a project-qualified `uid` of the form `project::sample`. Keying
on the bare name would merge verdicts across unrelated runs; roughly a hundred
names recur across projects in this archive.

### 3.6 What counts as a reference sample

```
referenceSample = nFail == 0 & nEvaluable > 0
referenceTier   = "full" if both pctChimeras and meanCoverage are present,
                  otherwise "partial"
```

The gate is **what a sample failed**, not which files it happens to have. An
earlier rule also required both Picard files present. That sounds conservative
and is not: a sample missing one of its two files then contributed to no range
at all, not even the ranges built from the file it did have. On an archive
where alignment and coverage metrics were collected into different subtrees,
that rule reduced 401 usable coverage samples to a coverage range built from 2.

`referenceTier` keeps the distinction visible rather than resolving it, and the
full/partial ratio is the single most useful indicator of how well vetted a
given background is. Section 6.5.

---

## 4. Outputs

Six files, all overwritten on every run.

| File | Content | Names | Committed |
|---|---|---|---|
| `backgroundStats.tsv` | per-metric n, median, MAD, q01/q05/q25/q75/q95/q99, min, max over reference samples; long form, NA dropped per metric | no | yes |
| `backgroundCoverageStats.tsv` | the same for `meanCoverage`, split by `sampleType`, since tumors are sequenced far deeper than normals | no | yes |
| `backgroundMetricCoverage.tsv` | per-metric `nAvailable` / `nTotal` / `pctAvailable` across **all** samples, not just reference ones | no | yes |
| `backgroundSamples.tsv` | every sample, every metric, verdict, provenance columns | yes | no |
| `backgroundFlagged.tsv` | samples with `nFail > 0` or `nWarn > 0` | yes | no |
| `backgroundImportAudit.tsv` | one row per file discovered, with `parsed` and `disposition` | yes | no |

`.gitignore` denies `data/background/*` by default and allows the three
name-free aggregates by name, so any further file the importer starts emitting
stays out of the repository until it has been reviewed. That is the PHI guard,
not an oversight.

Two details that look like bugs and are not:

- `backgroundFlagged.tsv` filters on `nFail`/`nWarn`, not on `verdict`.
  `sampleVerdict()` returns INCOMPLETE the moment any metric is missing, and
  since `pctReadUsed` is absent from nearly every archived file that describes
  nearly the whole archive. Filtering on `verdict == "WARN"` reported zero
  warned samples across the entire archive — in the one file whose purpose is
  to list them. The console counters use the counts for the same reason.
- Ranges are computed per metric over non-NA values, so different rows of
  `backgroundStats.tsv` rest on different, and sometimes very different,
  numbers of samples. The `n` column is not decoration.

The console summary prints, in order: the sample and reference counts, the
input-file audit, the samtools stage mix, the worst historical samples (naming
the metrics that actually failed, since which metrics exist depends on which
files the archive holds), the ranges with per-metric archive coverage, and a
collected warnings block. Warnings are accumulated through `warnImport()`
rather than `warning()` because R prints deferred warnings after the summary,
where they are routinely missed.

---

## 5. Current state of the committed background

The files under `data/background/` were built from the **QCDataV2** archive
(confirmed from the paths in the audit). Headline numbers:

```
samples            621        projects  17
reference samples  579        of which "full" tier:  2
below threshold     42
```

Per-metric availability across all 621 samples:

| Metric | n in range | % of archive |
|---|---|---|
| supplementaryRate, pctProperlyPaired, insertSizeAverage, interChromRate | 470 | 82.4 |
| meanCoverage, pctExcOverlap, pctExcTotal, pctExcDupe | 374 | 64.6 |
| pctChimeras, pctSoftclip, pctImproperPairs | 77 | 12.4 |
| pctReadUsed | 11 | 1.8 |

**This background is wide but weakly vetted, and that should be stated wherever
its numbers are quoted.** Only 2 of 579 reference samples carry both Picard
files. The archive holds 77 usable alignment-metrics samples against 374
coverage samples, so the great majority of samples entering the coverage,
insert-size and supplementary ranges were never checked for chimeras or
soft-clipping at all — they pass the reference gate because those metrics are
MISSING rather than because they are clean. The chimera range itself rests on
77 samples, and `pctReadUsed` on 11.

The samtools stage mix in this build is 413 `md`, 99 `recal`, 109 with no
samtools row.

Two project labels in the audit look wrong on inspection: `Umich10_Umich10_T`
(a sample name that reached the label position) and `ReMap_260130`, which
carries 271 of the 506 files and is a remapping batch rather than a cohort.
Neither is fatal — the label only needs to be *consistent* across a sample's
files to join correctly — but `ReMap_260130` means half the archive is one
undifferentiated bucket, which weakens any per-project reasoning about it.

The earlier `QCData` archive has substantially more alignment metrics (order
400 asm files against QCDataV2's 88 discovered / 77 parsed). Whether QCDataV2
was meant to be the smaller archive is an open question for whoever assembled
it; if it was not, rebuilding from `QCData`, or from the two combined, would
produce a far better vetted background.

---

## 6. What to check before trusting a NEW background

Ordered by likelihood of trouble.

### 6.1 Verify the derived project labels first

`SCAFFOLD_DIR` is a heuristic list. Patterns like `b[0-9]+`, `r_?[0-9]+`,
`set[0-9]+`, `batch[0-9]+`, `normal[0-9]+`, `down[a-z0-9_]*` will happily eat a
cohort genuinely named that way (over-walking: two cohorts merge into one
label). Conversely, a new wrapper directory not on the list stops the walk
early (under-walking: one cohort splits into several labels). Over-walking is
the more dangerous of the two, because it can merge unrelated runs under one
key.

```r
read_tsv("<out>/backgroundImportAudit.tsv") |> count(kind, project) |> print(n = 100)
```

Symptoms: a label that is a sample name, a filename, or a bare batch token; two
known cohorts collapsed into one row; one known cohort spread over several.

### 6.2 mtime is not provenance

Duplicate resolution keeps the newest mtime. `rsync -a` and `cp -p` preserve
mtimes; plain `cp` does not. An archive assembled by plain copy can present the
older data as the newest. Inspect every collision:

```r
read_tsv("<out>/backgroundImportAudit.tsv") |> filter(disposition == "supersededDuplicate")
```

On QCDataV2 this correctly kept `results/r_002/mapping` (2026-06-11) over
`Map/out` (2026-05-26) for 11 Proj_17495_I samples.

### 6.3 Adjudicate cross-project name collisions

When the importer warns that a name appears under more than one project label,
decide whether it is one run filed twice or two genuinely separate mappings.
The Picard header settles it: it records the input BAM, the reference and the
run date.

```bash
grep -m1 -E 'INPUT|Started on' <sample>.asm.txt
grep -m1 -E 'INPUT|Started on' <sample>.wgs.txt
```

If input BAM, reference or date disagree, they are separate mappings and must
stay separate rows. This is exactly how the name-only join was ruled out.

### 6.4 Stage mixing in the samtools metrics

Ranges blend `.recal` and `.md`. Measured on the samples carrying both,
`.recal` runs below `.md` by at most 0.055 points of supplementary rate, 0.52
of properly-paired percent, and 1.3 of insert size — around 3% of the 1.0%
supplementary threshold, so it moves nothing near a boundary in the current
archive. A new archive with a different stage mix will shift the ranges. Check
the stage counts in the console block and the `samtoolsStage` column.

Do not "fix" this by reverting to `.recal` only. Every elevated-supplementary
cohort in the current archive exists as `.md` rows alone, so a `.recal`-only
background contains no failing sample at all, and the threshold ends up
validated against evidence that excludes the thing it is meant to catch.

### 6.5 Watch the full/partial ratio

With `referenceTier == "partial"`, a sample with only coverage metrics enters
the coverage range without ever being checked for chimeras. Read the two
reference-sample lines in the console summary together:

```
Reference samples   579
                      2   of those carry both Picard files
```

A ratio like that one is a legitimate import and a weak background. Say so in
any report that quotes the ranges.

### 6.6 Contamination fraction

The design assumption is that defective cohorts are excluded *by the
thresholds*. A defective cohort that passes every threshold silently widens the
ranges, and robust statistics are not sufficient protection when the
contaminated fraction is large and concentrated in one cohort. Review
`backgroundFlagged.tsv` against expectation: 42 failing samples in the current
build.

### 6.7 Thin metrics

A range built from a handful of samples is not a range. Always read
`backgroundMetricCoverage.tsv` alongside `backgroundStats.tsv`. `pctReadUsed`
was absent archive-wide because older Picard lacks `MEAN_ALIGNED_READ_LENGTH`;
QCDataV2 supplied the first 11 samples for it, so this resolves as new projects
accumulate rather than needing a fix.

### 6.8 Smaller traps

- Symlinked directories are not followed. An archive assembled from symlinked
  cohort directories imports as empty.
- Manifest files are ignored by design; the importer scans the tree only.
- `docs/METHODS.md` quotes concrete numbers from the background. Regenerating
  it invalidates them; check them again after any rebuild.

---

## 7. Known weaknesses and what to improve

Roughly in order of how much they matter.

1. **The committed background is built from the weaker archive.** 2 of 579
   reference samples are fully vetted. The first improvement is not a code
   change: confirm which archive is meant to be authoritative and rebuild from
   it. Combining `QCData` and `QCDataV2` would need the duplicate and
   cross-project machinery to be re-checked, since the same cohorts appear in
   both, but the machinery exists for exactly that case.

2. **`wgsTriage.R` reports `max(refN)` as the background size.** One line
   in the console summary quotes the largest per-metric `n` as "reference
   ranges from N previously mapped samples". With the current build that
   prints 470 while the chimera comparison beside it rests on 77. It should
   quote the per-metric `n`, or a range, not the maximum.

3. **No provenance record in the outputs.** Nothing in `backgroundStats.tsv`
   records which archive it came from, when, or with which version of the
   importer — the archive identity above had to be recovered from the audit
   file, which is gitignored and therefore absent from a clean checkout. A
   header comment or a small `backgroundProvenance.tsv` carrying the archive
   path, run date, `WGSTRIAGE_VERSION` and the sample counts would close this.

4. **The quantiles and the per-class coverage table are unused.** They are
   computed, and two of the three committed files are never read by anything.
   Either the consumer should use the spread (a fold-change against the median
   says nothing about how tight the distribution is; q05-q95 or the MAD would),
   or the outputs should be trimmed to what is actually consumed. The former is
   the better change: `foldOf()` currently treats a metric with a 6x spread the
   same as one spanning 0.09 to 0.16.

5. **`referenceTier` is recorded but not acted on.** Nothing weights or
   segregates partial samples, and no output surfaces the tier breakdown per
   metric. At minimum, `backgroundMetricCoverage.tsv` could carry an
   `nFullTier` column so the thinness of the vetting is visible per metric
   rather than only in aggregate.

6. **`SCAFFOLD_DIR` is a heuristic with no test.** The failure mode is silent
   and it changes the join key. A small fixture tree exercising over-walking
   and under-walking, asserting the expected labels, would lock this down; it
   is the piece most likely to break on a new archive layout.

7. **Two samtools readers now exist.** `readMultiqcSamtoolsAny` in the
   importer and `readMultiqcSamtools` in qcLib differ in stage handling and in
   tolerance of missing columns, so the background and a current cohort are not
   read on identical terms. This was a deliberate scope constraint during the
   rewrite, not a design decision. The qcLib reader should probably adopt the
   tolerant behaviour and the local copy be deleted.

8. **Robustness is proven but not tested.** A hostile-archive generator
   (zero-byte, truncated, binary, unreadable, permission-denied, duplicate,
   dotted names, paths with spaces, multiqc with no `.recal`/`.md`, missing
   columns, empty) was used to demonstrate graceful degradation, but it was a
   scratch script and is not in the repository. Promoting it into `tests/`
   would keep that behaviour from regressing.

9. **The importer is around 840 lines in one file.** Discovery, parsing,
   joining, evaluation, aggregation and console rendering are all in it. It is
   readable, but the rendering block in particular has no reason to live beside
   the join logic.

---

## 8. Rebuilding

```bash
Rscript bin/wgsTriageBackground.R <QCDir> [--out <OutDir>]
Rscript bin/wgsTriageBackground.R --help
```

Default output is `data/background` under the repository root, resolved from
the script's own location. Expect one to two minutes over a roughly 1.4 GB
archive.

After a rebuild:

1. Work through section 6, starting with the project labels.
2. Re-run the threshold tests: `Rscript tests/testThresholds.R`.
3. Confirm the consumer still works against the new background:
   `Rscript bin/wgsTriage.R tests/fixtures/miniCohort --background <OutDir> --out <TmpDir>`.
4. Re-check the concrete numbers quoted in `docs/METHODS.md` and in section 5
   of this document.
