# Methods

Thresholds as of `wgsTriage` 0.9.0. Any report produced by a different version
may have been judged against different numbers; `--help` prints the version of
the tool that wrote it.

How `wgsTriage` decides that a sample is unusable, and why each number is what
it is. Every threshold here was set against an archive of previously mapped
projects rather than chosen in the abstract, and the observed ranges that
justify them are given alongside.

Read this before changing a threshold.

---

## 1. What is read

Three files per cohort, all of them already written by the Map stage:

```
<MapDir>/out/metrics/<sample>/<sample>.asm.txt    Picard CollectAlignmentSummaryMetrics
<MapDir>/out/metrics/<sample>/<sample>.wgs.txt    Picard CollectWgsMetrics
<MapDir>/sbam/multiqc/multiqc_data/multiqc_samtools_stats.txt
```

Nothing is recomputed and nothing is written back. The assessment is therefore
free in compute terms and cannot perturb the pipeline it inspects.

Three parsing decisions matter enough to state, because each one has already
caused a wrong answer:

**Picard files are parsed structurally, never by line offset.** The parser
locates the `## METRICS CLASS` marker, takes the header row after it, and stops
at the first blank line before the optional histogram block. Picard's schema
drifted between the version that produced the archive and the one running now:
`MEAN_ALIGNED_READ_LENGTH` is absent from all 421 archived alignment-metrics
files and present in current output. A fixed `skip=` would have silently read
the wrong column.

**The `PAIR` row is the one that is used.** `FIRST_OF_PAIR` and
`SECOND_OF_PAIR` each describe half the data. Reading either instead is a quiet
way to halve every count.

**Only per-sample multiqc rows are used.** `multiqc_general_stats.txt` is
per read group, keyed `SAMPLE-FLOWCELL_LANE`; only `multiqc_samtools_stats.txt`
carries clean per-sample rows. The `.recal` rows are taken because that is the
final BAM state and the same object Picard measured, so the two sources are
directly comparable.

### Why qualimap is excluded

Qualimap is not read at all. On degraded samples it overstates usable coverage
by up to 2.6x, and in Proj_16840_N it silently reported on only 12 of the 16
samples. An analysis built on it produced a conclusion that had to be
retracted. The failure mode is the dangerous one: it does not error, it reports
a smaller cohort as though it were the whole cohort.

This is also why cohort completeness is counted from the number of sample
directories under `out/metrics`, not from the number of samples successfully
parsed. A sample whose mapping died leaves an empty directory. Counting only
what parsed is precisely how a report covers 12 of 16 samples without anyone
noticing.

---

## 2. The seven filter thresholds

| Metric | What it physically measures | Source | Warn | Fail |
|---|---|---|---|---|
| `pctChimeras` | Read pairs whose two ends align to distant loci. Directly indicates broken or artefactually joined template molecules. | Picard asm | above 1.0% | above 5.0% |
| `supplementaryRate` | Reads aligning in multiple pieces. The same physical defect as chimerism, measured by a different tool. | samtools | above 0.5% | above 1.0% |
| `pctSoftclip` | Bases the aligner declined to align at read ends. Rises when read ends do not match the reference. | Picard asm | above 1.0% | above 5.0% |
| `pctReadUsed` | Fraction of each read the aligner actually used. Falls when reads are only partly alignable. | Picard asm | below 98% | below 95% |
| `pctProperlyPaired` | Pairs in the expected orientation and separation. | samtools | below 98.5% | below 97.0% |
| `pctExcOverlap` | Bases discarded because mates overlap. High values mean inserts shorter than twice the read length. | Picard wgs | above 10% | none |
| `pctExcTotal` | Bases discarded from coverage for any reason. The summary measure of how much sequencing was wasted. | Picard wgs | above 35% | none |

The two `pctExc` metrics warn but never fail on their own. They describe
consequences rather than causes, and a high value is informative only alongside
one of the others.

These seven are the fixed thresholds. An eighth check, the coverage floor,
varies by sample class and so is defined separately; see section 5. A report
states the total, which is eight when every source is present.

---

## 3. How the thresholds were set

The derivation is in `docs/NORMAL_BAM_QC_REPORT.md`, which is the original
analysis of the Proj_16840_N failure and the source of the separation argument
below. `docs/DECISIONS.md` records what was chosen and what remains
unvalidated. This section summarises both; those two are what to read if a
number here is challenged.

### The chimera threshold, and why it moved

The original analysis proposed failing above 1.0% chimeric pairs, on the
grounds that affected and unaffected samples separated by a factor of 45 with
nothing in between. That is true of the 27 samples it examined. It is not true
across the 421 archived samples carrying alignment metrics:

| Band | Samples | Reading |
|---|---|---|
| below 0.5% | 283 | the normal population |
| 0.5 to 1% | 107 | still normal, upper tail |
| 1 to 5% | 9 | mildly degraded, cause unknown |
| above 5% | 22 | the catastrophic class |

There is a real, sparse middle. A 1.0% fail threshold rejects 31 of 421
historical samples, or 7.4% of all previous work. A 5.0% threshold rejects 22,
or 5.2%, and surfaces the nine borderline samples as warnings for review.

FAIL is therefore set at 5.0% and 1.0% became a warning. Nothing that mattered
was lost: every sample in the Proj_16840_N failure ran between 8.5% and 17.3%
and fails outright at either setting. The reasoning is that a filter which
rejects 7% of all historical work gets switched off, and a switched-off filter
is exactly the failure this tool exists to prevent.

This is the single most arguable number in the tool. To be noisier rather than
more permissive, set `fail` for `pctChimeras` and `pctSoftclip` back to `1.0`
in `R/qcLib.R`; it is one number each.

### The observed clean ranges

The remaining thresholds are anchored on what clean samples actually do:

| Metric | Observed in clean samples | Threshold | Margin |
|---|---|---|---|
| `supplementaryRate` | 0.094 to 0.389% | fail above 1.0% | about 2.6x above the observed maximum |
| `pctProperlyPaired` | 97.52 to 99.60% | fail below 97.0% | defective samples ran 72.6 to 96.8% |
| `pctChimeras` | below 0.5% for 283 of 421 | fail above 5.0% | about 5x |

`supplementaryRate` is still the tightest metric available anywhere in the
input. The clean range quoted here is wider than the 0.094 to 0.160% reported
in earlier versions of this document, and the margin correspondingly smaller.
That is a correction, not a drift: the earlier figure was computed over 98
samples that came almost entirely from one cohort, because the importer read
only multiqc rows suffixed `.recal` and discarded the `.md` rows that most
cohorts emit instead. The range across all 11 cohorts is genuinely wider.

### Where the background is thin, and the tool says so

Two honest gaps, both reported in the output rather than hidden:

`pctReadUsed` has no historical background at all. Every archived alignment
metrics file predates `MEAN_ALIGNED_READ_LENGTH`, so this threshold runs on its
fixed value with no distribution behind it. This resolves itself as new
projects accumulate.

The samtools background is no longer thin: 512 of 621 samples, 82.4%, now carry
samtools metrics. It used to be 98 of 454, and the difference is a fixed importer
rather than new data: the reader kept only multiqc rows suffixed `.recal` and
dropped the `.md` rows that most cohorts emit instead.

This matters for more than sample count. Every cohort with an elevated
supplementary rate exists in this archive as `.md` rows alone, so the old
samtools background contained no defective sample at all, and those thresholds
rested on the clean range plus a margin with no observed second population to
point at. Clean samples now run 0.094 to 0.389% and the archive holds samples
up to 15.4%, so the separation is visible in the data rather than assumed.

One caveat that comes with the fix: `.recal` and `.md` describe the same library
one pipeline step apart, and on the 99 samples carrying both, `.recal` runs
lower by up to 0.055 points of supplementary rate, 0.52 of properly paired
percent and 1.3 of insert size. The reference range therefore blends two BAM
stages. The offset is around 3% of the 1.0% supplementary threshold and moves no
sample near a verdict boundary, and `backgroundSamples.tsv` records which stage
each sample came from in `samtoolsStage`.

Reference ranges are computed from the 579 samples that failed none of the
thresholds they could be evaluated on; only 2 of those carry both Picard files.
Requiring both files, as this once did, excluded a sample from every reference
range whenever it was missing one of them, including the ranges built from the
file it did have. Including defective samples would widen the range enough to
admit the next bad cohort, so the gate is what a sample failed. Robust
statistics alone are not sufficient protection when the contaminated fraction is
above 5% and concentrated in a single cohort.

One caveat while the background is mid-rebuild: the committed archive (QCDataV2)
is rich in samtools metrics but thin on alignment metrics, only 77 samples. The
chimera figures above, the 421-sample bands and the ReMap finding, come from the
richer QCData archive rather than the shipped stats. The two converge when the
background is rebuilt; see `docs/BACKGROUND.md` section 7 and the punch list.

---

## 4. Corroboration

This is the section that answers "how do you know?"

`pctChimeras` comes from Picard. `supplementaryRate` comes from samtools. They
measure substantially the same physical defect, fragmented or artefactually
joined template, but they are computed by different tools, on different passes
over the data, from different intermediate files. They are not two views of one
calculation.

When both fire on the same sample, the verdict does not rest on a quirk of one
tool's definition. In Proj_16840_N they agreed on every failed sample. That
agreement is what makes the result defensible to someone upstream who would
prefer to hear otherwise.

The corollary is enforced in the output. When a cohort has no multiqc data, the
samtools thresholds cannot be evaluated, and the report states explicitly that
every verdict rests on a single tool. A reduced check set is never presented as
a full one.

---

## 5. Escalation and pair rules

### Three warnings become a failure

A sample raising three or more warnings is failed even when no single threshold
fires. These metrics are not independent: genuine degradation pushes all of
them at once. Without this rule a sample can sit just below every threshold on
every metric and still be unusable.

The number three is a judgement call, not a derived value. Across the whole
archive it promotes exactly one sample, which carried 3.78% chimeras with five
simultaneous warnings. It is `WARN_ESCALATION` in `R/qcLib.R`.

### Verdict precedence

Evaluated in this order: any failing threshold gives `FAIL`; three or more
warnings gives `FAIL`; any unreadable metric gives `INCOMPLETE`; any warning
gives `WARN`; otherwise `PASS`.

`INCOMPLETE` outranks `WARN` deliberately. A metric that could not be read has
not been checked, and reporting it as clean is how a partial report becomes a
wrong report.

### Pair rules

Tumor and normal are inferred from sample names by pattern. Names that carry no
recognisable marker return `unknown` and are displayed as `?`. They are never
silently treated as normal, and unpaired samples are named in the report.

A pair fails if either member fails, or if the two insert size distributions
differ by more than 1.5x. Facets requires comparable insert size distributions
between tumor and normal; when they diverge it produces a result that is wrong
without crashing, which is the worst available outcome.

Two consequences worth stating plainly. A pair can fail when both members
individually pass. And when every normal in a cohort fails, no pair can be
analysed as matched, including pairs whose tumor is entirely clean. That is
what happened in Proj_16840_N: 8 of 8 pairs were unusable.

The 1.5x tolerance is a starting point and has not been confirmed against the
Facets authors. Insert size is taken from samtools as a mean over
properly-paired reads; a median would be preferable for a heavy-tailed metric,
and the full distribution is available in the samtools `IS` block if that
becomes necessary. The mean already separates the cohorts cleanly: defective
samples ran 108 to 274 bases against 342 to 411 for clean ones.

### Coverage floors

30x for normals, 80x for tumors, and 80x for samples whose class cannot be read
from the name. The unclassifiable case takes the stricter of the two floors: an
unreadable name is a reason to look, not a reason to apply the more forgiving
threshold.

The floor is a warning, not a failure. It counts toward the three-warning
escalation like any other, so coverage can contribute to a `FAIL` but never
causes one alone. It is stratified by class because tumors are sequenced far
deeper than normals and a single number would be wrong for one of them.

It cannot live in `THRESHOLDS`, which carries one fixed value per metric, so it
is computed before the verdicts and appended to the threshold results as a
`meanCoverage` row. The same dropping rule applies as to the other checks: a
cohort with no coverage figures at all is not told once per sample that
coverage is missing.

Adopting these floors moved no verdicts. In Proj_16840_N seven samples fall
below floor and all seven already fail on other thresholds. Across the 401
archived samples that carry a coverage figure, none gains a third warning from
the coverage check, so none escalates.

The floors fire on 12.2% of those 401 samples, against 1.0% for the 25x/50x
pair they replaced. The unclassifiable group is the largest, 266 of 401: the
80x floor flags 32 of them (12.0%), where a 30x floor would flag one (0.4%).

---

## 6. Validation

| Cohort | n | Result |
|---|---|---|
| Proj_16840_N | 16 | 10 fail: all 8 normals plus the MDA001 and MDA002 tumors |
| Proj_17495_I | 11 | 0 fail, the negative control |
| Proj_17608 | 96 | 0 fail |
| ReMap_260130 | 268 | 23 fail |

Reproduce with:

```
Rscript bin/wgsTriage.R <MapDir> --out <OutDir> --project <Name>
```

Proj_16840_N reproduces the original finding exactly, including which samples
fail and in what order. Proj_17495_I and Proj_17608, 107 samples between them,
produce no false positives.

The 23 failures in ReMap_260130 are not a validation result. They are a finding:
that cohort predates Proj_16840_N, shows the same defect class at 7 to 22%
chimeric pairs, and whether it was analysed and released is an open question.

Re-run the first two after moving the repository. The only thing that changes in
a move is paths, and paths are what break.

---

## 7. What this cannot see

Stating the limits is what makes the rest credible.

**Sample swaps and mislabelling.** Detecting these requires genotype
concordance across samples, which needs more than mapping metrics. This is a
common real failure mode and the tool is blind to it.

**Contamination between samples.** Same reason. Cross-sample contamination
needs allele-level evidence.

**Anything requiring the reads themselves.** The tool reads summary metrics, not
BAMs. It cannot assess anything the Map stage did not already measure.

**Whether a defect matters for your specific analysis.** The thresholds describe
data integrity. Whether a given level of degradation invalidates a particular
downstream question is a judgement the tool does not make, which is why it
always exits 0 and leaves the decision to a person.

**Cause.** The tool characterises the data. It does not and cannot say whether a
defect originated in the sample, the library preparation, the sequencing, or the
mapping. Every output is worded to avoid implying otherwise.
