# Defective BAMs in Proj_16840_N — findings and pre-QC specification

**Written:** 2026-07-18 (revised — data source corrected, see section 0)
**Subject:** ten of sixteen BAMs delivered for Proj_16840_N are unusable for
structural-variant calling; how to detect this class of defect before compute is
committed; and how to report it so a non-quantitative reader can act on it.

**Audience:** two. Sections 1-7 are internal (BIC). Sections 8-9 are templates
for outward-facing communication and are written in plain language.

---

## 0. Revision note — read this first

An earlier version of this report drew its metrics from **qualimap**. That was
wrong in two ways, and both errors mattered.

**Coverage was overstated.** Qualimap's mean coverage counts every aligned base,
including duplicates, overlapping mate pairs, and low-MAPQ reads. Picard's
`MEAN_COVERAGE` reports usable depth after those exclusions. Because
`PCT_EXC_TOTAL` reaches 41-73% in the affected samples, the two differ by up to
2.6-fold:

| Sample | Qualimap mean | Picard MEAN_COVERAGE | PCT_EXC_TOTAL |
|---|---|---|---|
| APTL_MDA008_N01 | 64.0x | **25.0x** | 63.8% |
| APTL_MDA010_N01 | 38.8x | **16.8x** | 59.9% |
| APTL_MDA008_T01 | 132.0x | **110.5x** | 22.6% |

**Qualimap covered only 12 of 16 samples, and the gap hid a finding.** MDA001 and
MDA002 have Picard metrics but no qualimap output. **Both of their tumors are
affected** (12.91% and 9.75% chimeric pairs). The earlier conclusion that the
defect was confined to normals was therefore false.

**Process lesson, worth carrying forward:** verify that a metrics source covers
the entire cohort before drawing a pattern from it. The tumor/normal split looked
clean and mechanistically tidy, and it was an artifact of a partial denominator.

**All figures in this report now come from Picard**
(`CollectAlignmentSummaryMetrics`, `CollectWgsMetrics`) at
`Proj_16840_N/Map/out/metrics/*/`. Insert size figures, where cited, remain from
qualimap — see 5.5.

---

## 1. Summary

Ten of sixteen samples in Proj_16840_N carry a technical artifact: all eight
normals, plus the tumors of MDA001 and MDA002. Six tumors (MDA007-012) are clean.

The affected samples show 8.5-17.3% chimeric read pairs against 0.09-0.19% in
unaffected samples and in an 11-sample reference cohort. This killed the Delly
structural-variant caller 90 times out of 90 in run `20260714_130053_12559`.

Because all eight normals are affected, **no tumor/normal pair in the project is
usable for somatic SV calling**, including the six with clean tumors.

**Most useful single metric:** Picard `PCT_CHIMERAS`. It is already computed for
every sample as part of routine processing. Nothing new needed to detect this.

---

## 2. The evidence

### 2.1 Affected samples (n=10)

| Sample | T/N | PCT_CHIMERAS | PCT_SOFTCLIP | Aligned len | MEAN_COV | PCT_EXC_OVERLAP | PCT_EXC_TOTAL |
|---|---|---|---|---|---|---|---|
| APTL_MDA009_N01 | N | **17.28** | 13.69 | 130.3 | 32.53 | 20.90 | 54.31 |
| APTL_MDA002_N01 | N | **17.08** | 32.41 | 102.1 | 20.47 | 27.92 | 63.71 |
| APTL_MDA001_N01 | N | **14.78** | 17.14 | 125.1 | 27.93 | 27.54 | 50.56 |
| APTL_MDA011_N01 | N | **13.76** | 11.00 | 134.4 | 30.57 | 22.82 | 41.29 |
| APTL_MDA012_N02 | N | **13.05** | 13.46 | 130.7 | 28.88 | 25.66 | 44.90 |
| APTL_MDA010_N01 | N | **13.00** | 16.04 | 125.8 | **16.82** | 20.49 | 59.93 |
| APTL_MDA001_T01 | **T** | **12.91** | 22.74 | 116.6 | 25.41 | 16.82 | 71.83 |
| APTL_MDA008_N01 | N | **12.02** | 20.72 | 119.5 | 25.03 | 22.97 | 63.81 |
| APTL_MDA002_T01 | **T** | **9.75** | 30.27 | 105.2 | 23.19 | 19.98 | 72.88 |
| APTL_MDA007_N01 | N | **8.51** | 7.91 | 139.0 | 39.04 | 16.05 | 53.86 |

### 2.2 Unaffected samples, same project (n=6)

| Sample | T/N | PCT_CHIMERAS | PCT_SOFTCLIP | Aligned len | MEAN_COV | PCT_EXC_OVERLAP | PCT_EXC_TOTAL |
|---|---|---|---|---|---|---|---|
| APTL_MDA008_T01 | T | 0.19 | 0.17 | 150.7 | 110.50 | 2.88 | 22.60 |
| APTL_MDA012_T01 | T | 0.17 | 0.17 | 150.7 | 109.68 | 1.92 | 22.03 |
| APTL_MDA009_T01 | T | 0.14 | 0.16 | 150.8 | 114.39 | 2.55 | 21.74 |
| APTL_MDA011_T01 | T | 0.14 | 0.17 | 150.7 | 108.27 | 2.74 | 23.64 |
| APTL_MDA010_T01_D | T | 0.14 | 0.17 | 150.7 | 109.73 | 2.44 | 22.09 |
| APTL_MDA007_T01_D2 | T | 0.11 | 0.19 | 150.7 | 105.15 | 6.50 | 24.91 |

### 2.3 Reference cohort — Proj_17495_I (n=11, all clean)

| | PCT_CHIMERAS | PCT_SOFTCLIP | MEAN_COV | PCT_EXC_OVERLAP | PCT_EXC_TOTAL |
|---|---|---|---|---|---|
| Normals (n=5) | 0.11 - 0.13 | 0.20 - 0.22 | 43.7 - 46.9 | 1.74 - 3.43 | 10.5 - 12.5 |
| Tumors (n=6) | 0.09 - 0.12 | 0.15 - 0.16 | 109.3 - 117.0 | 1.58 - 2.79 | 11.2 - 12.7 |

This cohort completed all 45 Delly tasks with a peak memory of 11.4 GB.

### 2.4 Read length — previously unverified, now resolved

`MEAN_READ_LENGTH` is **151** for every sample in both projects. The earlier
version of this report assumed 150 and flagged it as unconfirmed. The assumption
was close enough that no conclusion changes.

`MEAN_ALIGNED_READ_LENGTH` is the more informative figure: 150.7-150.8 in
unaffected samples (essentially the whole read aligns) versus 102.1-139.0 in
affected ones. Between 12 and 49 bases per read are being discarded.

---

## 3. Effect on the pipeline

75 of 90 Delly tasks were OOM-killed, 13 hit the 2 h wall-clock limit, 2 are
indeterminate. Every OOM died inside `Paired-end and split-read scanning`, the
phase that accumulates chimeric evidence. SLURM confirms:

```
error: Detected 1 oom_kill event in StepId=1738063.batch.
```

`.command.trace` is zero bytes in all 90, which is why `peak_rss` is empty in the
Nextflow trace — the process was killed before the trace was flushed. Memory was
raised 16 -> 32 GB across retries without success, against a 22.6 GB maximum
observed across 565 successful Delly tasks on healthy data.

### 3.1 Why more memory would not have rescued the run

The somatic filter discards any variant with support in the matched normal:

```
-e "FORMAT/DV[1] > 0 | FORMAT/RV[1] > 0"
```

With all eight normals affected, spurious normal support is expected at a large
number of loci, removing genuine somatic events. A run with generous memory would
likely have completed and reported very few SVs, with no indication of a problem.

**The loud failure mode here is the fortunate case. The silent one is the risk.**

### 3.2 Impact beyond Delly — copy number

**Facets requires closely matched insert size distributions between tumor and
normal.** These pairs differ roughly two-fold. Copy number output from these
matched pairs should not be trusted, and unlike Delly there is no crash to signal
it — Facets runs to completion and produces a result.

The established workaround, per the Facets authors, is to substitute an
**unmatched** normal with a comparable insert profile. A well-matched unmatched
normal is preferable to a badly mismatched matched normal.

Consequences for pre-QC design: insert size concordance **within each pair** must
be its own gate (5.5), and reports must carry pair-level as well as sample-level
verdicts.

---

## 4. What is *not* wrong — confounders ruled out

| Hypothesis | Verdict | Evidence |
|---|---|---|
| BAM size | **No** | 242 GB here vs 245 GB in 17495_I, which completed |
| Read group count | **No** | 140-160 here vs 240 in 17495_I |
| Delly version / container | **No** | identical `cmopipeline-delly-bcftools-0.0.3.img` |
| Sequencing depth | **No** | affected samples span 16.8-39.0x; depth does not track the artifact |
| Genuine biology (real SVs) | **No** | see 5.6 |
| Memory config alone | **Partially** | the 160 -> 16 GB cut was real, but 32 GB also failed and healthy data peaks at 22.6 GB |
| Aligner flag difference (`-B 3`) | **No — withdrawn** | `-B 3` on tumors is intentional pipeline design to accommodate higher tumor mutation rates. All prior normals were processed identically. Not causal. An earlier version raised this in error. |

---

## 5. Detection — metrics and thresholds

All primary metrics come from Picard files **already produced by the pipeline**.
No new computation is required for any sample that has been through Map.

### 5.1 Primary gate: PCT_CHIMERAS

Source: `*.asm.txt`, `PAIR` row.

| | Value |
|---|---|
| Unaffected (n=17, both projects) | 0.09 - 0.19% |
| Affected (n=10) | 8.51 - 17.28% |
| Separation | **45x** |
| **Threshold** | **FAIL above 1.0%** |

1.0% sits ~5x above every unaffected sample and ~8x below every affected one.
Anything from 0.5% to 5% produces identical verdicts on all 27 samples.

### 5.2 Confirmatory gate: PCT_SOFTCLIP

Source: `*.asm.txt`, `PAIR` row.

| | Value |
|---|---|
| Unaffected | 0.15 - 0.22% |
| Affected | 7.91 - 32.41% |
| **Threshold** | **FAIL above 1.0%** |

Independent of 5.1 and equally clean. Both firing together is strong confirmation;
only one firing warrants a closer look.

### 5.3 Confirmatory gate: aligned read fraction

`MEAN_ALIGNED_READ_LENGTH / MEAN_READ_LENGTH`.

| | Value |
|---|---|
| Unaffected | 0.998 - 0.999 |
| Affected | 0.676 - 0.921 |
| **Threshold** | **FAIL below 0.95** |

The most directly interpretable of the three: it is the fraction of each
sequenced base that survives alignment.

### 5.4 Warning: PCT_EXC_OVERLAP and PCT_EXC_TOTAL

| Metric | Unaffected | Affected | Threshold |
|---|---|---|---|
| PCT_EXC_OVERLAP | 1.58 - 6.50% | 16.05 - 27.92% | WARN > 10% |
| PCT_EXC_TOTAL | 10.5 - 24.9% | 41.3 - 72.9% | WARN > 35% |

Warnings rather than gates because the 16840_N clean tumors sit at 22-25%
`PCT_EXC_TOTAL` against 10-13% in the reference cohort — elevated but usable.
`PCT_EXC_OVERLAP` is the metric most directly tied to the mechanism.

### 5.5 Cross-source check: delivered vs usable coverage ratio

Requires the coverage figure reported at delivery alongside Picard
`MEAN_COVERAGE`. Delivered figures were obtained for four samples:

| Sample | Delivered | Picard MEAN_COVERAGE | Ratio | PCT_EXC_TOTAL |
|---|---|---|---|---|
| APTL_MDA002_T01 | 52x | 23.19 | **2.24** | 72.88 |
| APTL_MDA001_T01 | 54x | 25.41 | **2.13** | 71.83 |
| APTL_MDA002_N01 | 39x | 20.47 | **1.91** | 63.71 |
| APTL_MDA010_N01 | 28x | 16.82 | **1.66** | 59.93 |

**Ranking by ratio reproduces ranking by `PCT_EXC_TOTAL` exactly.** This is
arithmetic, not coincidence: any metric applying fewer exclusions than Picard
diverges from it by a factor bounded above by `1 / (1 - PCT_EXC_TOTAL)`.

| Cohort | PCT_EXC_TOTAL | Max possible ratio |
|---|---|---|
| 17495_I reference | 11.6% | 1.13x |
| 16840_N clean tumors | 23.5% | 1.31x |
| 16840_N affected | 59.9 - 72.9% | 2.50 - 3.69x |

Observed ratios (1.66-2.24) fall below those ceilings, so the delivery metric
applies *some* exclusions — but materially fewer than Picard. We cannot determine
the method from four numbers and do not need to.

**Gate: WARN above 1.5x.**

**What this check does and does not add.** It provides no detection sensitivity
beyond `PCT_EXC_TOTAL`, which already captures the same signal from our own data
and is gated in 5.4. Its value is different and still substantial:

1. It expresses the defect in a quantity both parties measure, which makes the
   discrepancy checkable rather than assertable.
2. It is the clearest available explanation of *why* delivery QC did not flag
   these samples — the reported metric is inflated by the very defect it should
   reveal, and inflated most in the worst samples.

**Caveats.** n=4, one project, delivered figures rounded to whole numbers. This is
the only check requiring an external input, so it cannot be fully automated
unless delivery figures arrive machine-readable. Treat it as an analyst-facing
diagnostic rather than a hard pipeline gate.

**Worth obtaining:** delivery coverage for the six clean tumors. If those fall
within ~30% of our measurements as predicted, it establishes the divergence as
sample-specific rather than a systematic method difference, which materially
strengthens the case.

### 5.6 Pair-level gate: tumor/normal insert size concordance

**Distinct from all the above and applied per pair, not per sample.** Two samples
can each pass individually and still be unusable together.

Facets requires closely matched insert distributions; where they diverge, the
Facets authors' guidance is that an unmatched normal with a comparable profile is
preferable to a badly mismatched matched normal.

| | Value |
|---|---|
| Healthy pairs (17495_I) | tumor and normal within ~30 bp |
| 16840_N pairs | 132-212 vs 324-397, ~2x divergence |
| **Proposed check** | **FAIL if pair insert medians differ by more than ~1.5x** |

**Two caveats.** The 1.5x figure is a starting point, not validated — confirm the
tolerance against Facets documentation or the authors, since it is theirs to
define. And **insert size is not in the two Picard files currently produced**;
`CollectInsertSizeMetrics` is not being run. Either add it to the Map stage or
take the figure from qualimap. Adding it to Picard is preferable, so the whole
gate reads from one authoritative source.

### 5.7 Why this is not biology

Recorded because it is the first hypothesis anyone raises, including for a
tumor/normal label swap.

- **The affected set does not follow the T/N split.** Two tumors are affected
  alongside all eight normals, and for MDA001 and MDA002 both pair members are
  affected. No labelling error produces that.
- **Soft-clipping and mate overlap are not biological.** Rearrangements do not
  shorten the alignable portion of a read or cause mates to overlap.
- **Magnitude.** 8.5-17.3% of read pairs implies on the order of a million or
  more junctions per genome. Chromothripsis involves hundreds to low thousands.

### 5.8 Candidate additional checks — NOT validated

- **Usable coverage floor.** MDA010_N01 at 16.8x is low for a normal regardless
  of the artifact. A floor around 25-30x for normals seems defensible but is
  untested here.
- **PCT_PF_READS_IMPROPER_PAIRS.** Tracks the primary signal closely (0.37-0.54%
  unaffected vs 5.68-12.09% affected) but adds nothing over 5.1.
- **PCT_ADAPTER.** Near zero in most affected samples, so it does *not*
  discriminate — worth recording because it is the metric one would naively
  expect to fire.
- **Tumor/normal sample swap or mismatch.** Valuable and a common real failure,
  but requires genotype concordance rather than alignment metrics. Conpair is
  already in the container set.

---

## 6. Why the existing QC did not catch this

Every number in section 2 was computed by our own pipeline and written to disk
before the run failed. Nobody looked. That is a reporting failure, not a
measurement failure.

Faults to avoid in the replacement:

1. **No reference range.** `PCT_CHIMERAS 0.120206` is rendered identically to
   every other field. Nothing marks 12% as catastrophic and 0.1% as normal.
2. **No verdict.** Nothing states whether a sample is usable.
3. **Constants presented as findings.** Qualimap opens with reference genome size
   and contig count — identical for every sample ever run against the same
   reference. **A metric that cannot vary cannot be QC.** Leading with one trains
   the reader to skim, and skimming is what buried this.
4. **Flat hierarchy.** Hundreds of numbers, uniform styling, no severity order.
5. **Per-sample only, and incomplete.** No cohort view. The most diagnostic fact
   available — which samples are affected and which are not — required assembling
   16 files by hand. And qualimap silently covered only 12 of them (section 0).
6. **Jargon without translation.** "Chimeric" and "supplementary alignment" mean
   nothing to the audience that most needs to act.

### 6.1 Design rules for the replacement

1. Verdict first, in words, before any number.
2. Only show numbers that can be wrong.
3. Every number carries its reference range in the same row.
4. Express severity as a multiple — "90x higher than normal" lands where "12.0%"
   does not.
5. Cohort view first, per-sample detail second.
6. Sort worst first, never alphabetically.
7. Passing samples collapse to one line.
8. State the required action, not just the metric.
9. **Assert cohort completeness.** Print "16 of 16 samples checked." A report that
   silently covers a subset is worse than no report.

---

## 7. Specification for the pre-QC tool

### 7.1 Two placements, two mechanisms

**Primary — after Map, before variant calling.** Reads the Picard metrics the
pipeline already produces. Zero new computation, authoritative numbers, catches
everything in section 5. This is where the gate belongs, and it still sits ahead
of all the expensive work.

**Secondary — before Map, for externally delivered BAMs.** Where BAMs arrive
pre-aligned and we want a verdict before committing even to Map, Picard metrics
do not yet exist. A sampling check using indexed region queries gives an adequate
approximation: one 1 Mb window on a 172 GB BAM returned 400,000 reads in **1.1
seconds**, yielding a supplementary rate within ~1 point of the whole-BAM figure.
Use four or more windows on different chromosomes.

Given the 45x separation, the fast approximation is sufficient to block a run.
Prefer the Picard path wherever it is available.

### 7.2 Inputs

BAM mapping TSV **and the pairing file** — the latter is required for 5.5 and was
missing from the earlier design.

### 7.3 Metrics

| Metric | Level | Source | Gate |
|---|---|---|---|
| PCT_CHIMERAS | sample | `*.asm.txt` PAIR | FAIL > 1.0% |
| PCT_SOFTCLIP | sample | `*.asm.txt` PAIR | FAIL > 1.0% |
| Aligned read fraction | sample | `*.asm.txt` PAIR | FAIL < 0.95 |
| PCT_EXC_OVERLAP | sample | `*.wgs.txt` | WARN > 10% |
| PCT_EXC_TOTAL | sample | `*.wgs.txt` | WARN > 35% |
| MEAN_COVERAGE | sample | `*.wgs.txt` | record; floor TBD (5.8) |
| Delivered/usable ratio | sample | delivery figure + `*.wgs.txt` | WARN > 1.5x (5.5) |
| Insert concordance | **pair** | needs `CollectInsertSizeMetrics` | FAIL > ~1.5x (5.6) |
| Cohort completeness | run | count vs mapping file | FAIL on any missing sample |

### 7.4 Statistical note

**Report median insert size, never the mean.** The qualimap means for affected
samples are 350,000-664,000 bp where the medians are 132-212 bp. A few percent of
chimeric pairs mapping megabases apart destroys the mean while the median stays
robust. This generalises to any heavy-tailed metric, and matters for
outward-facing reports, where "average fragment length 494,733" would be read as
a real measurement.

### 7.5 Outputs

1. **Exit code** — for the run script.
2. **Cohort summary** — one line per sample, worst first, plus pair verdicts.
3. **Per-failure explanation card** — plain language, for the sequencing core.
4. **Machine-readable TSV** — for trend tracking. Six months of accumulated
   `PCT_CHIMERAS` would turn the section 5 thresholds from n=27 estimates into
   properly grounded ranges.

---

## 8. Output template — cohort summary

```
==================================================================
  BAM PRE-FLIGHT QC
  Project 16840_N  |  16 of 16 samples checked  |  2026-07-14
==================================================================

  RESULT:  10 of 16 samples FAILED.
           8 of 8 pairs unusable.  Pipeline run BLOCKED.

  All 8 normals failed, plus the tumors of MDA001 and MDA002.
  Because every normal is affected, no pair can be analysed --
  including the 6 pairs whose tumors are clean.

------------------------------------------------------------------
  FAILED SAMPLES                                       worst first
------------------------------------------------------------------
  SAMPLE              T/N   SPLIT READS      READ USED   COVERAGE
                            (normal: 0.1%)   (norm 100%) (usable)
  APTL_MDA009_N01      N     17.3%   173x       86%        32x
  APTL_MDA002_N01      N     17.1%   171x       68%        20x
  APTL_MDA001_N01      N     14.8%   148x       83%        28x
  APTL_MDA011_N01      N     13.8%   138x       89%        31x
  APTL_MDA012_N02      N     13.1%   131x       87%        29x
  APTL_MDA010_N01      N     13.0%   130x       83%        17x  LOW
  APTL_MDA001_T01      T     12.9%   129x       77%        25x
  APTL_MDA008_N01      N     12.0%   120x       79%        25x
  APTL_MDA002_T01      T      9.8%    98x       70%        23x
  APTL_MDA007_N01      N      8.5%    85x       92%        39x

------------------------------------------------------------------
  PASSED SAMPLES
------------------------------------------------------------------
  APTL_MDA007_T01_D2, APTL_MDA008_T01, APTL_MDA009_T01,
  APTL_MDA010_T01_D, APTL_MDA011_T01, APTL_MDA012_T01
  All within normal ranges (105-114x usable coverage).

------------------------------------------------------------------
  PAIR CHECKS
------------------------------------------------------------------
  All 8 pairs FAIL: normal sample failed sample-level QC.
  Additionally, insert size divergence ~2x in all pairs --
  copy number calling unreliable, consider unmatched normal.

==================================================================
```

Note what is absent: reference genome size, contig counts, GC content, and every
other number that could not have been wrong.

---

## 9. Output template — per-failure card

For the sequencing core. No technical terms without definition.

```
------------------------------------------------------------------
SAMPLE: APTL_MDA008_N01                              RESULT: FAILED
------------------------------------------------------------------

THE PROBLEM
  A large fraction of the sequencing reads in this sample cannot
  be placed on the genome as whole, intact reads.

WHAT WE MEASURED
                              This sample     Normal range
  Reads split across two
    genomic locations           12.0%          under 0.2%
  Portion of each read that
    could be used                79%           over 99%
  Usable coverage after
    quality filtering            25x           45x or more
  DNA fragment length           132 bases      350-400 bases

WHAT THIS MEANS
  Each read is 151 bases long. In a normal sample essentially all
  151 align to the genome. In this sample only about 119 do; the
  rest are discarded.

  The fragment length measurement points to why: at 132 bases, the
  fragments are shorter than a single 151-base read.

WHY IT MATTERS TO US
  Our software identifies rearranged chromosomes by finding reads
  split across two locations. That is exactly what these reads look
  like. The software was overwhelmed and shut down every time.

  Even had it finished, the results would be unusable: false
  signals from this sample would mask real findings in the matched
  tumor.

WHAT WE CANNOT TELL FROM THE DATA
  We can measure what the sequence data looks like, but not what
  produced it. We are not drawing conclusions about sample
  handling or preparation.

WHAT WE NEED
  This sample cannot be used for structural variant analysis as
  delivered, and re-processing the existing data will not change
  that. We would like to discuss whether new data can be generated.

FOR CONTEXT
  16 samples were submitted. 10 show this pattern, including all 8
  normals and 2 of the 8 tumors. 6 tumors are entirely normal, so
  the issue is not with the sequencing run as a whole.
------------------------------------------------------------------
```

### 9.1 Why this framing

- **"Split across two genomic locations" replaces "chimeric."** The technical term
  appears nowhere.
- **"Portion of each read that could be used"** is the most intuitive of the three
  primary metrics and needs no statistical background.
- **The 6 clean tumors are cited** to foreclose "your pipeline is broken." Same
  run, same pipeline, 16 samples, 6 entirely normal.
- **No causal claim.** The earlier draft speculated about input DNA quantity and
  degradation. Removed — we have no visibility there.

---

## 10. Known gaps

### 10.1 Read length — RESOLVED

`MEAN_READ_LENGTH` is 151 across both projects. Previously assumed and flagged.

### 10.2 The sampling fallback is timed but not validated

One 1 Mb window on one BAM at 1.1 seconds with good agreement. Multi-window
aggregation, healthy-control behaviour, and window-choice sensitivity are all
untested. Validation is read-only and takes minutes. **This matters less than it
did**, since the Picard path (7.1) is now primary and needs no validation at all.

### 10.3 Insert size is not in the current Picard outputs

`CollectInsertSizeMetrics` is not run. The pair-level gate (5.5) currently depends
on qualimap. Adding it to the Map stage would put the whole gate on one source.

### 10.4 Thresholds rest on n=27 across two projects

The 45x separation makes them forgiving, but they are two cohorts. Emit the TSV
(7.5) so ranges become empirical.

### 10.5 Cause is not established and is outside our visibility

We can characterise the data. We cannot say what produced it, and should not
speculate in any outward-facing document.

### 10.6 MDA001 and MDA002 were not in the SV run

The Delly run covered six pairs (MDA007-012). MDA001 and MDA002 have metrics but
were not part of that run. Their status in the project should be confirmed.

---

## 11. Immediate recommendations

1. **Do not re-run Proj_16840_N** until the normals are resolved. No resource
   configuration produces usable somatic SV calls.
2. **Correct the record on scope:** ten affected samples, not six, and two of them
   are tumors.
3. **Build the section 7 tool against the Picard path first.** It is nearly free
   and needs no validation.
4. **Add `CollectInsertSizeMetrics` to the Map stage** (10.3).
5. **Confirm the status of MDA001 and MDA002** (10.6).
6. **Flag APTL_MDA010_N01 separately** — 16.8x usable coverage is marginal for a
   normal irrespective of the artifact.
