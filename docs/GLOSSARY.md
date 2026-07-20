# Glossary

Every column and term emitted by `wgsTriage` 0.9.0, with its source.

Threshold values and their derivation are in `docs/METHODS.md`.
`tests/testGlossary.R` checks that every output column appears here.

---

## 1. The HTML report, `wgsTriage.html`

The report sent to the data provider. Sections in the order they appear.

| Section | Contents | Membership |
|---|---|---|
| Banner and counts | Cohort verdict and six headline counts. | Always. |
| Notes | Incomplete cohort, reduced check set. | Cohort is incomplete, or a check's entire source is absent. |
| `Cohort` | One row per sample: name, class, verdict, one cell per applied check, coverage. Each cell shows the measured value above the reference median. Ordered failures, incomplete, warnings, passes. | Every sample submitted. |
| `Tumor / normal pairs` | Patient, tumor, normal, both insert sizes, ratio, verdict, reason. | All inferred pairs. States when none could be inferred. |
| `Per-sample detail for the data provider` | One card per sample: the checks it failed, each value against the reference median, and the effect on downstream analysis. | Verdict `FAIL`. |
| `Thresholds applied` | Each check with its plain-language description, fail and warn values, reference median, and the number of background samples behind that median. | The checks evaluated for this cohort. |
| `Where each number comes from` | Report label, internal name, source tool and file, source field(s), transform. Generated from the checks applied plus `PROVENANCE` in `R/qcLib.R`. | The checks evaluated, plus `insertSizeAverage`, `insertRatio`, `sampleType` and `patient`. |

Column headings in the `Cohort` table are the report labels listed in section 4.

A check whose entire source is absent for a cohort is dropped rather than
reported as `MISSING` once per sample, and is named in the reduced-check-set
note.

Below 980px every table renders as one block per row, each cell labelled from
its `data-label` attribute. Tables do not scroll horizontally.

---

## 2. The console report, `wgsTriage.txt`

The same evaluation as section 1, as text.

| Section | Contents | Membership |
|---|---|---|
| Header | Project, samples checked of samples submitted, date. | Always. |
| `RESULT` | The cohort verdict, and pairs unusable of pairs found. | Always. |
| Incomplete-cohort warning | Count of samples with incomplete metrics. | Samples checked is below samples submitted. |
| Reduced-check-set note | Which checks were not evaluated. | A check's entire source is absent for the cohort. |
| `FAILED SAMPLES` | Sample, class, split-read rate with fold, read used, coverage, count of thresholds failed. Ordered by chimera rate, highest first. | Verdict `FAIL`. |
| `BELOW COVERAGE FLOOR` | Sample, class, coverage, the floor it fell below. Ordered by coverage, lowest first. | `lowCoverage` is `TRUE`, whatever the verdict. |
| `SAMPLES NEEDING REVIEW` | Sample, class, `verdictReason`. | Verdict `WARN`. |
| `SAMPLES NOT FULLY ASSESSED` | Sample, class, `verdictReason`. | Verdict `INCOMPLETE`. |
| `PASSED SAMPLES` | Names, wrapped, and the coverage range across them. | Verdict `PASS`. |
| `PAIR CHECKS` | Patient, `pairVerdict`, `pairReason`; then unpaired samples; then the insert size divergence block. | All inferred pairs, failures first. The divergence block lists pairs with `insertRatio` above 1.5. |
| Footer | Number of samples behind the reference ranges, or a statement that no background was loaded. | Always. |

The tumor is named alongside the patient where one patient contributes more
than one pair.

---

## 3. Source files

```
<MapDir>/out/metrics/<sample>/<sample>.asm.txt    Picard CollectAlignmentSummaryMetrics
<MapDir>/out/metrics/<sample>/<sample>.wgs.txt    Picard CollectWgsMetrics
<MapDir>/sbam/multiqc/multiqc_data/multiqc_samtools_stats.txt
```

Referred to below as **asm**, **wgs** and **samtools**. The asm file is read
from its `CATEGORY=PAIR` row, the samtools table from its `.recal` rows.

`data/background/backgroundStats.tsv` supplies the reference medians. Section 7.

---

## 4. `wgsTriage_samples.tsv`

One row per sample directory under `out/metrics`. A sample whose metrics could
not be read appears with empty values rather than being omitted.

### Identity

| Column | Definition |
|---|---|
| `project` | The `--project` argument, or the directory containing `<MapDir>`. |
| `sample` | File basename with the `.asm.txt` / `.wgs.txt` extension and any `.smap`, `.md` or `.recal` suffix removed. |
| `sampleType` | Sample class read from the sample name by regex: `N`, `T`, or `unknown` when no trailing class token is present. |
| `patient` | The sample name with its trailing class token removed. `APTL_MDA012_N02` and `APTL_MDA012_T01` both reduce to `APTL_MDA012`. Used to match samples into pairs. |

### Verdict

| Column | Definition |
|---|---|
| `verdict` | Sample verdict: `PASS`, `WARN`, `FAIL` or `INCOMPLETE`. Section 6. |
| `verdictReason` | One sentence naming the checks that produced `verdict`. |
| `nFail` | Number of checks with status `FAIL`. |
| `nWarn` | Number of checks with status `WARN`. |
| `nMissing` | Number of checks with status `MISSING`. |
| `failedMetrics` | Internal names of the failing checks, comma-separated. Empty when none failed. |
| `warnedMetrics` | Internal names of the warning checks, comma-separated. Empty when none warned. |

### Gated metrics

The seven fixed checks, in report order. Report label is the column heading in
the HTML `Cohort` table.

| Column | Report label | Source | Field(s) | Transform |
|---|---|---|---|---|
| `pctChimeras` | Split pairs | asm | `PCT_CHIMERAS` | x 100 |
| `supplementaryRate` | Multi-mapped | samtools | `supplementary_alignments`, `raw_total_sequences` | **derived**: `supplementary_alignments / raw_total_sequences` x 100 |
| `pctSoftclip` | Clipped bases | asm | `PCT_SOFTCLIP` | x 100 |
| `pctReadUsed` | Read used | asm | `MEAN_ALIGNED_READ_LENGTH`, `MEAN_READ_LENGTH` | **derived**: `MEAN_ALIGNED_READ_LENGTH / MEAN_READ_LENGTH` x 100 |
| `pctProperlyPaired` | Proper pairs | samtools | `reads_properly_paired_percent` | none, already a percentage |
| `pctExcOverlap` | Overlap loss | wgs | `PCT_EXC_OVERLAP` | x 100 |
| `pctExcTotal` | Total loss | wgs | `PCT_EXC_TOTAL` | x 100 |

`pctReadUsed` is the percentage of each delivered read that survived alignment.
No Picard field reports it directly.

`pctReadUsed` has no reference median. Every archived asm file was written by a
Picard version predating `MEAN_ALIGNED_READ_LENGTH`, so the metric is judged
against its fixed threshold alone.

### Coverage

| Column | Definition |
|---|---|
| `meanCoverage` | wgs, `MEAN_COVERAGE`. Depth after Picard's quality filtering, not raw depth. |
| `medianCoverage` | wgs, `MEDIAN_COVERAGE`. Not gated. |
| `coverageFloor` | The floor for this sample's class: 30x for `N`, 80x for `T`, 80x for `unknown`. |
| `lowCoverage` | `TRUE` when `meanCoverage` is below `coverageFloor`. |

The coverage floor is the eighth check. It is defined in `bin/wgsTriage.R`
rather than in `THRESHOLDS`, since its value varies by sample class. It warns
and never fails. A warning here counts toward `WARN_ESCALATION`.

### Ungated metrics

| Column | Definition |
|---|---|
| `insertSizeAverage` | samtools, `insert_size_average`. A mean over properly-paired reads, not a median. Gated at the pair level only. |
| `pctImproperPairs` | asm, `PCT_PF_READS_IMPROPER_PAIRS` x 100. |
| `totalReads` | asm, `TOTAL_READS`. |
| `meanReadLength` | asm, `MEAN_READ_LENGTH`. Read length as delivered. |
| `meanAlignedLength` | asm, `MEAN_ALIGNED_READ_LENGTH`. Read length aligned. Absent from older Picard output, in which case it and `pctReadUsed` are empty. |
| `alignedFrac` | `meanAlignedLength / meanReadLength`. The same quantity as `pctReadUsed`, unscaled: `pctReadUsed` is `alignedFrac` x 100. |
| `chimeraFold` | `pctChimeras` divided by the reference median for `pctChimeras`. Empty when no background is loaded. |
| `suppFold` | `supplementaryRate` divided by the reference median for `supplementaryRate`. Empty when no background is loaded or the cohort has no samtools data. |

Both fold columns express a value as a multiple of the reference median. For a
metric whose direction is `high` the fold is value over median; for one whose
direction is `low` it is median over value. A fold above 1 is worse than the
reference. Both metrics above are `high`.

---

## 5. `wgsTriage_pairs.tsv`

One row per tumor/normal pair, matched on `patient` between a sample classified
`T` and one classified `N`. A patient with several tumors contributes several
rows. Samples classified `unknown` form no pair.

| Column | Definition |
|---|---|
| `patient` | The patient stem the two samples were matched on. |
| `sampleTumor` | Name of the tumor sample. |
| `sampleNormal` | Name of the normal sample. |
| `sampleTypeTumor` | Class of the tumor sample. `T` by construction. |
| `sampleTypeNormal` | Class of the normal sample. `N` by construction. |
| `insertSizeAverageTumor` | The tumor's `insertSizeAverage`. |
| `insertSizeAverageNormal` | The normal's `insertSizeAverage`. |
| `insertRatio` | Larger of the two insert size averages divided by the smaller. Always at least 1. |
| `insertStatus` | `FAIL` when `insertRatio` exceeds `INSERT_RATIO_FAIL`, 1.5. `MISSING` when either insert size is absent. Otherwise `PASS`. |
| `verdictTumor` | The tumor sample's `verdict`, unchanged. |
| `verdictNormal` | The normal sample's `verdict`, unchanged. |
| `pairVerdict` | `FAIL` when either sample failed or `insertStatus` is `FAIL`. `INCOMPLETE` when either sample is `INCOMPLETE`. Otherwise `PASS`. |
| `pairReason` | The condition that produced `pairVerdict`. Empty when the pair passes. |

A sample-level failure takes precedence in `pairReason` over insert divergence.
Both reports list divergent pairs in a separate block.

---

## 6. Verdict terms

Assigned per check, then rolled up per sample.

| Term | Definition |
|---|---|
| `PASS` | Per check: the value was read and is within range. Per sample: every check passed. |
| `WARN` | Per check: the value crossed the warn threshold but not the fail threshold. Per sample: at least one warning, no failures, none missing. |
| `FAIL` | Per check: the value crossed the fail threshold. Per sample: at least one failing check, or `WARN_ESCALATION` below. |
| `MISSING` | Per check only. The value could not be read: the field, the file or the sample's metrics were absent. |
| `INCOMPLETE` | Per sample only. No check failed and at least one was `MISSING`. |

`MISSING` does not roll up to `PASS`.

`WARN_ESCALATION` is the promotion of a sample with no failing check and three
or more warnings to `FAIL`. The count is `WARN_ESCALATION` in `R/qcLib.R`; the
reasoning is in `docs/METHODS.md` section 5.

Cohort verdict, shown in the HTML banner and the console `RESULT` line:

| Term | Condition |
|---|---|
| `BLOCK` | At least one sample has verdict `FAIL`. |
| `INCOMPLETE` | No sample failed, at least one is `INCOMPLETE`. |
| `REVIEW` | No sample failed or is incomplete, at least one has verdict `WARN`. |
| `CLEAR` | Every sample passed. |

The tool exits 0 in all four cases.

---

## 7. The background

Built by `bin/wgsTriageBackground.R` over previously mapped projects and
shipped as aggregates in `data/background/`.

**`backgroundStats.tsv`** is the file the report reads. One row per metric,
with columns `n`, `median`, `mad`, `q01`, `q05`, `q25`, `q75`, `q95`, `q99`,
`min` and `max`. The report uses `median` as the reference shown under every
value in the cohort table and as the denominator of `chimeraFold` and
`suppFold`, and `n` as the number of samples behind that median.

**Which samples are included.** A reference sample is one that could be
evaluated on at least one check and failed none of the checks it could be
evaluated on. Samples with a failing check are excluded. A sample is not
required to carry both Picard files: one with only coverage metrics contributes
to the coverage ranges.

**Keying.** Background rows are keyed on `project` plus `sample`. `sample`
alone is not unique across the archive; 100 names occur in more than one
project.

**`backgroundCoverageStats.tsv`** holds the same statistics for `meanCoverage`
split by `sampleType`. It is the only metric given a per-class reference.

When no background is loaded, the checks run on their fixed thresholds, no
reference appears beside any value, `chimeraFold` and `suppFold` are empty, and
both reports state that no background was loaded.
