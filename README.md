# wgsTriage 0.9.0

A post-mapping QC filter for WGS cohorts. It reads only what the Map stage
already wrote to disk, computes nothing new, and renders a verdict per sample
and per tumor/normal pair. A defective cohort is caught in about a minute
instead of after three days of compute.

Everything needed to catch the Proj_16840_N failure was already on disk before
the pipeline started. Nothing looked at it. That is the whole reason this
exists.

Triage is the right word deliberately. Compute is a scarce resource, and
deciding which cases can benefit from it comes before treatment rather than
after. It also stays neutral about cause: this tool characterises the data, it
does not say what produced it.

## Layout

```
bin/wgsTriage.R             assess one project, write the report
bin/wgsTriageBackground.R   build reference ranges from previously mapped projects
R/qcLib.R                   shared parsing and filter threshold logic
data/background/            generated reference ranges, aggregate only
docs/METHODS.md             how the thresholds were set, and why
docs/DECISIONS.md           what was chosen, and what is still unvalidated
docs/NORMAL_BAM_QC_REPORT.md  the original analysis the thresholds derive from
TODO_260719.md              outstanding work, kept at the root to stay visible
tests/testThresholds.R      smoke test over synthetic fixtures
```

## Usage

Both scripts answer `--help` with a full description, and report the version in
the first line of that output. They resolve their own location, so they can be
run from any working directory.

Assess a project:

```
Rscript bin/wgsTriage.R <MapDir> [--background <BgDir>] [--out <OutDir>] [--project <Name>]
```

Rebuild the background, needed only when new projects land:

```
Rscript bin/wgsTriageBackground.R <QCDir> [--out <OutDir>]
```

`<MapDir>` is a Map stage output directory. `<QCDir>` is an archive of
previously mapped projects, searched recursively. Reference ranges ship in
`data/background/`, so the first command works without running the second.

## Inputs read

```
<MapDir>/out/metrics/<sample>/<sample>.asm.txt    Picard CollectAlignmentSummaryMetrics
<MapDir>/out/metrics/<sample>/<sample>.wgs.txt    Picard CollectWgsMetrics
<MapDir>/sbam/multiqc/multiqc_data/multiqc_samtools_stats.txt
```

Read only. Nothing is modified and nothing is recomputed, which is what makes
the tool free to adopt: it cannot break a pipeline it only reads from.

Qualimap is deliberately not read. It overstates usable coverage by up to 2.6x
on degraded samples, and in Proj_16840_N it silently covered only 12 of 16
samples, producing a conclusion that had to be retracted.

If the multiqc file is absent the report still runs on Picard alone, says so
plainly, and names the thresholds it could not evaluate.

## Filter thresholds

| Metric | Source | Warn | Fail |
|---|---|---|---|
| `pctChimeras` | Picard asm, PAIR row | above 1.0% | above 5.0% |
| `supplementaryRate` | samtools | above 0.5% | above 1.0% |
| `pctSoftclip` | Picard asm, PAIR row | above 1.0% | above 5.0% |
| `pctReadUsed` | Picard asm, PAIR row | below 98% | below 95% |
| `pctProperlyPaired` | samtools | below 98.5% | below 97.0% |
| `pctExcOverlap` | Picard wgs | above 10% | none |
| `pctExcTotal` | Picard wgs | above 35% | none |

A sample tripping three or more warnings at once is failed even when no single
threshold fires, because these metrics are not independent and move together
under genuine degradation.

Coverage floors are 30x for normals and 80x for tumors. A sample whose class
cannot be read from its name is held to the tumor floor. Falling below the
floor raises a warning, which counts toward the three above; coverage alone
never fails a sample.

Pair checks are separate. A pair fails if either member fails, or if the tumor
and normal insert size distributions differ by more than 1.5x. Two samples can
each pass on their own and still be unusable together.

Verdicts are `FAIL`, `WARN`, `PASS` and `INCOMPLETE`. `INCOMPLETE` means a
metric could not be read, and it is never reported as `PASS`: a sample that was
not checked has not passed.

See `docs/METHODS.md` for how every number above was chosen.

## Outputs

| File | Purpose |
|---|---|
| `wgsTriage.txt` | Console transcript, verdict first, worst first |
| `wgsTriage.html` | Standalone report with reference ranges and per-sample cards |
| `wgsTriage_samples.tsv` | One row per sample, for trend tracking |
| `wgsTriage_pairs.tsv` | One row per inferred tumor/normal pair |

The HTML contains plain-language failure cards intended to be sent to the data
provider. They carry no technical vocabulary and make no claim about cause.

All four carry real sample names and are gitignored.

**The script always exits 0.** It is advisory: it reports loudly and leaves the
run/do-not-run decision to a person. To make it blocking, change the final
`quit(status = 0)` in `bin/wgsTriage.R` to key off `nFail`. That is one line.

## Validation

| Cohort | n | Result |
|---|---|---|
| Proj_16840_N | 16 | 10 fail: all 8 normals plus the MDA001 and MDA002 tumors |
| Proj_17495_I | 11 | 0 fail, the negative control |
| Proj_17608 | 96 | 0 fail |
| ReMap_260130 | 268 | 23 fail |

Proj_16840_N reproduces the original finding exactly. Proj_17495_I produces no
false positives. Re-run both after moving anything: paths are what break in a
move.

## A note on the background

The archive is not clean. 23 of 454 historical samples fall below these
thresholds, all in `ReMap_260130`, at 7 to 22% chimeric read pairs. That cohort
predates Proj_16840_N and whether it was analysed and released is an open
question, not a settled one.

Reference ranges are therefore computed only from samples that pass. Including
the defective ones would widen the range enough to admit the next bad cohort,
and robust statistics alone are not sufficient protection when the contaminated
fraction is this large and concentrated in one place.

Only aggregate statistics are committed. The per-sample background carries
names and stays out of the repository until it is anonymised.

See `docs/DECISIONS.md` for open questions and unvalidated assumptions.

## Version

0.9.0. Defined once as `WGSTRIAGE_VERSION` in `R/qcLib.R` and reported by both
scripts under `--help`. The version tracks the thresholds: a report can be
traced back to the threshold set that produced it, so bump it whenever a
threshold in `R/qcLib.R` moves.
