# Decisions and open questions

First pass, written 2026-07-18. Everything below is a judgement call I made to
keep moving. The ones under "Needs your call" are the ones I would not want left
to me.

---

## Needs your call

### 1. The chimera FAIL threshold moved from 1.0% to 5.0%

This is the biggest change from `NORMAL_BAM_QC_REPORT.md` section 5.1 and the
one most worth arguing about.

Section 5.1 proposed FAIL above 1.0%, on the basis that affected and unaffected
samples separated by 45x with nothing between them. That is true of the 27
samples it examined. It is not true of the 421 sample archive:

| Band | Samples |
|---|---|
| below 0.5% | 283 |
| 0.5 to 1% | 107 |
| 1 to 5% | 9 |
| above 5% | 22 |

There is a real, sparse middle. A 1.0% gate fails 31 of 421 historical samples
(7.4%); a 5.0% gate fails 22 (5.2%) and leaves 9 as warnings.

I moved FAIL to 5.0% and made 1.0% a WARN. Every sample in the Proj_16840_N
disaster ran 8.5 to 17.3% and still fails outright, so nothing that mattered is
lost. My reasoning is that a gate firing on 7% of all historical work gets
switched off, and a switched-off gate is exactly the failure this exists to
prevent.

**If you would rather be noisier than permissive, change `fail` for
`pctChimeras` and `pctSoftclip` back to 1.0 in `R/qcLib.R`. It is one number
each.** The nine borderline samples are listed in
`background/backgroundFlagged.tsv`.

### 2. Twenty-three historical samples fail, all in `ReMap_260130`

Chimera rates of 7 to 22%, i.e. the same defect class as Proj_16840_N, in a
cohort that predates it. Names include `Umich3_N`, `Umich4_N`, `Umich11_N`,
`22-4426_T`, `23-71689_T`, `16-8625_T`, `CTCL91_66_T`.

I do not know whether that cohort was analysed and released. If it was, results
derived from those samples are suspect on the same grounds as Proj_16840_N.
**This is the finding I would look at first tomorrow.** Full list in
`background/backgroundFlagged.tsv`.

### 3. Warning escalation threshold is a guess

Three simultaneous warnings promotes a sample to FAIL. The rationale is sound
(these metrics are not independent and move together under real degradation) but
the number 3 is arbitrary. It currently promotes exactly one archive sample
(`Umich10_T` in `ReMap_260130`, 3.78% chimeras with five warnings). Change
`WARN_ESCALATION` in `R/qcLib.R` if you disagree.

### 4. Advisory, not blocking

Per your answer, the script always exits 0. Given how strongly the brief was
worded about stopping bad runs, flagging this so the choice is deliberate. To
make it blocking, change the final `quit(status = 0)` in `bin/wgsTriage.R` to key
off `nFail`. One line.

---

## Decisions I made without asking

### Insert size comes from samtools, not Picard

Section 10.3 lists "insert size is not in the current Picard outputs" as a known
gap and recommends adding `CollectInsertSizeMetrics` to the Map stage. **That is
not necessary.** `multiqc_samtools_stats.txt` already carries
`insert_size_average` per sample, and it separates this cohort cleanly (defective
108 to 274 bases, clean 342 to 411). The pair gate works today with no change to
the Map stage.

Caveat: it is a **mean**, and section 7.4 is explicit that insert size should be
reported as a median because the metric is heavy-tailed. Samtools computes it
over properly-paired reads only, so it is far less distorted than the qualimap
means that prompted that warning, but it is not what section 7.4 asked for. The
full distribution is in the samtools stats `IS` block if you want a real median;
I did not parse it because the mean already discriminates.

### The samtools metrics are a second independent gate

Not in the original spec and worth keeping. `supplementary_alignments /
raw_total_sequences` measures essentially the same defect as `PCT_CHIMERAS` but
comes from a different tool on a different pass. It is also the tightest metric
available: 0.094 to 0.160% across every clean sample in the archive, so the 1.0%
gate sits 6x above the observed maximum, versus roughly 5x for chimeras.

`reads_properly_paired_percent` splits the cohort at 99% (defective 87.9 to 94.3,
clean 99.45 to 99.58) and is included too.

When Picard and samtools agree, the verdict is defensible to someone upstream who
does not want to hear it. That was worth the extra parsing.

### Tumor/normal inference is name-based and fails loudly

Regex on the `_N01` / `_T01` / `-N` / trailing `N` patterns. Works on every
sample in Proj_16840_N and Proj_17495_I. Returns `unknown` rather than guessing,
and unknown is displayed as `?`, never silently treated as normal. Unpaired
samples are named in the report.

This will mis-handle unusual naming. If it becomes a problem the fix is an
optional pairing-file argument, which I left out under YAGNI.

### Background keyed on project plus sample

100 sample names occur in more than one project in the archive, usually from
remapping under a new name. Keying verdicts on the bare sample name merged them
and assigned one project's failure to another project's sample. It produced two
wrong verdicts before I caught it. Now keyed on `project::sample`.

Worth knowing if you write anything else against `backgroundSamples.tsv`:
**`sample` is not a unique key in that file.**

### Reference ranges exclude failing samples

Robust statistics alone are not enough when the contaminated fraction is 5%+ and
concentrated in one cohort. Ranges come from the 346 samples that pass the gates
and have both Picard files.

### Qualimap not read at all

Per section 0 of the QC report. Nothing in this tool touches it.

---

## Known gaps

### No background for `pctReadUsed`

Every archived `.asm.txt` was produced by an older Picard that does not emit
`MEAN_ALIGNED_READ_LENGTH`. All 421 lack it; the current project has it. So the
read-used gate runs on its fixed threshold with no historical context, and the
report says so. This will resolve itself as new projects accumulate.

The metric is gated as a percentage (`pctReadUsed`, fail below 95%) rather than
as the fraction section 5.3 describes. Same quantity, but a bare `0.86` sitting
in a column of percentages gets read as `0.86%`, which inverts its meaning. The
raw `alignedFrac` is still written to the TSV.

### Samtools background is thin

Only 98 of 454 archived samples (21.6%) have multiqc data, and all 98 are clean,
so there is no observed defective distribution for those metrics. The thresholds
are anchored on the clean range plus a wide margin rather than on separation
between two observed populations.

### Coverage floors are chosen, not derived

Values settled 2026-07-20: 30x normal, 80x tumor, 80x where the class cannot be
read from the name. They now warn on the verdict and count toward the
three-warning escalation, replacing the 25x/50x advisory block I invented.

What remains a gap is the basis. These are the depths NS considers adequate,
not thresholds separating two observed populations the way the chimera
threshold does. The archive can say how often a floor fires but not whether the
figure is right, because it holds no record of which samples were later found
too shallow to analyse.

Adopting them moved no verdicts: seven samples in Proj_16840_N fall below floor
and all seven already failed, and no archived sample gains a third warning from
the check. They fire on 12.2% of the 401 archived samples carrying a coverage
figure, against 1.0% before. See `docs/METHODS.md` section 5.

### The 1.5x insert ratio is unconfirmed

Section 5.6 says this needs checking against Facets documentation or its authors.
Still true. It currently flags 4 of 8 pairs in Proj_16840_N.

### No sample-swap detection

Genotype concordance needs Conpair, which needs more than mapping metrics. Out of
scope given the "only mapping QC files" constraint, but it is a common real
failure and the one class of problem this gate cannot see.

### Not wired into the pipeline

`bin/wgsTriage.R` is standalone. Where it hooks into `runTempoWGSBam.sh` or as a
Nextflow process was listed as open decision 1 in the handoff and I have not
touched it. The script takes a Map directory and exits 0, so it can be dropped in
anywhere between the two stages.

---

## Odd things noticed in passing

- `Umich10_T` reports 661x mean coverage across three projects. Either it is not
  WGS or something is wrong with it. Not investigated.
- The archive contains `Proj_16840_C` and `UMich` copies of `Umich10_T` that have
  `.wgs.txt` but no `.asm.txt`, so they cannot be gated on the primary metric.
  They show as INCOMPLETE, which is the honest answer.
- `Proj_17608` and `ReMap_260130` share 96 sample names, consistent with a
  remap rather than distinct data.
