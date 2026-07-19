# PLAN_01 — first version of the repository

What to create, what to copy, what to leave behind. Target is a repo that is
complete enough to defend and small enough to read in one sitting.

---

## 1. Name

**`wgsTriage`**

Triage needs no explanation to a clinical audience, and it carries the argument
without having to make it: scarce resources go to the cases that can benefit,
and deciding that comes before treatment rather than after. That is precisely
the case for this tool. Compute is the scarce resource, three days is what
gets spent on a case that could never have benefited, and the assessment
already existed in files nobody read.

It also stays neutral about cause. Triage sorts by condition, not by blame,
which keeps the tool on ground it can defend — we can characterise the data, we
cannot say what produced it.

Contains `wgs`, and matches the camelCase convention already used across
`getWGSStats.R`, `nfTraceReport.R`, `qcReport01.R`.

**Note on verdict vocabulary.** The name offers one, and I am deliberately not
taking it. Verdicts stay `FAIL` / `WARN` / `PASS` / `INCOMPLETE` rather than
becoming triage categories or colours, because "FAILED" states the consequence
and "Red" requires a legend. The metaphor belongs in the name and in how you
explain the tool, not in the report a tired analyst reads at 7pm.

---

## 2. Renames from the current working copy

Fresh repo is the right moment to fix the names. `preflight` was always the
wrong word for the tool.

| Now | In the repo |
|---|---|
| `preflightQC.R` | `bin/wgsTriage.R` |
| `importQCBackground.R` | `bin/wgsTriageBackground.R` |
| `R/qcLib.R` | `R/qcLib.R` (unchanged) |

Only two internal references need updating: `source(path(scriptDir, "R", "qcLib.R"))`
in each script resolves relative to the script, so moving both into `bin/` means
changing that to `path(scriptDir, "..", "R", "qcLib.R")`. Nothing else in the
code refers to a filename.

---

## 3. File manifest

`COPY` = move the existing file. `NEW` = write from scratch. `GEN` = output of
running the importer.

```
wgsTriage/
├── README.md                              NEW
├── .gitignore                             NEW
├── bin/
│   ├── wgsTriage.R                        COPY  (preflightQC.R)
│   └── wgsTriageBackground.R              COPY  (importQCBackground.R)
├── R/
│   └── qcLib.R                            COPY
├── data/background/
│   ├── backgroundStats.tsv                GEN   1.8 KB, no names
│   ├── backgroundCoverageStats.tsv        GEN   346 B, no names
│   └── backgroundMetricCoverage.tsv       GEN   357 B, no names
├── TODO_260719.md                         COPY  root, not docs/ -- see below
├── docs/
│   ├── METHODS.md                         NEW
│   ├── DECISIONS.md                       COPY
│   └── NORMAL_BAM_QC_REPORT.md            COPY
└── tests/
    ├── testThresholds.R                   NEW
    └── fixtures/miniCohort/               NEW
        ├── out/metrics/CLEAN_N01/CLEAN_N01.{asm,wgs}.txt
        ├── out/metrics/DEFECT_N01/DEFECT_N01.{asm,wgs}.txt
        └── sbam/multiqc/multiqc_data/multiqc_samtools_stats.txt
```

Seventeen files. That is the whole repo.

`TODO_260719.md` stays at the repo root by preference, not by accident. A TODO
list filed under `docs/` stops being read. Do not tidy it away.

### Why the background can ship now

This was the one thing I expected to block v1, and it does not. The three files
under `data/background/` are pure aggregates — columns are
`metric, n, median, mad, q01…q99, min, max`. **No sample names, no project
names, 2.5 KB total.** They are safe to commit today.

The two files that carry names, `backgroundSamples.tsv` (455 rows) and
`backgroundFlagged.tsv` (24 rows), stay out until the anonymisation work in
`TODO_260719.md` section 3 is done. Nothing in `wgsTriage.R` reads them, so the
gate is fully functional without them: reference ranges, fold-change severity,
and the "reference ranges from N previously mapped samples" line all come from
`backgroundStats.tsv`.

Consequence worth stating plainly: **v1 ships a working gate with full
historical context and no PHI in the repository.** The anonymisation is needed
to publish the per-sample evidence, not to run the tool.

---

## 4. Contents of each NEW file

### `README.md`

Purpose, usage, outputs, gates table, validation summary. Most of it exists in
the current `README.md` — carry it over and update paths for `bin/`. Must state
in the first paragraph that the tool reads only what Map already wrote and
computes nothing new, because that is what makes it free to adopt.

### `.gitignore`

Guards against the accidental commit that undoes the PHI decision. As written,
this differs from the first draft of the plan in two ways, both deliberate.

The background rule denies by default and allows only the three name-free
aggregates, rather than naming the two files to exclude. Listing the exclusions
holds only as long as those stay the only name-bearing outputs; inverting it
means anything the importer starts emitting later is ignored until someone
consciously allows it. The failure mode becomes a safe file left out, which is
visible, instead of a name-bearing file committed, which is not.

The reports are excluded by filename as well as by directory, because they
carry sample names and land wherever `--out` points, which is not always a
`preflight*` directory.

```
# Background. Deny by default, allow only the three name-free aggregates.
data/background/*
!data/background/backgroundStats.tsv
!data/background/backgroundCoverageStats.tsv
!data/background/backgroundMetricCoverage.tsv

# Any sample-name crosswalk, until anonymised
*crosswalk*

# Generated reports, by directory and by filename
preflight*/
wgsTriage_out/
preflightQC.txt
preflightQC.html
preflightQC_samples.tsv
preflightQC_pairs.tsv

# Input archive, large and name-bearing
QCData/

# R
.Rhistory
.RData
.Rproj.user/
```

### `docs/METHODS.md`

The document that makes the thresholds defensible. This is the one to have open
in the conversation with the boss. Sections:

1. **What is read** — the three input files, and why qualimap is excluded
   (overstates usable coverage 2.6x, silently covered 12 of 16 samples).
2. **The seven gates** — metric, source, warn, fail, and one line on what each
   physically measures.
3. **How the thresholds were set** — the 421-sample band table, why FAIL moved
   from 1.0% to 5.0%, and the observed clean range for each metric. Assembled
   from `DECISIONS.md` section 1 and `NORMAL_BAM_QC_REPORT.md` section 5.
4. **Corroboration** — Picard and samtools are independent tools measuring the
   same defect. This is the paragraph that answers "how do you know?"
5. **Escalation and pair rules** — the 3-warning promotion, the insert-size pair
   gate, why a pair can fail when both members pass.
6. **Validation** — the four-cohort table below, with the command to reproduce it.
7. **What this cannot see** — sample swaps, contamination, anything needing
   genotypes. Stating the limits is what makes the rest credible.

Validation table for section 6:

| Cohort | n | Result |
|---|---|---|
| Proj_16840_N | 16 | 10 fail: all 8 normals plus the MDA001 and MDA002 tumors |
| Proj_17495_I | 11 | 0 fail (negative control) |
| Proj_17608 | 96 | 0 fail |
| ReMap_260130 | 268 | 23 fail |

### `tests/testThresholds.R` and `tests/fixtures/miniCohort/`

One smoke test. A gate that blocks multi-day compute should not be able to
silently stop gating, and the thresholds are the crown jewels — a regression
there is invisible until it costs another three days.

Fixtures are two synthetic samples built by trimming real Picard files: keep the
header, the `## METRICS CLASS` line, the column header and the `PAIR` row, drop
the histogram block. Rename the samples to `CLEAN_N01` and `DEFECT_N01`, which
also keeps real names out of the repo. About 4 KB total.

The test asserts:

- `DEFECT_N01` returns `FAIL`, `CLEAN_N01` returns `PASS`
- cohort completeness reports 2 of 2
- exit status is 0
- the run works with `--background` pointing at `data/background/`

That is sufficient. Do not build a test suite; build this one test.

---

## 5. What does NOT go in the repo

| Item | Size | Why not |
|---|---|---|
| `QCData/` | 1.4 GB | Input archive. Stays on disk; the importer takes a path. |
| `data/background/backgroundSamples.tsv` | 455 rows | Carries `project` + `sample`, 60 of them DMP or accession style. Blocked on anonymisation. |
| `data/background/backgroundFlagged.tsv` | 24 rows | Same. Also the evidence list for the ReMap_260130 problem — keep it, locally. |
| `preflight_16840_N/`, `preflight_17495_I/` | 88 KB | Generated output. Regenerate on demand. |
| `DELLY_FAILURE_REPORT.*`, `delly_nf_traces.txt`, `execution_trace__DELLY.txt` | 340 KB | Incident forensics for one project, superseded. Belongs with the project, not the tool. |
| `NORMAL_BAM_QC_REPORT_PI.md`, `NORMAL_BAM_QC_EMAIL_SUMMARY.md` | 15 KB | Audience-specific communications, not tool documentation. |
| `handoff-adagio-preflight-qc-report.md` | 8.7 KB | Superseded by this repo existing. |

**Resolved, 2026-07-19.** `NORMAL_BAM_QC_REPORT.md` carries sample IDs
(`APTL_MDA008_N01` and similar). These were confirmed to be properly anonymised
study codes that will appear in the final paper, not accession numbers like
`C-0KNW20` or `S16-8625`. They are safe to commit and safe to publish. Do not
scrub them. The file goes into `docs/` because the thresholds are indefensible
without their derivation.

This does not relax anything else: `backgroundSamples.tsv`,
`backgroundFlagged.tsv` and the generated `preflightQC*` reports still carry
uncontrolled names and stay out of the repo, per `.gitignore`.

---

## 6. Order of operations

1. Create the repo, `git init`, branch `feat/init-01` — not master.
2. Write `.gitignore` **first**, before anything generated exists and before the
   first `git add`. This is the step that prevents the mistake.
3. Copy the four `COPY` files into the layout above; apply the two renames and
   fix the `source()` path in both scripts.
4. Regenerate the background into `data/background/` and confirm only the three
   aggregate files are tracked.
5. Write `README.md` and `docs/METHODS.md`.
6. Build the fixtures and `tests/testThresholds.R`; confirm it passes.
7. Run the gate against Proj_16840_N and Proj_17495_I from the new location to
   confirm nothing broke in the move — 10 fail and 0 fail respectively.
8. Commit. Draft the message in `/tmp` first, conventional commits, no emoji.
9. Create `master` from the finished branch, then merge. See below.
10. Add the GitHub remote and push `master`. See below.

Step 7 is not optional. The only thing that changed is paths, and paths are
exactly what breaks in a move.

`.gitignore` is step 2 rather than step 4 by correction: in the original order
the background was generated before the ignore rules existed, which opens
exactly the window the rule is meant to close.

### 6.1 Where master comes from

`git init` does not create a branch. It points `HEAD` at an unborn
`init.defaultBranch` (here `master`), and that branch only materialises on its
first commit. Because step 1 switches to `feat/init-01` before committing,
**`master` is never created.** After the initial commit the repo has exactly one
branch. This is expected, not a mistake, but it has to be closed out explicitly
or there is no trunk to merge into or push.

Create `master` once the branch is actually presentable — after step 7, when the
tests pass and validation reproduces. That way master's first commit is a
coherent v1 rather than a half-built repo:

```
git branch master          # master now points at the current commit
git checkout master
git merge --ff-only feat/init-01
```

From then on the normal loop applies, per `~/.claude/CLAUDE.md`: master is the
trunk, work happens on short-named branches (`feat/updates-01`), and those merge
back to master. Never commit directly to master.

### 6.2 Going to GitHub

Push `master`, and let GitHub take its name from what you push:

```
git remote add origin git@github.com:<org>/wgsTriage.git
git push -u origin master
```

One trap worth avoiding. GitHub has defaulted new repositories to `main` since
2020, while this repo and `init.defaultBranch` both say `master`. If the GitHub
repo is created through the web UI with a README or any initialising file, it
arrives with a `main` branch and pushing `master` leaves two unrelated trunks in
one repository. Create the GitHub repo **empty** — no README, no `.gitignore`,
no licence — and push into it. `gh repo create <org>/wgsTriage --private
--source=. --push` does this correctly in one step.

Before the first push, confirm the PHI decision survived. `.gitignore` protects
the working tree, not history, and a push is not reversible in the way a local
commit is:

```
git ls-files | grep -E 'backgroundSamples|backgroundFlagged|preflightQC|crosswalk'
```

That must return nothing.

---

## 7. Deliberately deferred

Not in v1, and each is a line in `TODO_260719.md` rather than a gap in the plan:

- Anonymised per-sample background (unblocks publishing the evidence, not the tool)
- Manifest support, `--manifest` and `--make-manifest`
- Pipeline integration into `runTempoWGSBam.sh` or Nextflow
- Recomputing `.asm.txt` under current Picard to give `pctReadUsed` a background
- Example reports checked in as documentation — they carry real sample names, so
  they wait for the anonymisation work

---

## 8. One thing to say out loud in that conversation

The tool finds 23 defective samples in `ReMap_260130`, a cohort that predates
Proj_16840_N. That is not a comfortable fact, but it is the argument. The
measurement was always available in files the pipeline already wrote; what was
missing was anything that looked at them and said no. This repo is that thing,
it validates cleanly against a known-good cohort, and it would have caught both
incidents before compute was committed.
