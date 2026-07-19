#!/usr/bin/env Rscript
##
## Post-mapping QC filter. Runs after Map, before anything expensive.
##
## Run with --help for usage. Reads only what the Map stage already produced
## and computes nothing new. Writes four report files, all name-bearing and
## all gitignored.
##
## Always exits 0. This is advisory by decision: it reports loudly and leaves
## the run/do-not-run call to a human. If that stops working, make it exit
## non-zero on FAIL, which is a one line change in the last block.
##

suppressPackageStartupMessages({
    library(tidyverse)
    library(fs)
    library(glue)
})

##
## Repo root, derived from this script's own location in bin/ and expressed
## relative to the working directory. Deliberately not here::here(): that
## resolves from the working directory, which for this tool is the data
## directory rather than the repo, and it fails silently when it guesses wrong.
##
scriptPath <- commandArgs(trailingOnly = FALSE) |>
    str_subset("^--file=") |>
    str_remove("^--file=")
repoRoot <- if (length(scriptPath) > 0) {
    path_rel(path_dir(path_dir(path_real(scriptPath))))
} else {
    "."
}

source(path(repoRoot, "R", "qcLib.R"))

##
## Arguments.
##
args <- commandArgs(trailingOnly = TRUE)
defaultBackground <- path(repoRoot, "data", "background")

getOpt <- function(flag, fallback) {
    i <- which(args == flag)
    if (length(i) > 0 && length(args) > i[1]) args[i[1] + 1] else fallback
}

usage <- function() {
    glue("
wgsTriage.R -- post-mapping QC filter. Reads only what the Map stage already
wrote and renders a verdict per sample and per tumor/normal pair.

Usage:
  Rscript bin/wgsTriage.R <MapDir> [--background <BgDir>] [--out <OutDir>] [--project <Name>]
  Rscript bin/wgsTriage.R --help

Arguments:
  <MapDir>               Map stage output to assess. Read only, nothing is
                         modified and nothing is recomputed. Expects to find:
                           <MapDir>/out/metrics/<sample>/<sample>.asm.txt
                           <MapDir>/out/metrics/<sample>/<sample>.wgs.txt
                           <MapDir>/sbam/multiqc/multiqc_data/multiqc_samtools_stats.txt

Options:
  --background <BgDir>   Reference ranges built by bin/wgsTriageBackground.R.
                         Default: <repoRoot>/data/background
                         If missing, the filters still run but on fixed
                         thresholds only, with no out-of-range detection.
  --out <OutDir>         Directory for the report. Default: ./preflight
  --project <Name>       Label shown in the report.
                         Default: the directory containing <MapDir>.
  -h, --help             Show this message and exit.

Writes into <OutDir>:
  preflightQC.txt          the console report, as text
  preflightQC.html         the same report, formatted
  preflightQC_samples.tsv  per-sample metrics and verdicts
  preflightQC_pairs.tsv    per tumor/normal pair verdicts

All four carry real sample names and are gitignored.

Always exits 0, including on FAIL. This is advisory by decision: it reports
loudly and leaves the run/do-not-run call to a human.
")
}

if (any(args %in% c("-h", "--help"))) {
    cat(usage(), "\n", sep = "")
    quit(save = "no", status = 0)
}

flags <- c("--background", "--out", "--project")
flagValues <- map_chr(flags, \(f) getOpt(f, NA_character_))
consumed <- c(flags, flagValues[!is.na(flagValues)])
positional <- args[!args %in% consumed]

if (length(positional) == 0) {
    stop(glue("
No Map directory given.

  Usage: Rscript bin/wgsTriage.R <MapDir> [--background <BgDir>] [--out <OutDir>] [--project <Name>]

  <MapDir> is the Map stage output to assess.
  Run with --help for the full description and the list of outputs.
"), call. = FALSE)
}

mapDir <- positional[1]
if (!dir_exists(mapDir)) {
    stop(glue("Map directory not found: {mapDir}\n  Run with --help for usage."),
         call. = FALSE)
}

backgroundDir <- getOpt("--background", defaultBackground)
outDir <- getOpt("--out", "preflight")
projectName <- getOpt("--project", path_file(path_real(path(mapDir, ".."))))
dir_create(outDir)

metricsDir <- path(mapDir, "out", "metrics")
multiqcFile <- path(mapDir, "sbam", "multiqc", "multiqc_data", "multiqc_samtools_stats.txt")

if (!dir_exists(metricsDir)) stop(glue("Metrics directory not found: {metricsDir}"))

##
## Cohort completeness.
##
## The denominator is every sample directory under out/metrics, not every
## sample we managed to parse. A sample whose mapping died leaves an empty
## directory, and counting only what parsed is precisely how the original
## qualimap analysis reported on 12 of 16 samples without noticing.
##
expectedFromDirs <- dir_ls(metricsDir, type = "directory", fail = FALSE) |> path_file()

picard <- collectPicardSamples(metricsDir)
if (is.null(picard) || nrow(picard) == 0) stop(glue("No Picard metrics parsed from {metricsDir}"))
picard <- picard |> mutate(project = projectName)

samtools <- readMultiqcSamtools(multiqcFile)
haveSamtools <- !is.null(samtools) && nrow(samtools) > 0

expectedSamples <- union(expectedFromDirs, picard$sample)
if (haveSamtools) expectedSamples <- union(expectedSamples, samtools$sample)
expectedSamples <- sort(expectedSamples)

dat <- tibble(sample = expectedSamples) |>
    left_join(picard, by = "sample") |>
    mutate(project = projectName,
           sampleType = classifySampleType(sample),
           patient = patientStem(sample))

if (haveSamtools) dat <- dat |> left_join(samtools, by = "sample")

dat <- dat |>
    mutate(hasAsm = !is.na(pctChimeras),
           hasWgs = !is.na(meanCoverage),
           hasSamtools = if (haveSamtools) !is.na(supplementaryRate) else FALSE)

##
## Drop filter thresholds whose entire source is absent for this cohort. A missing multiqc
## run is a cohort level fact worth stating once, not sixteen per-sample
## MISSING verdicts that bury the real signal.
##
usableThresholds <- THRESHOLDS |>
    filter(metric %in% names(dat)) |>
    filter(map_lgl(metric, \(m) any(!is.na(dat[[m]]))))

droppedThresholds <- THRESHOLDS |> filter(!metric %in% usableThresholds$metric)
droppedSources <- setdiff(THRESHOLDS$source, usableThresholds$source)

thresholdResults <- dat |>
    select(sample, all_of(usableThresholds$metric)) |>
    pivot_longer(-sample, names_to = "metric", values_to = "value") |>
    left_join(usableThresholds, by = "metric") |>
    mutate(status = case_when(
        is.na(value) ~ "MISSING",
        direction == "high" & !is.na(fail) & value > fail ~ "FAIL",
        direction == "low"  & !is.na(fail) & value < fail ~ "FAIL",
        direction == "high" & !is.na(warn) & value > warn ~ "WARN",
        direction == "low"  & !is.na(warn) & value < warn ~ "WARN",
        .default = "PASS"))

verdicts <- sampleVerdict(thresholdResults)

dat <- dat |> left_join(verdicts, by = "sample")

##
## Coverage floor. Advisory only and deliberately not a filter threshold: section 5.8 flags
## these figures as untested, and a normal at 17x is a judgement call about the
## analysis being attempted rather than a defect in the data. It is annotated
## everywhere the verdict appears so it cannot be missed.
##
dat <- dat |>
    mutate(coverageFloor = COVERAGE_WARN[sampleType],
           lowCoverage = !is.na(meanCoverage) & meanCoverage < coverageFloor,
           ## Short form keeps the fixed-width console columns aligned. "?" is
           ## deliberately not "N": an unclassifiable name is its own category.
           tn = recode(sampleType, unknown = "?"))

##
## Background reference ranges.
##
statsFile <- path(backgroundDir, "backgroundStats.tsv")
haveBackground <- file_exists(statsFile)

refStats <- if (haveBackground) {
    read_tsv(statsFile, show_col_types = FALSE, progress = FALSE)
} else {
    tibble(metric = character(), n = integer(), median = numeric(),
           q05 = numeric(), q95 = numeric())
}

refMedian <- set_names(refStats$median, refStats$metric)
refN <- set_names(refStats$n, refStats$metric)

## Severity as a multiple of the clean median. Design rule 4: "90x higher than
## normal" lands where "12.0%" does not.
foldOf <- function(metric, value) {
    ref <- refMedian[metric]
    if (is.na(ref) || is.na(value) || ref == 0) return(NA_real_)
    dir <- THRESHOLDS$direction[THRESHOLDS$metric == metric]
    if (length(dir) == 0) return(NA_real_)
    if (dir == "high") value / ref else ref / value
}

dat <- dat |>
    mutate(chimeraFold = map2_dbl("pctChimeras", pctChimeras, foldOf),
           suppFold = if (haveSamtools) map2_dbl("supplementaryRate", supplementaryRate, foldOf) else NA_real_)

##
## Pair level checks.
##
pairs <- if (haveSamtools) {
    evaluatePairs(dat, verdicts)
} else {
    dat |>
        select(sample, patient, sampleType) |>
        left_join(verdicts |> select(sample, verdict), by = "sample") |>
        filter(sampleType %in% c("N", "T")) |>
        (\(p) inner_join(filter(p, sampleType == "T"), filter(p, sampleType == "N"),
                         by = "patient", suffix = c("Tumor", "Normal")))() |>
        mutate(insertSizeAverageTumor = NA_real_, insertSizeAverageNormal = NA_real_,
               insertRatio = NA_real_, insertStatus = "MISSING",
               pairVerdict = case_when(
                   verdictTumor == "FAIL" | verdictNormal == "FAIL" ~ "FAIL",
                   .default = "PASS"),
               pairReason = case_when(
                   verdictTumor == "FAIL" & verdictNormal == "FAIL" ~ "both samples failed QC",
                   verdictNormal == "FAIL" ~ "normal failed QC",
                   verdictTumor == "FAIL" ~ "tumor failed QC",
                   .default = ""))
}

unpaired <- dat |>
    filter(!patient %in% pairs$patient) |>
    pull(sample)

##
## Counts.
##
nExpected <- length(expectedSamples)
nChecked <- sum(dat$hasAsm & dat$hasWgs)
nFail <- sum(dat$verdict == "FAIL")
nWarn <- sum(dat$verdict == "WARN")
nPass <- sum(dat$verdict == "PASS")
nIncomplete <- sum(dat$verdict == "INCOMPLETE")

failedSamples <- dat |> filter(verdict == "FAIL") |> arrange(desc(pctChimeras))
warnSamples <- dat |> filter(verdict == "WARN") |> arrange(desc(pctChimeras))
passSamples <- dat |> filter(verdict == "PASS") |> arrange(sample)
incompleteSamples <- dat |> filter(verdict == "INCOMPLETE") |> arrange(sample)

nPairFail <- sum(pairs$pairVerdict == "FAIL")
nPairs <- nrow(pairs)

overallVerdict <- case_when(
    nFail > 0 ~ "BLOCK",
    nIncomplete > 0 ~ "INCOMPLETE",
    nWarn > 0 ~ "REVIEW",
    .default = "CLEAR")

##
## Plain language explanation for each filter threshold, used in both outputs.
## Section 9.1: the technical term appears nowhere the sequencing core will read.
##
PLAIN <- c(
    pctChimeras = "Reads split across two separate genomic locations",
    supplementaryRate = "Reads reported at more than one location",
    pctSoftclip = "Portion of each read discarded at alignment",
    pctReadUsed = "Portion of each read that could be used",
    pctProperlyPaired = "Read pairs landing where they should",
    pctExcOverlap = "Bases dropped because mates overlapped",
    pctExcTotal = "Bases dropped for all reasons")

ACTION <- c(
    pctChimeras = "Structural variant calling will fail or produce false positives. Do not run SV analysis on this sample.",
    supplementaryRate = "Confirms the split-read finding from a second independent tool.",
    pctSoftclip = "Alignment is discarding a large fraction of each read. Usable data is far below the delivered figure.",
    pctReadUsed = "Most of each sequenced read is being thrown away. Effective yield is a fraction of what was billed.",
    pctProperlyPaired = "Paired-end structure is degraded. Insert size and copy number estimates will be unreliable.",
    pctExcOverlap = "Fragments are shorter than the read length, so mates overlap and are counted once.",
    pctExcTotal = "Usable coverage is far below raw coverage. Any externally reported depth is overstated.")

##
## Console report.
##
w <- 76
rule <- strrep("=", w)
thin <- strrep("-", w)
L <- character()
add <- function(...) L <<- c(L, glue(..., .envir = parent.frame()))
addRaw <- function(x) L <<- c(L, x)

addRaw(rule)
addRaw("  BAM PRE-FLIGHT QC")
add("  Project {projectName}  |  {nChecked} of {nExpected} samples checked  |  {Sys.Date()}")
addRaw(rule)
addRaw("")

verdictLine <- switch(overallVerdict,
    BLOCK = glue("  RESULT:  {nFail} of {nExpected} samples FAILED."),
    INCOMPLETE = glue("  RESULT:  {nIncomplete} of {nExpected} samples could not be assessed."),
    REVIEW = glue("  RESULT:  {nWarn} of {nExpected} samples need review."),
    CLEAR = glue("  RESULT:  all {nExpected} samples passed."))
addRaw(verdictLine)

if (nPairs > 0) add("           {nPairFail} of {nPairs} tumor/normal pairs unusable.")
addRaw(switch(overallVerdict,
    BLOCK = "           DO NOT RUN THE PIPELINE ON THE FAILED SAMPLES.",
    INCOMPLETE = "           RESOLVE THE MISSING METRICS BEFORE RUNNING.",
    REVIEW = "           Review the warnings below before committing compute.",
    CLEAR = "           Cleared to proceed."))
addRaw("")

if (nExpected != nChecked) {
    add("  WARNING: {nExpected - nChecked} sample(s) have incomplete metrics and were not")
    addRaw("  fully assessed. A partial report is worse than no report.")
    addRaw("")
}

if (nrow(droppedThresholds) > 0) {
    add("  NOTE: {nrow(droppedThresholds)} filter threshold(s) not evaluated, no data in this cohort:")
    add("  {str_c(droppedThresholds$metric, collapse = ', ')}")
    if ("samtools" %in% droppedSources) {
        addRaw("  Without samtools metrics there is no independent confirmation of the")
        addRaw("  Picard result. Verdicts rest on a single tool.")
    }
    addRaw("")
}

if (nFail > 0) {
    nFailN <- sum(failedSamples$sampleType == "N")
    nFailT <- sum(failedSamples$sampleType == "T")
    nTotalN <- sum(dat$sampleType == "N")
    add("  {nFailN} of {nTotalN} normals failed and {nFailT} tumor(s) failed.")
    if (nTotalN > 0 && nFailN == nTotalN && nTotalN > 0) {
        addRaw("  Every normal is affected, so no pair can be analysed as matched --")
        addRaw("  including pairs whose tumor is clean.")
    }
    addRaw("")

    addRaw(thin)
    addRaw("  FAILED SAMPLES                                            worst first")
    addRaw(thin)
    chimRef <- if (haveBackground) sprintf("%.2f%%", refMedian["pctChimeras"]) else "n/a"
    addRaw(glue("  SAMPLE                T/N  SPLIT-READS  READ-USED  COVERAGE  THRESHOLDS"))
    addRaw(glue("                             norm {chimRef}    norm 100%   usable   failed"))
    failedSamples |>
        mutate(line = sprintf("  %-21s %-4s %5.1f%% %5s %8s %8s%-4s %3d",
                              str_trunc(sample, 21),
                              tn,
                              pctChimeras,
                              if_else(is.na(chimeraFold), "", sprintf("(%.0fx)", chimeraFold)),
                              if_else(is.na(pctReadUsed), "n/a", sprintf("%.0f%%", pctReadUsed)),
                              if_else(is.na(meanCoverage), "n/a", sprintf("%.0fx", meanCoverage)),
                              if_else(lowCoverage, " LOW", ""),
                              nFail)) |>
        pull(line) |>
        walk(addRaw)
    addRaw("")
}

if (any(dat$lowCoverage)) {
    lowCov <- dat |> filter(lowCoverage) |> arrange(meanCoverage)
    addRaw("  Coverage advisory (separate from the filter thresholds above):")
    lowCov |>
        mutate(line = sprintf("    %-21s %-4s %.0fx usable, below the %.0fx floor for %s",
                              str_trunc(sample, 21), tn, meanCoverage,
                              coverageFloor, sampleType)) |>
        pull(line) |>
        walk(addRaw)
    addRaw("  These floors are not yet validated. Treat as a prompt to check, not a verdict.")
    addRaw("")
}

if (nWarn > 0) {
    addRaw(thin)
    addRaw("  SAMPLES NEEDING REVIEW")
    addRaw(thin)
    warnSamples |>
        mutate(line = sprintf("  %-21s %-4s %s", str_trunc(sample, 21), tn,
                              str_trunc(verdictReason, 45))) |>
        pull(line) |>
        walk(addRaw)
    addRaw("")
}

if (nIncomplete > 0) {
    addRaw(thin)
    addRaw("  SAMPLES NOT FULLY ASSESSED")
    addRaw(thin)
    incompleteSamples |>
        mutate(line = sprintf("  %-21s %-4s %s", str_trunc(sample, 21), tn,
                              str_trunc(verdictReason, 45))) |>
        pull(line) |>
        walk(addRaw)
    addRaw("")
}

if (nPass > 0) {
    addRaw(thin)
    addRaw("  PASSED SAMPLES")
    addRaw(thin)
    str_c(passSamples$sample, collapse = ", ") |>
        str_wrap(width = w - 4) |>
        str_split_1("\n") |>
        walk(\(x) addRaw(str_c("  ", x)))
    covRange <- passSamples |> filter(!is.na(meanCoverage)) |> pull(meanCoverage)
    if (length(covRange) > 0) {
        add("  All within normal ranges ({sprintf('%.0f', min(covRange))}-{sprintf('%.0f', max(covRange))}x usable coverage).")
    }
    addRaw("")
}

if (nPairs > 0) {
    addRaw(thin)
    addRaw("  PAIR CHECKS")
    addRaw(thin)
    ## Most patients contribute one pair, so the patient name alone reads best.
    ## Name the tumor only where a patient has several, which is the only case
    ## where the bare patient name would appear twice with different verdicts.
    pairs |>
        mutate(nForPatient = n(), .by = patient) |>
        mutate(label = if_else(nForPatient > 1,
                               str_c(patient, " [", sampleTumor, "]"), patient)) |>
        arrange(pairVerdict != "FAIL", patient) |>
        mutate(line = sprintf("  %-34s %-6s %s", str_trunc(label, 34),
                              pairVerdict, str_trunc(pairReason, 31))) |>
        pull(line) |>
        walk(addRaw)
    if (length(unpaired) > 0) {
        add("  Unpaired samples: {str_c(unpaired, collapse = ', ')}")
    }

    ## Insert divergence is reported separately because a sample-level failure
    ## takes precedence in the pair reason and would hide it. It is the one
    ## finding here with no loud failure mode: Facets runs to completion on
    ## mismatched pairs and returns copy number that is quietly wrong.
    divergent <- pairs |> filter(!is.na(insertRatio), insertRatio > INSERT_RATIO_FAIL)
    if (nrow(divergent) > 0) {
        addRaw("")
        add("  Insert size divergence in {nrow(divergent)} of {nPairs} pairs (over {INSERT_RATIO_FAIL}x):")
        divergent |>
            arrange(desc(insertRatio)) |>
            mutate(line = sprintf("    %-24s %.0f vs %.0f bases  (%.2fx)",
                                  str_trunc(patient, 24), insertSizeAverageTumor,
                                  insertSizeAverageNormal, insertRatio)) |>
            pull(line) |>
            walk(addRaw)
        addRaw("  Copy number calling on these pairs is unreliable and will not")
        addRaw("  announce itself. Consider a well-matched unmatched normal.")
    }
    addRaw("")
}

addRaw(thin)
addRaw("  WHAT TO DO")
addRaw(thin)
if (nFail > 0) {
    addRaw("  1. Do not commit compute to the failed samples.")
    addRaw("  2. Send the per-sample cards in the HTML report to the data provider.")
    addRaw("  3. Re-processing the existing data will not fix this. New data is needed.")
} else if (nWarn > 0) {
    addRaw("  1. Review the flagged samples against the reference ranges in the HTML.")
    addRaw("  2. Proceed if the deviations are understood and acceptable.")
} else {
    addRaw("  Nothing. The cohort is within expected ranges. Proceed.")
}
addRaw("")

if (haveBackground) {
    add("  Reference ranges from {max(refN, na.rm = TRUE)} previously mapped samples.")
} else {
    addRaw("  No background loaded. Thresholds are fixed values, no historical comparison.")
    addRaw("  Run bin/wgsTriageBackground.R to enable out-of-range detection.")
}
addRaw(rule)

consoleText <- as.character(L)
walk(consoleText, \(x) cat(x, "\n", sep = ""))
write_lines(consoleText, path(outDir, "preflightQC.txt"))

##
## Machine readable output, per section 7.5. Accumulating these is what turns
## the thresholds from a two-cohort estimate into a grounded range.
##
tsvOut <- dat |>
    select(project, sample, sampleType, patient, verdict, verdictReason,
           nFail, nWarn, nMissing, failedMetrics, warnedMetrics,
           any_of(c("pctChimeras", "pctSoftclip", "alignedFrac", "pctReadUsed", "supplementaryRate",
                    "pctProperlyPaired", "pctExcOverlap", "pctExcTotal",
                    "meanCoverage", "medianCoverage", "coverageFloor", "lowCoverage",
                    "insertSizeAverage", "pctImproperPairs", "totalReads",
                    "meanReadLength", "meanAlignedLength", "chimeraFold", "suppFold"))) |>
    arrange(verdict != "FAIL", desc(pctChimeras))

write_tsv(tsvOut, path(outDir, "preflightQC_samples.tsv"))
write_tsv(pairs, path(outDir, "preflightQC_pairs.tsv"))

##
## HTML report.
##
esc <- function(x) {
    x |> str_replace_all("&", "&amp;") |> str_replace_all("<", "&lt;") |>
        str_replace_all(">", "&gt;")
}

statusClass <- function(s) str_c("s", str_to_lower(s))

## One row of the cohort table: value, reference range, and verdict together.
## Design rule 3 -- every number carries its reference range in the same row.
metricCells <- function(sampleName) {
    thresholdResults |>
        filter(sample == sampleName) |>
        arrange(match(metric, THRESHOLDS$metric)) |>
        mutate(cell = pmap_chr(list(metric, value, status, units), \(m, v, s, u) {
            ref <- refMedian[m]
            refTxt <- if (is.na(ref)) "no background" else sprintf("norm %.2f%s", ref, u)
            valTxt <- if (is.na(v)) "n/a" else sprintf("%.2f%s", v, u)
            glue('<td class="{statusClass(s)}"><span class="v">{valTxt}</span>',
                 '<span class="r">{refTxt}</span></td>')
        })) |>
        pull(cell) |>
        str_c(collapse = "")
}

cohortRows <- dat |>
    arrange(match(verdict, c("FAIL", "INCOMPLETE", "WARN", "PASS")), desc(pctChimeras)) |>
    mutate(row = pmap_chr(list(sample, sampleType, verdict, meanCoverage, lowCoverage, coverageFloor),
        \(s, t, v, cov, low, floor) {
            covTxt <- if (is.na(cov)) "n/a" else sprintf("%.0fx", cov)
            covRef <- if (is.na(floor)) "" else sprintf("floor %.0fx", floor)
            glue('<tr><td class="name">{esc(s)}</td><td>{t}</td>',
                 '<td class="{statusClass(v)} verdict">{v}</td>{metricCells(s)}',
                 '<td class="{if (isTRUE(low)) "swarn" else "spass"}">',
                 '<span class="v">{covTxt}</span><span class="r">{covRef}</span></td></tr>')
        })) |>
    pull(row) |>
    str_c(collapse = "\n")

metricHeader <- usableThresholds |>
    mutate(h = glue('<th>{esc(label)}<span class="sub">{esc(PLAIN[metric])}</span></th>')) |>
    pull(h) |>
    str_c(collapse = "") |>
    str_c('<th>Usable coverage<span class="sub">Depth after quality filtering; advisory floor only</span></th>')

pairRows <- if (nPairs > 0) {
    pairs |>
        arrange(pairVerdict != "FAIL", patient) |>
        mutate(row = pmap_chr(list(patient, sampleTumor, sampleNormal, insertSizeAverageTumor,
                                   insertSizeAverageNormal, insertRatio, pairVerdict, pairReason),
            \(p, st, sn, it, inn, ir, v, why) {
                glue('<tr><td class="name">{esc(p)}</td><td>{esc(st)}</td><td>{esc(sn)}</td>',
                     '<td>{if (is.na(it)) "n/a" else sprintf("%.0f", it)}</td>',
                     '<td>{if (is.na(inn)) "n/a" else sprintf("%.0f", inn)}</td>',
                     '<td>{if (is.na(ir)) "n/a" else sprintf("%.2fx", ir)}</td>',
                     '<td class="{statusClass(v)} verdict">{v}</td><td>{esc(why)}</td></tr>')
            })) |>
        pull(row) |>
        str_c(collapse = "\n")
} else {
    '<tr><td colspan="8">No tumor/normal pairs could be inferred from the sample names.</td></tr>'
}

##
## Per-failure cards, in plain language, for the sequencing core.
##
failureCards <- if (nFail > 0) {
    failedSamples |>
        pull(sample) |>
        map_chr(\(s) {
            failing <- thresholdResults |> filter(sample == s, status == "FAIL")
            rows <- failing |>
                mutate(r = pmap_chr(list(metric, value, units), \(m, v, u) {
                    ref <- refMedian[m]
                    refTxt <- if (is.na(ref)) "no background available" else sprintf("%.2f%s", ref, u)
                    ## A multiplier is only worth printing when it is the
                    ## headline. "(1x)" next to a failing value reads as
                    ## reassurance and undercuts the row it sits in.
                    fold <- foldOf(m, v)
                    foldTxt <- if (is.na(fold) || fold < 1.5) "" else sprintf(" (%.0fx)", fold)
                    glue('<tr><td>{esc(PLAIN[m])}</td><td class="num">{sprintf("%.2f%s", v, u)}{foldTxt}</td>',
                         '<td class="num">{refTxt}</td></tr>')
                })) |>
                pull(r) |>
                str_c(collapse = "\n")
            actions <- failing |> pull(metric) |> (\(m) unique(ACTION[m]))() |>
                map_chr(\(a) glue("<li>{esc(a)}</li>")) |> str_c(collapse = "\n")
            row <- dat |> filter(sample == s)
            glue('
<div class="card">
  <h3>{esc(s)} <span class="sfail verdict">FAILED</span></h3>
  <p class="lead">This sample did not pass {nrow(failing)} of the {nrow(usableThresholds)} checks applied.</p>
  <table class="cardtable">
    <thead><tr><th>What we measured</th><th>This sample</th><th>Normal</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
  <p class="lead">Usable coverage after quality filtering: <b>{if (is.na(row$meanCoverage)) "n/a" else sprintf("%.0fx", row$meanCoverage)}</b>.</p>
  <p class="what">What this means</p>
  <ul>{actions}</ul>
  <p class="caveat">We can measure what the sequence data looks like, but not what
  produced it. No conclusion is drawn here about sample handling or preparation.</p>
</div>')
        }) |>
        str_c(collapse = "\n")
} else {
    '<p class="lead">No samples failed. Nothing to report to the data provider.</p>'
}

bannerClass <- switch(overallVerdict, BLOCK = "block", INCOMPLETE = "block",
                      REVIEW = "review", CLEAR = "clear")
bannerText <- switch(overallVerdict,
    BLOCK = glue("{nFail} of {nExpected} samples FAILED -- do not run the pipeline on them"),
    INCOMPLETE = glue("{nIncomplete} of {nExpected} samples could not be assessed"),
    REVIEW = glue("{nWarn} of {nExpected} samples need review before committing compute"),
    CLEAR = glue("All {nExpected} samples passed -- cleared to proceed"))

thresholdRows <- usableThresholds |>
    mutate(r = glue('<tr><td>{esc(label)}</td><td>{esc(PLAIN[metric])}</td>',
                    '<td class="num">{if_else(is.na(fail), "--", as.character(fail))}{units}</td>',
                    '<td class="num">{if_else(is.na(warn), "--", as.character(warn))}{units}</td>',
                    '<td class="num">{if_else(is.na(refMedian[metric]), "n/a", sprintf("%.3f", refMedian[metric]))}</td>',
                    '<td class="num">{if_else(is.na(refN[metric]), "0", as.character(refN[metric]))}</td></tr>')) |>
    pull(r) |>
    str_c(collapse = "\n")

html <- glue('<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pre-flight QC -- {esc(projectName)}</title>
<style>
:root {{ color-scheme: light dark; }}
* {{ box-sizing: border-box; }}
body {{ font: 15px/1.55 -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
  margin: 0; padding: 2rem 1.25rem 4rem; background: #fbfbfc; color: #1a1a1c; }}
.wrap {{ max-width: 1180px; margin: 0 auto; }}
h1 {{ font-size: 1.5rem; margin: 0 0 .25rem; }}
h2 {{ font-size: 1.15rem; margin: 2.5rem 0 .75rem; padding-bottom: .35rem;
  border-bottom: 1px solid #e2e2e6; }}
h3 {{ font-size: 1rem; margin: 0 0 .5rem; }}
.meta {{ color: #63636b; font-size: .875rem; margin-bottom: 1.25rem; }}
.banner {{ padding: 1rem 1.25rem; border-radius: 8px; font-weight: 600; font-size: 1.05rem;
  margin: 1rem 0 1.5rem; border-left: 5px solid; }}
.banner.block {{ background: #fdeaea; border-color: #c0392b; color: #7f1d13; }}
.banner.review {{ background: #fdf5e3; border-color: #d19a0a; color: #7a5901; }}
.banner.clear {{ background: #eaf7ee; border-color: #2e8b4f; color: #1c5c33; }}
.counts {{ display: flex; flex-wrap: wrap; gap: 1.5rem; margin: 0 0 1rem; padding: 0; list-style: none; }}
.counts li {{ font-size: .875rem; color: #63636b; }}
.counts b {{ display: block; font-size: 1.5rem; color: #1a1a1c; font-weight: 650; }}
.scroll {{ overflow-x: auto; -webkit-overflow-scrolling: touch; }}
table {{ border-collapse: collapse; width: 100%; font-size: .8125rem; }}
th, td {{ padding: .45rem .55rem; text-align: left; border-bottom: 1px solid #ececf0;
  vertical-align: top; white-space: nowrap; }}
th {{ background: #f2f2f5; font-weight: 600; font-size: .75rem; }}
th .sub {{ display: block; font-weight: 400; color: #71717a; font-size: .6875rem;
  white-space: normal; max-width: 12rem; }}
td.name {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-weight: 600; }}
td.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
.v {{ display: block; font-variant-numeric: tabular-nums; font-weight: 600; }}
.r {{ display: block; font-size: .6875rem; color: #71717a; }}
.verdict {{ font-weight: 700; font-size: .75rem; letter-spacing: .02em; }}
.sfail {{ background: #fdeaea; color: #a5281a; }}
.swarn {{ background: #fdf5e3; color: #8a6400; }}
.spass {{ background: transparent; }}
.smissing, .sincomplete {{ background: #eeeef2; color: #52525b; }}
.card {{ border: 1px solid #e2e2e6; border-radius: 8px; padding: 1.1rem 1.25rem;
  margin-bottom: 1rem; background: #fff; }}
.cardtable {{ margin: .5rem 0 .85rem; font-size: .8125rem; }}
.cardtable th {{ background: transparent; border-bottom: 1px solid #d4d4d8; }}
.lead {{ margin: .35rem 0; }}
.what {{ font-weight: 600; margin: .75rem 0 .25rem; }}
.card ul {{ margin: .25rem 0 .5rem; padding-left: 1.2rem; }}
.card li {{ margin-bottom: .25rem; }}
.caveat {{ font-size: .8125rem; color: #71717a; border-top: 1px solid #ececf0;
  padding-top: .6rem; margin-top: .75rem; }}
.note {{ background: #f2f2f5; border-radius: 6px; padding: .75rem 1rem; font-size: .875rem;
  margin: 1rem 0; }}
code {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: .8125rem; }}
@media (prefers-color-scheme: dark) {{
  body {{ background: #121214; color: #e8e8ea; }}
  h2 {{ border-color: #2c2c32; }}
  .meta, .counts li, .r, th .sub, .caveat {{ color: #a1a1aa; }}
  .counts b {{ color: #e8e8ea; }}
  th {{ background: #1d1d21; }}
  th, td {{ border-color: #2c2c32; }}
  .banner.block {{ background: #3b1614; border-color: #e05545; color: #ffb4a8; }}
  .banner.review {{ background: #3a2e10; border-color: #d19a0a; color: #f5cf6a; }}
  .banner.clear {{ background: #12301d; border-color: #3fa564; color: #92dca9; }}
  .sfail {{ background: #3b1614; color: #ffb4a8; }}
  .swarn {{ background: #3a2e10; color: #f5cf6a; }}
  .smissing, .sincomplete {{ background: #26262b; color: #a1a1aa; }}
  .card {{ background: #18181b; border-color: #2c2c32; }}
  .cardtable th {{ border-color: #3f3f46; }}
  .note {{ background: #1d1d21; }}
}}
</style>
</head>
<body>
<div class="wrap">

<h1>BAM pre-flight QC</h1>
<p class="meta">Project <b>{esc(projectName)}</b> &nbsp;|&nbsp; {nChecked} of {nExpected} samples checked
 &nbsp;|&nbsp; {Sys.Date()} &nbsp;|&nbsp; source: <code>{esc(mapDir)}</code></p>

<div class="banner {bannerClass}">{bannerText}</div>

<ul class="counts">
  <li><b>{nExpected}</b>samples submitted</li>
  <li><b>{nChecked}</b>fully assessed</li>
  <li><b>{nFail}</b>failed</li>
  <li><b>{nWarn}</b>need review</li>
  <li><b>{nPass}</b>passed</li>
  <li><b>{nPairFail} / {nPairs}</b>pairs unusable</li>
</ul>

{if (nExpected != nChecked) glue(\'<div class="note"><b>Cohort is incomplete.</b> {nExpected - nChecked} sample(s) are missing metrics and were not fully assessed. A report covering a subset is how the original investigation reached a wrong conclusion.</div>\') else ""}
{if (nrow(droppedThresholds) > 0) glue(\'<div class="note"><b>Reduced check set.</b> {nrow(droppedThresholds)} filter threshold(s) were not evaluated because this cohort carries no data for them: <code>{str_c(droppedThresholds$metric, collapse = ", ")}</code>.{if ("samtools" %in% droppedSources) " Without samtools metrics there is no independent confirmation of the Picard result, so every verdict below rests on a single tool." else ""}</div>\') else ""}

<h2>Cohort</h2>
<p class="meta">Worst first. Each cell shows the measured value above the clean-sample reference.</p>
<div class="scroll">
<table>
<thead><tr><th>Sample</th><th>T/N</th><th>Verdict</th>{metricHeader}</tr></thead>
<tbody>
{cohortRows}
</tbody>
</table>
</div>

<h2>Tumor / normal pairs</h2>
<p class="meta">Pairing is inferred from sample names. Two samples can each pass on their
own and still be unusable together: Facets needs comparable insert size distributions,
and unlike the SV caller it fails silently rather than crashing.</p>
<div class="scroll">
<table>
<thead><tr><th>Patient</th><th>Tumor</th><th>Normal</th><th>Insert T</th><th>Insert N</th>
<th>Ratio</th><th>Verdict</th><th>Reason</th></tr></thead>
<tbody>
{pairRows}
</tbody>
</table>
</div>

<h2>Per-sample detail for the data provider</h2>
{failureCards}

<h2>Thresholds applied</h2>
<div class="scroll">
<table>
<thead><tr><th>Metric</th><th>Plain language</th><th>Fail</th><th>Warn</th>
<th>Clean median</th><th>Background n</th></tr></thead>
<tbody>
{thresholdRows}
</tbody>
</table>
</div>
<p class="meta">A sample tripping {WARN_ESCALATION} or more warnings at once is failed even when no single
threshold fires, since these metrics move together under genuine degradation.
Coverage floors ({COVERAGE_WARN[["N"]]}x normal, {COVERAGE_WARN[["T"]]}x tumor) are advisory and not yet validated.</p>

</div>
</body>
</html>
', .open = "{", .close = "}")

write_lines(html, path(outDir, "preflightQC.html"))

cat("\n")
cat(sprintf("  Wrote %s\n", path(outDir, "preflightQC.txt")))
cat(sprintf("  Wrote %s\n", path(outDir, "preflightQC.html")))
cat(sprintf("  Wrote %s\n", path(outDir, "preflightQC_samples.tsv")))
cat(sprintf("  Wrote %s\n", path(outDir, "preflightQC_pairs.tsv")))

## Advisory by decision: always exit 0. See the header comment.
quit(status = 0)
