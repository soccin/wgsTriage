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
wgsTriage.R {WGSTRIAGE_VERSION} -- post-mapping QC filter. Reads only what the
Map stage already wrote and renders a verdict per sample and per tumor/normal
pair.

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
  --out <OutDir>         Directory for the report. Default: ./wgsTriage_out
  --project <Name>       Label shown in the report.
                         Default: the directory containing <MapDir>.
  -h, --help             Show this message and exit.

Writes into <OutDir>:
  wgsTriage.txt          the console report, as text
  wgsTriage.html         the same report, formatted
  wgsTriage_samples.tsv  per-sample metrics and verdicts
  wgsTriage_pairs.tsv    per tumor/normal pair verdicts

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
outDir <- getOpt("--out", "wgsTriage_out")
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
    pctChimeras = "Structural variant calling will fail or produce false positives.",
    supplementaryRate = "Confirms the split-read finding from a second independent tool.",
    pctSoftclip = "Alignment is discarding a large fraction of each read. Usable data is below the delivered figure.",
    pctReadUsed = "Most of each sequenced read is discarded at alignment. Effective yield is a fraction of the raw yield.",
    pctProperlyPaired = "Paired-end structure is degraded. Insert size and copy number estimates will be unreliable.",
    pctExcOverlap = "Fragments are shorter than the read length, so mates overlap and are counted once.",
    pctExcTotal = "Usable coverage is below raw coverage. Externally reported depth overstates what is usable.")

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
addRaw("")

if (nExpected != nChecked) {
    add("  WARNING: {nExpected - nChecked} sample(s) have incomplete metrics and were not")
    addRaw("  fully assessed.")
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
    addRaw("  These floors are not yet validated.")
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
        addRaw("  Copy number calling on these pairs is unreliable and does not")
        addRaw("  fail visibly.")
    }
    addRaw("")
}

if (haveBackground) {
    add("  Reference ranges from {max(refN, na.rm = TRUE)} previously mapped samples.")
} else {
    addRaw("  No background loaded. Thresholds are fixed values, no historical comparison.")
    addRaw("  Run bin/wgsTriageBackground.R to enable out-of-range detection.")
}
addRaw(rule)

consoleText <- as.character(L)
walk(consoleText, \(x) cat(x, "\n", sep = ""))
write_lines(consoleText, path(outDir, "wgsTriage.txt"))

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

write_tsv(tsvOut, path(outDir, "wgsTriage_samples.tsv"))
write_tsv(pairs, path(outDir, "wgsTriage_pairs.tsv"))

##
## HTML report.
##
esc <- function(x) {
    x |> str_replace_all("&", "&amp;") |> str_replace_all("<", "&lt;") |>
        str_replace_all(">", "&gt;")
}

##
## The Map directory is usually given relative ("../../Map"), which identifies
## nothing once the report has been mailed away from where it was produced.
## Resolve it and trim to the first "Users" component, where the informative
## part of these paths begins. Display only: nothing is resolved from this
## string, so the relative path stays the one the tool actually reads.
##
displayPath <- function(p) {
    parts <- path_split(path_real(p))[[1]]
    i <- which(parts == "Users")
    if (length(i) > 0) path_join(parts[i[1]:length(parts)]) else path_real(p)
}

sourcePath <- displayPath(mapDir)

statusClass <- function(s) str_c("s", str_to_lower(s))

## One row of the cohort table: value, reference range, and verdict together.
## Design rule 3 -- every number carries its reference range in the same row.
metricCells <- function(sampleName) {
    thresholdResults |>
        filter(sample == sampleName) |>
        arrange(match(metric, THRESHOLDS$metric)) |>
        mutate(cell = pmap_chr(list(metric, value, status, units, label), \(m, v, s, u, lab) {
            ref <- refMedian[m]
            refTxt <- if (is.na(ref)) "no background" else sprintf("norm %.2f%s", ref, u)
            valTxt <- if (is.na(v)) "n/a" else sprintf("%.2f%s", v, u)
            glue('<td class="{statusClass(s)}" data-label="{esc(lab)}"><span class="v">{valTxt}</span>',
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
            glue('<tr><td class="name" data-label="Sample">{esc(s)}</td>',
                 '<td data-label="T/N">{t}</td>',
                 '<td class="{statusClass(v)} verdict" data-label="Verdict">{v}</td>{metricCells(s)}',
                 '<td class="{if (isTRUE(low)) "swarn" else "spass"}" data-label="Coverage">',
                 '<span class="v">{covTxt}</span><span class="r">{covRef}</span></td></tr>')
        })) |>
    pull(row) |>
    str_c(collapse = "\n")

metricHeader <- usableThresholds |>
    mutate(h = glue('<th>{esc(label)}<span class="sub">{esc(PLAIN[metric])}</span></th>')) |>
    pull(h) |>
    str_c(collapse = "") |>
    str_c('<th>Coverage<span class="sub">Depth after quality filtering; advisory floor only</span></th>')

pairRows <- if (nPairs > 0) {
    pairs |>
        arrange(pairVerdict != "FAIL", patient) |>
        mutate(row = pmap_chr(list(patient, sampleTumor, sampleNormal, insertSizeAverageTumor,
                                   insertSizeAverageNormal, insertRatio, pairVerdict, pairReason),
            \(p, st, sn, it, inn, ir, v, why) {
                glue('<tr><td class="name" data-label="Patient">{esc(p)}</td>',
                     '<td data-label="Tumor">{esc(st)}</td><td data-label="Normal">{esc(sn)}</td>',
                     '<td data-label="Insert T">{if (is.na(it)) "n/a" else sprintf("%.0f", it)}</td>',
                     '<td data-label="Insert N">{if (is.na(inn)) "n/a" else sprintf("%.0f", inn)}</td>',
                     '<td data-label="Ratio">{if (is.na(ir)) "n/a" else sprintf("%.2fx", ir)}</td>',
                     '<td class="{statusClass(v)} verdict" data-label="Verdict">{v}</td>',
                     '<td data-label="Reason">{esc(why)}</td></tr>')
            })) |>
        pull(row) |>
        str_c(collapse = "\n")
} else {
    '<tr><td colspan="8" data-label="Pairs">No tumor/normal pairs could be inferred from the sample names.</td></tr>'
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
  <p class="what">Effect on downstream analysis</p>
  <ul>{actions}</ul>
  <p class="caveat">We can measure what the sequence data looks like, but not what
  produced it. No conclusion is drawn here about sample handling or preparation.</p>
</div>')
        }) |>
        str_c(collapse = "\n")
} else {
    '<p class="lead">No samples failed.</p>'
}

bannerClass <- switch(overallVerdict, BLOCK = "block", INCOMPLETE = "block",
                      REVIEW = "review", CLEAR = "clear")
bannerText <- switch(overallVerdict,
    BLOCK = glue("{nFail} of {nExpected} samples FAILED"),
    INCOMPLETE = glue("{nIncomplete} of {nExpected} samples could not be assessed"),
    REVIEW = glue("{nWarn} of {nExpected} samples outside reference range"),
    CLEAR = glue("All {nExpected} samples passed"))

thresholdRows <- usableThresholds |>
    mutate(r = glue('<tr><td data-label="Metric">{esc(label)}</td>',
                    '<td data-label="Plain language">{esc(PLAIN[metric])}</td>',
                    '<td class="num" data-label="Fail">{if_else(is.na(fail), "--", as.character(fail))}{units}</td>',
                    '<td class="num" data-label="Warn">{if_else(is.na(warn), "--", as.character(warn))}{units}</td>',
                    '<td class="num" data-label="Clean median">{if_else(is.na(refMedian[metric]), "n/a", sprintf("%.3f", refMedian[metric]))}</td>',
                    '<td class="num" data-label="Background n">{if_else(is.na(refN[metric]), "0", as.character(refN[metric]))}</td></tr>')) |>
    pull(r) |>
    str_c(collapse = "\n")

##
## Stylesheet. Held outside the glue template because it is entirely static and
## every literal brace would otherwise have to be doubled, which is how a
## stylesheet acquires a syntax error nobody can see.
##
## Colours are custom properties so light and dark differ only in the values of
## one block, not in a duplicate copy of every rule. Light is the default and
## prefers-color-scheme is deliberately not consulted: this report is a document
## that gets mailed to a data provider, and one that arrives dark because of the
## reader's operating system setting is a surprise. The toggle is there for
## anyone who wants dark.
##
styleBlock <- '
:root {
  color-scheme: light;
  --bg: #fbfbfc; --fg: #1a1a1c; --muted: #63636b; --muted2: #71717a;
  --line: #e2e2e6; --lineSoft: #ececf0; --thBg: #f2f2f5; --cardBg: #fff;
  --cardLine: #d4d4d8;
  --failBg: #fdeaea; --failFg: #a5281a; --failBd: #c0392b; --failText: #7f1d13;
  --warnBg: #fdf5e3; --warnFg: #8a6400; --warnBd: #d19a0a; --warnText: #7a5901;
  --passBg: #eaf7ee; --passBd: #2e8b4f; --passText: #1c5c33;
  --missBg: #eeeef2; --missFg: #52525b;
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --bg: #121214; --fg: #e8e8ea; --muted: #a1a1aa; --muted2: #a1a1aa;
  --line: #2c2c32; --lineSoft: #2c2c32; --thBg: #1d1d21; --cardBg: #18181b;
  --cardLine: #3f3f46;
  --failBg: #3b1614; --failFg: #ffb4a8; --failBd: #e05545; --failText: #ffb4a8;
  --warnBg: #3a2e10; --warnFg: #f5cf6a; --warnBd: #d19a0a; --warnText: #f5cf6a;
  --passBg: #12301d; --passBd: #3fa564; --passText: #92dca9;
  --missBg: #26262b; --missFg: #a1a1aa;
}
* { box-sizing: border-box; }
body { font: 15px/1.55 -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
  margin: 0; padding: 2rem 1.25rem 4rem; background: var(--bg); color: var(--fg); }
.wrap { max-width: 1180px; margin: 0 auto; }
h1 { font-size: 1.5rem; margin: 0 0 .25rem; }
h2 { font-size: 1.15rem; margin: 2.5rem 0 .75rem; padding-bottom: .35rem;
  border-bottom: 1px solid var(--line); }
h3 { font-size: 1rem; margin: 0 0 .5rem; }
.meta { color: var(--muted); font-size: .875rem; margin-bottom: 1.25rem; }
.head { display: flex; justify-content: space-between; align-items: flex-start; gap: 1rem; }
.themeToggle { flex: none; font: inherit; font-size: .8125rem; cursor: pointer;
  background: var(--thBg); color: var(--fg); border: 1px solid var(--line);
  border-radius: 6px; padding: .35rem .7rem; }
.themeToggle:hover { border-color: var(--muted); }
.banner { padding: 1rem 1.25rem; border-radius: 8px; font-weight: 600; font-size: 1.05rem;
  margin: 1rem 0 1.5rem; border-left: 5px solid; }
.banner.block { background: var(--failBg); border-color: var(--failBd); color: var(--failText); }
.banner.review { background: var(--warnBg); border-color: var(--warnBd); color: var(--warnText); }
.banner.clear { background: var(--passBg); border-color: var(--passBd); color: var(--passText); }
.counts { display: flex; flex-wrap: wrap; gap: 1.5rem; margin: 0 0 1rem; padding: 0; list-style: none; }
.counts li { font-size: .875rem; color: var(--muted); }
.counts b { display: block; font-size: 1.5rem; color: var(--fg); font-weight: 650; }
table { border-collapse: collapse; width: 100%; font-size: .8125rem; table-layout: auto; }
th, td { padding: .45rem .55rem; text-align: left; border-bottom: 1px solid var(--lineSoft);
  vertical-align: top; }
th { background: var(--thBg); font-weight: 600; font-size: .75rem; }
th .sub { display: block; font-weight: 400; color: var(--muted2); font-size: .6875rem;
  max-width: 9rem; }
/*
 * Sample names are one unbreakable token. break-all split them mid-name
 * (APTL_MDA009_ / N01) because it also tells the auto table layout the column
 * can shrink to a single character. nowrap makes the longest name set the
 * column width instead; min-width keeps it off the neighbouring column.
 */
td.name { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-weight: 600;
  white-space: nowrap; min-width: 9.5rem; }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
.v { display: block; font-variant-numeric: tabular-nums; font-weight: 600; white-space: nowrap; }
.r { display: block; font-size: .6875rem; color: var(--muted2); }
.verdict { font-weight: 700; font-size: .75rem; letter-spacing: .02em; }
.sfail { background: var(--failBg); color: var(--failFg); }
.swarn { background: var(--warnBg); color: var(--warnFg); }
.spass { background: transparent; }
.smissing, .sincomplete { background: var(--missBg); color: var(--missFg); }
.card { border: 1px solid var(--line); border-radius: 8px; padding: 1.1rem 1.25rem;
  margin-bottom: 1rem; background: var(--cardBg); }
.cardtable { margin: .5rem 0 .85rem; font-size: .8125rem; }
.cardtable th { background: transparent; border-bottom: 1px solid var(--cardLine); }
.lead { margin: .35rem 0; }
.what { font-weight: 600; margin: .75rem 0 .25rem; }
.card ul { margin: .25rem 0 .5rem; padding-left: 1.2rem; }
.card li { margin-bottom: .25rem; }
.caveat { font-size: .8125rem; color: var(--muted2); border-top: 1px solid var(--lineSoft);
  padding-top: .6rem; margin-top: .75rem; }
.note { background: var(--thBg); border-radius: 6px; padding: .75rem 1rem; font-size: .875rem;
  margin: 1rem 0; }
code { font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: .8125rem;
  word-break: break-all; }

/*
 * Narrow screens: every table becomes one block per row, each cell labelled
 * from its data-label. The tables are never allowed to scroll sideways --
 * a browser gives no hint that a table continues past the edge, so anything
 * out there is simply lost.
 */
@media (max-width: 980px) {
  thead { position: absolute; width: 1px; height: 1px; padding: 0; margin: -1px;
    overflow: hidden; clip: rect(0 0 0 0); white-space: nowrap; border: 0; }
  table, tbody, tr, td { display: block; width: 100%; }
  tr { border: 1px solid var(--line); border-radius: 8px; margin-bottom: .75rem;
    padding: .25rem .6rem; }
  td { border-bottom: 1px solid var(--lineSoft); padding: .45rem .1rem; }
  tr td:last-child { border-bottom: 0; }
  td::before { content: attr(data-label); display: block; font-size: .6875rem;
    font-weight: 600; color: var(--muted2); text-transform: uppercase;
    letter-spacing: .03em; }
  td.num { text-align: left; }
  /* The column no longer exists here, so nothing needs to be held open --
     and a nowrap name is exactly what would push a phone sideways. */
  td.name { white-space: normal; overflow-wrap: break-word; min-width: 0; }
  .counts { gap: 1rem; }
}

/*
 * Print. The failure cards are meant to be sent to a data provider, so the
 * status colours have to survive the trip: they are the whole point of the
 * page and a greyscale FAIL is just a number.
 */
@media print {
  :root { color-scheme: light; }
  body { background: #fff; padding: 0; font-size: 11pt; }
  .themeToggle { display: none; }
  .banner, .sfail, .swarn, .smissing, .sincomplete, th, .note {
    -webkit-print-color-adjust: exact; print-color-adjust: exact; }
  .card, tr, .note { break-inside: avoid; page-break-inside: avoid; }
  h2 { break-after: avoid; page-break-after: avoid; }
  a[href]::after { content: ""; }
}
'

##
## Theme toggle. Light unless the reader asks otherwise; the choice is
## remembered per browser. localStorage throws on a file:// page in some
## browsers, which is why both halves are wrapped -- a failure there costs the
## reader a remembered preference, not a working button.
##
scriptBlock <- '
(function () {
  var root = document.documentElement;
  var btn = document.getElementById("themeToggle");
  var saved = null;
  try { saved = localStorage.getItem("wgsTriageTheme"); } catch (e) {}
  if (saved === "dark" || saved === "light") root.setAttribute("data-theme", saved);
  function label() {
    btn.textContent = root.getAttribute("data-theme") === "dark" ? "Light mode" : "Dark mode";
  }
  label();
  btn.addEventListener("click", function () {
    var next = root.getAttribute("data-theme") === "dark" ? "light" : "dark";
    root.setAttribute("data-theme", next);
    try { localStorage.setItem("wgsTriageTheme", next); } catch (e) {}
    label();
  });
})();
'

html <- glue('<!doctype html>
<html lang="en" data-theme="light">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pre-flight QC -- {esc(projectName)}</title>
<style>
{styleBlock}
</style>
</head>
<body>
<div class="wrap">

<div class="head">
<h1>BAM pre-flight QC</h1>
<button type="button" id="themeToggle" class="themeToggle">Dark mode</button>
</div>
<p class="meta">Project <b>{esc(projectName)}</b> &nbsp;|&nbsp; source: <code>{esc(sourcePath)}</code><br>{nChecked} of {nExpected} samples checked<br>{Sys.Date()} &nbsp;|&nbsp; wgsTriage {WGSTRIAGE_VERSION}</p>

<div class="banner {bannerClass}">{bannerText}</div>

<ul class="counts">
  <li><b>{nExpected}</b>samples submitted</li>
  <li><b>{nChecked}</b>fully assessed</li>
  <li><b>{nFail}</b>failed</li>
  <li><b>{nWarn}</b>need review</li>
  <li><b>{nPass}</b>passed</li>
  <li><b>{nPairFail} / {nPairs}</b>pairs unusable</li>
</ul>

{if (nExpected != nChecked) glue(\'<div class="note"><b>Cohort is incomplete.</b> {nExpected - nChecked} sample(s) are missing metrics and were not fully assessed.</div>\') else ""}
{if (nrow(droppedThresholds) > 0) glue(\'<div class="note"><b>Reduced check set.</b> {nrow(droppedThresholds)} filter threshold(s) were not evaluated because this cohort carries no data for them: <code>{str_c(droppedThresholds$metric, collapse = ", ")}</code>.{if ("samtools" %in% droppedSources) " Without samtools metrics there is no independent confirmation of the Picard result, so every verdict below rests on a single tool." else ""}</div>\') else ""}

<h2>Cohort</h2>
<p class="meta">Worst first. Each cell shows the measured value above the clean-sample reference.</p>
<table>
<thead><tr><th>Sample</th><th>T/N</th><th>Verdict</th>{metricHeader}</tr></thead>
<tbody>
{cohortRows}
</tbody>
</table>

<h2>Tumor / normal pairs</h2>
<p class="meta">Pairing is inferred from sample names. Two samples can each pass on their
own and still be unusable together: Facets needs comparable insert size distributions,
and unlike the SV caller it fails silently rather than crashing.</p>
<table>
<thead><tr><th>Patient</th><th>Tumor</th><th>Normal</th><th>Insert T</th><th>Insert N</th>
<th>Ratio</th><th>Verdict</th><th>Reason</th></tr></thead>
<tbody>
{pairRows}
</tbody>
</table>

<h2>Per-sample detail for the data provider</h2>
{failureCards}

<h2>Thresholds applied</h2>
<table>
<thead><tr><th>Metric</th><th>Plain language</th><th>Fail</th><th>Warn</th>
<th>Clean median</th><th>Background n</th></tr></thead>
<tbody>
{thresholdRows}
</tbody>
</table>
<p class="meta">A sample tripping {WARN_ESCALATION} or more warnings at once is failed even when no single
threshold fires, since these metrics move together under genuine degradation.
Coverage floors ({COVERAGE_WARN[["N"]]}x normal, {COVERAGE_WARN[["T"]]}x tumor) are advisory and not yet validated.</p>

</div>
<script>{scriptBlock}</script>
</body>
</html>
', .open = "{", .close = "}")

write_lines(html, path(outDir, "wgsTriage.html"))

cat("\n")
cat(sprintf("  Wrote %s\n", path(outDir, "wgsTriage.txt")))
cat(sprintf("  Wrote %s\n", path(outDir, "wgsTriage.html")))
cat(sprintf("  Wrote %s\n", path(outDir, "wgsTriage_samples.tsv")))
cat(sprintf("  Wrote %s\n", path(outDir, "wgsTriage_pairs.tsv")))

## Advisory by decision: always exit 0. See the header comment.
quit(status = 0)
