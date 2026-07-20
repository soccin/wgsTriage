##
## Shared parsing and gating logic for post-mapping pre-flight QC.
## Sourced by bin/wgsTriage.R and bin/wgsTriageBackground.R.
##
## Inputs are restricted to what the Map stage already produces:
##   Map/out/metrics/<sample>/<sample>.asm.txt   Picard CollectAlignmentSummaryMetrics
##   Map/out/metrics/<sample>/<sample>.wgs.txt   Picard CollectWgsMetrics
##   Map/sbam/multiqc/multiqc_data/multiqc_samtools_stats.txt
##
## Qualimap is deliberately not read. It overstates usable coverage by up to
## 2.6x on degraded samples and silently covered only 12 of 16 samples in
## Proj_16840_N, which produced a materially wrong conclusion.
##

##
## Single source of truth for the version. Reported by --help so a report can
## be traced back to the threshold set that produced it.
##
WGSTRIAGE_VERSION <- "0.9.0"

suppressPackageStartupMessages({
    library(tidyverse)
    library(fs)
    library(glue)
})

##
## Filter threshold definitions.
##
## Starting point was NORMAL_BAM_QC_REPORT.md section 5, which reported a 45x
## separation with nothing in between and proposed FAIL above 1.0% chimeras.
## That held across the 27 samples it examined. It does not hold across the 421
## sample archive, which shows a sparse but real borderline band:
##
##   below 0.5%   283 samples      the normal population
##   0.5 to 1%    107 samples      still normal, upper tail
##   1 to 5%        9 samples      mildly degraded, cause unknown
##   above 5%      22 samples      the catastrophic class
##
## So FAIL moves to 5.0% and 1.0% becomes WARN. Every sample in the
## Proj_16840_N disaster ran 8.5 to 17.3% and still fails outright, while the
## nine borderline samples are surfaced for review rather than blocked. A filter
## that fires on 9% of all historical samples gets ignored, and an ignored filter
## is the failure mode this whole exercise exists to prevent.
##
## The samtools thresholds are the independent confirmation: supplementary alignment
## rate measures essentially the same physical defect as PCT_CHIMERAS but is
## computed by a different tool on a different pass over the data. It is also
## the tightest metric available, spanning only 0.094 to 0.160% across every
## clean sample in the archive, so a 1.0% threshold sits 6x above the observed
## maximum. Agreement between the two is what makes a FAIL verdict defensible
## to someone upstream who does not want to hear it.
##
## `label` is the column heading in the HTML report and nothing else, so it is
## kept short on purpose: a long heading widens its column past the width of the
## data under it, and eleven such columns push the table off the screen. The
## explanation belongs in the sub-heading, which is what PLAIN in bin/wgsTriage.R
## supplies.
THRESHOLDS <- tribble(
    ~metric,               ~label,           ~units, ~direction, ~fail, ~warn, ~source,
    "pctChimeras",         "Split pairs",       "%",     "high",   5.0,   1.0, "picard_asm",
    "supplementaryRate",   "Multi-mapped",      "%",     "high",   1.0,   0.5, "samtools",
    "pctSoftclip",         "Clipped bases",     "%",     "high",   5.0,   1.0, "picard_asm",
    "pctReadUsed",         "Read used",         "%",     "low",   95.0,  98.0, "picard_asm",
    "pctProperlyPaired",   "Proper pairs",      "%",     "low",   97.0,  98.5, "samtools",
    "pctExcOverlap",       "Overlap loss",      "%",     "high",    NA,  10.0, "picard_wgs",
    "pctExcTotal",         "Total loss",        "%",     "high",    NA,  35.0, "picard_wgs")

## A sample tripping this many warnings is treated as a failure even when no
## single threshold fires. These metrics are not independent: genuine degradation
## pushes all of them at once, and three simultaneous warnings is a pattern
## rather than a coincidence. Without this, a sample can sit just under every
## threshold on every metric and still be unusable.
WARN_ESCALATION <- 3

## Coverage floors are advisory and stratified by sample class, since tumors are
## sequenced far deeper than normals. Section 5.8 flags these as unvalidated.
COVERAGE_WARN <- c(N = 25, T = 50, unknown = 25)

## Tumor and normal insert size distributions must be comparable or Facets
## produces a result that is wrong without crashing. Section 5.6; the 1.5x
## tolerance is a starting point and has not been confirmed with the authors.
INSERT_RATIO_FAIL <- 1.5

##
## Picard metrics files carry '#' comment headers, then a '## METRICS CLASS'
## line, then a header row, then data rows, then a blank line, then optionally a
## '## HISTOGRAM' block. Parse structurally. Line offsets are not stable across
## Picard versions and a fixed skip= is how this breaks silently later.
##
readPicardMetrics <- function(path) {
    if (!file_exists(path)) return(NULL)
    lines <- read_lines(path, progress = FALSE)
    start <- which(str_detect(lines, "^## METRICS CLASS"))
    if (length(start) == 0) return(NULL)
    body <- lines[(start[1] + 1):length(lines)]
    blank <- which(str_trim(body) == "")
    if (length(blank) > 0) body <- body[seq_len(blank[1] - 1)]
    if (length(body) < 2) return(NULL)
    read_tsv(I(body), show_col_types = FALSE, progress = FALSE, na = c("", "NA", "?"))
}

## Sample name is the file basename with the metric extension and any
## intermediate pipeline suffix removed. Both layouts seen in the archive are
## handled: metrics/<sample>/<sample>.asm.txt and metrics/<sample>.md.wgs.txt
sampleFromMetricsPath <- function(path) {
    path_file(path) |>
        str_remove("\\.(asm|wgs)\\.txt$") |>
        str_remove("\\.(smap|md|recal)$")
}

## Project label is the directory enclosing the metrics tree, after stripping
## the conventional out/ and Map/ wrappers.
projectFromMetricsPath <- function(path) {
    str_remove(path, "/metrics/.*$") |>
        str_remove("/out$") |>
        str_remove("/Map$") |>
        path_file()
}

##
## Tumor/normal class inferred from the sample name. Naming is inconsistent
## across the archive (_N01, _N, -N, trailing N, optional _D suffix), so this
## returns "unknown" rather than guessing when no marker is present. Callers
## must treat "unknown" as a real category and not as normal.
##
classifySampleType <- function(sample) {
    case_when(
        str_detect(sample, "[._-]?N[0-9]*([._-]D[0-9]*)?$") ~ "N",
        str_detect(sample, "[._-]?T[0-9]*([._-]D[0-9]*)?$") ~ "T",
        .default = "unknown")
}

## Patient stem is the sample name with the tumor/normal token removed, used to
## group samples into pairs. APTL_MDA012_N02 and APTL_MDA012_T01 both reduce to
## APTL_MDA012.
patientStem <- function(sample) {
    str_remove(sample, "[._-]?[NT][0-9]*([._-]D[0-9]*)?$")
}

##
## Column accessor that tolerates schema drift.
##
## Picard's AlignmentSummaryMetrics schema changed between the version that
## produced the archive and the one running now: MEAN_ALIGNED_READ_LENGTH and
## the read length quantiles are absent from every historical file. Missing
## columns must surface as NA, which the filters then report as MISSING. Letting
## them error would make the importer refuse whole cohorts, and defaulting them
## to zero would silently fail every archived sample.
##
pickColumn <- function(dat, column) {
    if (column %in% names(dat)) as.numeric(dat[[column]][1]) else NA_real_
}

##
## Read one sample's Picard alignment summary metrics.
## Use the PAIR row. FIRST_OF_PAIR and SECOND_OF_PAIR each describe half the
## data and reading either one instead is a quiet way to halve every count.
##
readAsmMetrics <- function(path) {
    dat <- readPicardMetrics(path)
    if (is.null(dat) || !"CATEGORY" %in% names(dat)) return(NULL)
    row <- dat |> filter(CATEGORY == "PAIR")
    if (nrow(row) == 0) row <- dat |> slice(1)

    meanReadLength    <- pickColumn(row, "MEAN_READ_LENGTH")
    meanAlignedLength <- pickColumn(row, "MEAN_ALIGNED_READ_LENGTH")

    tibble(
        totalReads        = pickColumn(row, "TOTAL_READS"),
        pctChimeras       = pickColumn(row, "PCT_CHIMERAS") * 100,
        pctSoftclip       = pickColumn(row, "PCT_SOFTCLIP") * 100,
        pctAdapter        = pickColumn(row, "PCT_ADAPTER") * 100,
        pctImproperPairs  = pickColumn(row, "PCT_PF_READS_IMPROPER_PAIRS") * 100,
        pctReadsAligned   = pickColumn(row, "PCT_PF_READS_ALIGNED") * 100,
        meanReadLength    = meanReadLength,
        meanAlignedLength = meanAlignedLength,
        alignedFrac       = meanAlignedLength / meanReadLength,
        ## Expressed as a percentage so it reads on the same scale as every other
        ## number in the report. A bare 0.86 next to a column of percentages is
        ## read as 0.86%, which is the opposite of what it means.
        pctReadUsed       = meanAlignedLength / meanReadLength * 100,
        strandBalance     = pickColumn(row, "STRAND_BALANCE"))
}

## Read one sample's Picard WGS coverage metrics.
readWgsMetrics <- function(path) {
    dat <- readPicardMetrics(path)
    if (is.null(dat) || !"MEAN_COVERAGE" %in% names(dat)) return(NULL)
    row <- dat |> slice(1)

    tibble(
        meanCoverage   = pickColumn(row, "MEAN_COVERAGE"),
        medianCoverage = pickColumn(row, "MEDIAN_COVERAGE"),
        sdCoverage     = pickColumn(row, "SD_COVERAGE"),
        pctExcOverlap  = pickColumn(row, "PCT_EXC_OVERLAP") * 100,
        pctExcDupe     = pickColumn(row, "PCT_EXC_DUPE") * 100,
        pctExcMapq     = pickColumn(row, "PCT_EXC_MAPQ") * 100,
        pctExcTotal    = pickColumn(row, "PCT_EXC_TOTAL") * 100,
        pct30x         = pickColumn(row, "PCT_30X") * 100)
}

##
## Read the multiqc samtools stats table, keeping only sample level rows.
##
## multiqc mixes granularities in one file. The general stats table is per
## read group (SAMPLE-FLOWCELL_LANE), while samtools stats carries clean
## per-sample rows suffixed .md and .recal. Take .recal: it is the final BAM
## state and the same object Picard measured, so the two sources are comparable.
##
readMultiqcSamtools <- function(path) {
    if (!file_exists(path)) return(NULL)
    dat <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
    wanted <- c("supplementary_alignments", "raw_total_sequences",
                "insert_size_average", "reads_properly_paired_percent",
                "pairs_on_different_chromosomes", "reads_mapped_percent")
    if (!all(wanted %in% names(dat))) return(NULL)
    dat |>
        filter(str_detect(Sample, "\\.recal$")) |>
        transmute(
            sample             = str_remove(Sample, "\\.recal$"),
            supplementaryRate  = supplementary_alignments / raw_total_sequences * 100,
            interChromRate     = pairs_on_different_chromosomes / raw_total_sequences * 100,
            insertSizeAverage  = insert_size_average,
            pctProperlyPaired  = reads_properly_paired_percent,
            pctMapped          = reads_mapped_percent)
}

##
## Collect every sample under a metrics tree into one row per sample.
## root may be a Map directory or any ancestor; all metrics/ trees below it are
## picked up, which is what lets the importer walk the whole archive at once.
##
collectPicardSamples <- function(root) {
    asmFiles <- dir_ls(root, recurse = TRUE, glob = "*.asm.txt", fail = FALSE)
    wgsFiles <- dir_ls(root, recurse = TRUE, glob = "*.wgs.txt", fail = FALSE)
    if (length(asmFiles) == 0 && length(wgsFiles) == 0) return(NULL)

    asm <- tibble(path = as.character(asmFiles)) |>
        mutate(sample = sampleFromMetricsPath(path),
               project = projectFromMetricsPath(path),
               parsed = map(path, readAsmMetrics),
               hasAsm = !map_lgl(parsed, is.null)) |>
        filter(hasAsm) |>
        select(project, sample, asmPath = path, parsed) |>
        unnest(parsed)

    wgs <- tibble(path = as.character(wgsFiles)) |>
        mutate(sample = sampleFromMetricsPath(path),
               project = projectFromMetricsPath(path),
               parsed = map(path, readWgsMetrics),
               hasWgs = !map_lgl(parsed, is.null)) |>
        filter(hasWgs) |>
        select(project, sample, wgsPath = path, parsed) |>
        unnest(parsed)

    full_join(asm, wgs, by = c("project", "sample")) |>
        mutate(sampleType = classifySampleType(sample),
               patient = patientStem(sample)) |>
        relocate(project, sample, sampleType, patient)
}

##
## Evaluate the sample level filter thresholds. Returns one row per sample per metric so
## that the console, HTML and TSV outputs all render from the same evaluation
## rather than each recomputing its own verdict.
##
evaluateThresholds <- function(dat) {
    present <- THRESHOLDS |> filter(metric %in% names(dat))

    dat |>
        select(sample, all_of(present$metric)) |>
        pivot_longer(-sample, names_to = "metric", values_to = "value") |>
        left_join(present, by = "metric") |>
        mutate(
            status = case_when(
                is.na(value)                            ~ "MISSING",
                direction == "high" & !is.na(fail) & value > fail ~ "FAIL",
                direction == "low"  & !is.na(fail) & value < fail ~ "FAIL",
                direction == "high" & !is.na(warn) & value > warn ~ "WARN",
                direction == "low"  & !is.na(warn) & value < warn ~ "WARN",
                .default = "PASS"))
}

##
## Roll per-metric statuses up to one verdict per sample.
## MISSING is not PASS. A sample whose metrics could not be read has not been
## checked, and reporting it as clean is how a partial report becomes a wrong
## report.
##
sampleVerdict <- function(thresholdResults) {
    thresholdResults |>
        summarise(
            nFail = sum(status == "FAIL"),
            nWarn = sum(status == "WARN"),
            nMissing = sum(status == "MISSING"),
            failedMetrics = str_c(metric[status == "FAIL"], collapse = ","),
            warnedMetrics = str_c(metric[status == "WARN"], collapse = ","),
            .by = sample) |>
        mutate(
            escalated = nFail == 0 & nWarn >= WARN_ESCALATION,
            verdict = case_when(
                nFail > 0    ~ "FAIL",
                escalated    ~ "FAIL",
                nMissing > 0 ~ "INCOMPLETE",
                nWarn > 0    ~ "WARN",
                .default     = "PASS"),
            verdictReason = case_when(
                nFail > 0 ~ glue("failed {nFail} threshold(s): {failedMetrics}"),
                escalated ~ glue("{nWarn} simultaneous warnings: {warnedMetrics}"),
                nMissing > 0 ~ glue("{nMissing} metric(s) could not be read"),
                nWarn > 0 ~ glue("{nWarn} warning(s): {warnedMetrics}"),
                .default = "all thresholds within range"))
}

##
## Pair level checks. Two samples can each pass individually and still be
## unusable together, so these are computed after the sample level verdicts and
## reported alongside them rather than folded into them.
##
evaluatePairs <- function(dat, verdicts) {
    paired <- dat |>
        select(sample, patient, sampleType, insertSizeAverage) |>
        left_join(verdicts |> select(sample, verdict), by = "sample") |>
        filter(sampleType %in% c("N", "T"))

    tumors  <- paired |> filter(sampleType == "T")
    normals <- paired |> filter(sampleType == "N")

    inner_join(tumors, normals, by = "patient", suffix = c("Tumor", "Normal")) |>
        mutate(
            insertRatio = pmax(insertSizeAverageTumor, insertSizeAverageNormal) /
                          pmin(insertSizeAverageTumor, insertSizeAverageNormal),
            insertStatus = case_when(
                is.na(insertRatio)                ~ "MISSING",
                insertRatio > INSERT_RATIO_FAIL   ~ "FAIL",
                .default                          = "PASS"),
            pairVerdict = case_when(
                verdictTumor == "FAIL" | verdictNormal == "FAIL" ~ "FAIL",
                insertStatus == "FAIL"                           ~ "FAIL",
                verdictTumor == "INCOMPLETE" | verdictNormal == "INCOMPLETE" ~ "INCOMPLETE",
                .default                                          = "PASS"),
            pairReason = case_when(
                verdictTumor == "FAIL" & verdictNormal == "FAIL" ~ "both samples failed QC",
                verdictNormal == "FAIL"                          ~ "normal failed QC",
                verdictTumor == "FAIL"                           ~ "tumor failed QC",
                insertStatus == "FAIL" ~ glue("insert size differs {round(insertRatio, 2)}x"),
                .default                                         = ""))
}

