#!/usr/bin/env Rscript
# Download the per-library, donor-assigned snRNA-seq count matrices (87 files, ~41 GB)
# from Ling et al. (2024) Nature 627:604-611, NeMO nemo:dat-djf3v1i (processed/counts/).
# Usage: Rscript download_ling2024_counts.R [output_dir]   (default: ./ling2024_counts)
# Re-running skips files that are already complete.

base_url <- "https://data.nemoarchive.org/other/grant/broad/mccarroll/transcriptome/sncell/10x_v3.1/human/processed/counts/"
args     <- commandArgs(trailingOnly = TRUE)
out_dir  <- if (length(args)) args[1] else "ling2024_counts"

files <- c(
  "BA46_05-07-2019_rxn2_MyelinRemoved.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-07-2019_rxn3_MyelinRemoved.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-07-2019_rxn4_MyelinRemoved.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-07-2019_rxn5_MyelinRemoved.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-07-2019_rxn6_MyelinRemoved.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-07-2019_rxn7_MyelinRemoved.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-07-2019_rxn8_MyelinRemoved.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-21-2021_rxn1.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-21-2021_rxn2.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-21-2021_rxn3.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-21-2021_rxn4.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-21-2021_rxn5.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-21-2021_rxn6.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-21-2021_rxn7.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-21-2021_rxn8.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-28-2021_rxn1.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-28-2021_rxn2.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-28-2021_rxn3.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-28-2021_rxn4.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-28-2021_rxn5.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-28-2021_rxn6.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-28-2021_rxn7.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_05-28-2021_rxn8.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_07-10-2019_rxn1.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_07-10-2019_rxn2.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_07-10-2019_rxn3.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_07-10-2019_rxn4.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_07-10-2019_rxn5.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_07-10-2019_rxn6.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_07-10-2019_rxn7.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_07-10-2019_rxn8.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-09-2019_rxn1.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-09-2019_rxn2.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-09-2019_rxn3.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-09-2019_rxn4.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-09-2019_rxn5.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-09-2019_rxn6.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-09-2019_rxn7.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-09-2019_rxn8.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-27-2019_rxn1.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-27-2019_rxn2.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-27-2019_rxn3.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-27-2019_rxn4.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-27-2019_rxn5.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-27-2019_rxn6.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-27-2019_rxn7.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_08-27-2019_rxn8.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-02-2021_rxn1.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-02-2021_rxn2.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-02-2021_rxn3.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-02-2021_rxn4.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-02-2021_rxn5.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-02-2021_rxn6.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-02-2021_rxn7.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-02-2021_rxn8.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-03-2021_rxn1.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-03-2021_rxn2.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-03-2021_rxn3.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-03-2021_rxn4.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-03-2021_rxn5.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-03-2021_rxn6.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-03-2021_rxn7.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-03-2021_rxn8.intronic_0.300-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-18-2019_rxn1.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-18-2019_rxn2.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-18-2019_rxn3.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-18-2019_rxn4.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-18-2019_rxn5.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-18-2019_rxn6.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-18-2019_rxn7.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-18-2019_rxn8.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-22-2021_rxn1.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-22-2021_rxn2.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-22-2021_rxn3.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-22-2021_rxn4.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-22-2021_rxn5.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-22-2021_rxn6.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-22-2021_rxn7.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_09-22-2021_rxn8.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_10-16-2019_rxn1.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_10-16-2019_rxn2.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_10-16-2019_rxn3.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_10-16-2019_rxn4.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_10-16-2019_rxn5.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_10-16-2019_rxn6.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_10-16-2019_rxn7.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5",
  "BA46_10-16-2019_rxn8.intronic_0.200-1.000.10X.exonic+intronic.gmg.donors.h5"
)

options(timeout = 7200)  # large files; R's default 60 s is too short
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (i in seq_along(files)) {
  f    <- files[i]
  url  <- paste0(base_url, gsub("+", "%2B", f, fixed = TRUE))  # "+" must be URL-encoded
  dest <- file.path(out_dir, f)
  hdr  <- curlGetHeaders(url)
  size <- as.numeric(sub(".*:\\s*", "", tail(grep("^content-length", hdr, ignore.case = TRUE, value = TRUE), 1)))
  if (file.exists(dest) && file.size(dest) == size) { message(sprintf("[%d/%d] skip %s", i, length(files), f)); next }
  message(sprintf("[%d/%d] get  %s", i, length(files), f))
  part <- paste0(dest, ".part")
  ok <- tryCatch(download.file(url, part, mode = "wb", quiet = TRUE) == 0,
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (ok && file.size(part) == size) file.rename(part, dest)
  else { unlink(part); message("  FAILED (re-run the script to retry)") }
}