# ============================================================
# FlexSeq Tetraploid Dosage Recall — Complete Script
# Project ARK_133905
# ============================================================

# ---- 0. Install any missing packages ----
required_packages <- c("updog", "data.table", "tidyverse")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

library(updog)
library(data.table)
library(tidyverse)

# ---- 1. File locations (hardcoded, explicit — no reliance on getwd()) ----
input_dir  <- "C:/Users/RhysB/OneDrive/Desktop"
output_dir <- "C:/Users/RhysB/OneDrive/Desktop/CH2 GWAS"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ad_file      <- file.path(input_dir, "ARK_133905.AD.FORMAT")
ref_alt_file <- file.path(input_dir, "ARK_133905_ref_alt.txt")

# Confirm both input files actually exist before doing anything else
stopifnot(
  "AD.FORMAT file not found -- check the path" = file.exists(ad_file),
  "ref_alt.txt file not found -- check the path" = file.exists(ref_alt_file)
)
cat('[OK] Both input files found.\n')

# ---- 2. Load the AD matrix ----
ad_raw <- read.table(ad_file, header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, check.names = FALSE)
samples <- colnames(ad_raw)[-(1:2)]
n_samples <- length(samples)

cat('Markers:', nrow(ad_raw), '| Samples:', n_samples, '\n')

# ---- 3. Parse "ref,alt" text into two numeric matrices ----
ref_mat <- matrix(NA_real_, nrow = nrow(ad_raw), ncol = n_samples)
siz_mat <- matrix(NA_real_, nrow = nrow(ad_raw), ncol = n_samples)

for (j in seq_len(n_samples)) {
  col        <- ad_raw[[samples[j]]]
  split_vals <- strsplit(as.character(col), ',')
  ref_counts <- sapply(split_vals, function(x) {
    if (length(x) < 2 || any(x == '.')) return(NA_real_)
    as.numeric(x[1])
  })
  alt_counts <- sapply(split_vals, function(x) {
    if (length(x) < 2 || any(x == '.')) return(NA_real_)
    as.numeric(x[2])
  })
  ref_mat[, j] <- ref_counts
  siz_mat[, j] <- ref_counts + alt_counts
}

marker_keys <- paste0(ad_raw$CHROM, ':', ad_raw$POS)
rownames(ref_mat) <- marker_keys; rownames(siz_mat) <- marker_keys
colnames(ref_mat) <- samples;     colnames(siz_mat) <- samples

# ---- 4. Drop markers missing in more than 20% of samples ----
keep <- rowMeans(!is.na(siz_mat)) >= 0.8
ref_mat <- ref_mat[keep, ]
siz_mat <- siz_mat[keep, ]
cat('Markers passing 80% obs filter:', sum(keep), 'of', length(keep), '\n')

# ---- 5. Run updog's tetraploid dosage recall (the slow step) ----
cat('Detected CPU cores:', parallel::detectCores(), '\n')
cat('[UPDOG] Starting multidog() -- ploidy=4, model=norm...\n')
cat('[UPDOG] This can take tens of minutes to a few hours. Do not close RStudio.\n')

start_time <- Sys.time()

mout <- multidog(
  refmat  = ref_mat,
  sizemat = siz_mat,
  ploidy  = 4,
  model   = 'norm',
  nc      = 4   # adjust if parallel::detectCores() above is lower
)

end_time <- Sys.time()
cat('[UPDOG] Done! Took:', round(difftime(end_time, start_time, units = "mins"), 1), 'minutes\n')

# Save immediately -- this is the expensive part, protect it right away
saveRDS(mout, file.path(output_dir, "FlexSeq_2026_ARK133905_updog_raw_result.rds"))
cat('[OK] Raw updog result saved.\n')

# ---- 6. Quality filter: drop markers with >5% estimated error rate ----
snp_stats <- mout$snpdf
keep_snps <- snp_stats$snp[!is.na(snp_stats$oe_prop) & snp_stats$oe_prop <= 0.05]
cat('Markers passing error filter (<=5%):', length(keep_snps), 'of', nrow(snp_stats), '\n')

# ---- 7. Extract POSTERIOR MEAN dosage (continuous) -- matches CaptureSeq format ----
postmean_mat <- format_multidog(mout, varname = 'postmean')
postmean_mat <- postmean_mat[keep_snps, ]

# ---- 8. Rebuild full marker IDs using REF/ALT alleles ----
ref_alt <- read.table(ref_alt_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                      col.names = c("CHROM", "POS", "REF", "ALT"))
ref_alt$key <- paste0(ref_alt$CHROM, ':', ref_alt$POS)

marker_info <- data.frame(key = rownames(postmean_mat), stringsAsFactors = FALSE) %>%
  tidyr::separate(key, into = c('Chrom', 'Position'), sep = ':', remove = FALSE) %>%
  dplyr::mutate(Position = as.integer(Position)) %>%
  dplyr::left_join(ref_alt %>% dplyr::select(key, REF, ALT), by = "key") %>%
  dplyr::mutate(Marker = paste0(Chrom, ':', Position, '_', REF, '/', ALT)) %>%
  dplyr::select(Marker, Chrom, Position)

# ---- 9. Combine and save, matching CaptureSeq file's column structure ----
gwaspy_df <- dplyr::bind_cols(marker_info, as.data.frame(postmean_mat))

dosage_filename <- "FlexSeq_2026_ARK133905_TetraploidPostmean_GWASready.csv"
stats_filename  <- "FlexSeq_2026_ARK133905_UpdogQualityStats.csv"

data.table::fwrite(gwaspy_df, file.path(output_dir, dosage_filename), row.names = FALSE)
data.table::fwrite(snp_stats, file.path(output_dir, stats_filename), row.names = FALSE)

cat('\n[DONE]\n')
cat('  File:', dosage_filename, '\n')
cat('  Location:', output_dir, '\n')
cat('  Markers:', nrow(gwaspy_df), '| Samples:', ncol(gwaspy_df) - 3, '\n')
cat('  Marker ID example:', gwaspy_df$Marker[1], '\n')

