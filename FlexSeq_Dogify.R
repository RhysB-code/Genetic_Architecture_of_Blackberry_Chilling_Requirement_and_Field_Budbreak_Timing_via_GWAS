# ============================================================
# UPDOG RERUN v2 -- excluding 20 samples with >20% missing
# genotype calls (DP-based), per Mason's recommendation that
# these be removed BEFORE calling, not after, since their noise
# affects population-wide marker parameter estimates.
# ============================================================

library(updog)
library(data.table)
library(tidyverse)

work_dir <- "~/blackberry_flexseq/vcf"
ad_file      <- file.path(work_dir, "ARK_133905.AD.FORMAT")
ref_alt_file <- file.path(work_dir, "ARK_133905_ref_alt.txt")

# ---- The 20 excluded samples (RG_Sample_Code, as they appear
#      as column headers in the AD.FORMAT file) ----
excluded_samples <- c(
  "ARK_133905_P001_WA03",  # Immaculate    79.19%
  "ARK_133905_P002_WA03",  # APF-751T      76.40%
  "ARK_133905_P001_WG08",  # APF-722T      61.80%
  "ARK_133905_P001_WB09",  # A-2987T       60.20%
  "ARK_133905_P001_WD12",  # A-3008T       51.87%
  "ARK_133905_P002_WC11",  # APF-732TN     51.07%
  "ARK_133905_P001_WF11",  # APF-711T      49.91%
  "ARK_133905_P002_WF03",  # APF-789T      44.65%
  "ARK_133905_P002_WE10",  # APF-784T      40.13%
  "ARK_133905_P003_WA03",  # A-3033T       39.46%
  "ARK_133905_P002_WF07",  # APF-793T      36.80%
  "ARK_133905_P001_WD03",  # A-2491T       35.67%
  "ARK_133905_P003_WC03",  # A-3060TN      35.51%
  "ARK_133905_P001_WC01",  # A-2991T       31.76%
  "ARK_133905_P003_WC11",  # A-3068T       30.75%
  "ARK_133905_P001_WB01",  # A-2976T       30.02%
  "ARK_133905_P002_WH12",  # A-3029T       28.67%
  "ARK_133905_P001_WG11",  # APF-726T      25.07%
  "ARK_133905_P001_WC06",  # A-2993T       22.78%
  "ARK_133905_P002_WD02"   # APF-748T      20.10%
)

cat('Excluding', length(excluded_samples), 'samples (>20% missing, DP-based)\n\n')

ad_raw <- read.table(ad_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE)

all_samples <- colnames(ad_raw)[-(1:2)]
cat('Total samples before exclusion:', length(all_samples), '\n')

# ---- Confirm every excluded name actually exists in the file
#      before dropping anything -- catches a typo immediately ----
missing_from_file <- setdiff(excluded_samples, all_samples)
if (length(missing_from_file) > 0) {
  stop("These excluded sample names were not found in the AD.FORMAT file: ",
       paste(missing_from_file, collapse = ", "))
}

samples <- setdiff(all_samples, excluded_samples)
n_samples <- length(samples)
cat('Samples remaining after exclusion:', n_samples, '(expect 210)\n\n')

# ---- Parse "ref,alt" text into two numeric matrices, using only
#      the retained samples ----
ref_mat <- matrix(NA_real_, nrow = nrow(ad_raw), ncol = n_samples)
siz_mat <- matrix(NA_real_, nrow = nrow(ad_raw), ncol = n_samples)
for (j in seq_len(n_samples)) {
  split_vals <- strsplit(as.character(ad_raw[[samples[j]]]), ',')
  ref_mat[, j] <- sapply(split_vals, function(x) if (length(x) < 2 || any(x == '.')) NA_real_ else as.numeric(x[1]))
  siz_mat[, j] <- ref_mat[, j] + sapply(split_vals, function(x) if (length(x) < 2 || any(x == '.')) NA_real_ else as.numeric(x[2]))
}
marker_keys <- paste0(ad_raw$CHROM, ':', ad_raw$POS)
rownames(ref_mat) <- marker_keys; rownames(siz_mat) <- marker_keys
colnames(ref_mat) <- samples;     colnames(siz_mat) <- samples

# ---- 80% observation filter (marker-level), recomputed fresh
#      on the reduced 210-sample set ----
keep <- rowMeans(!is.na(siz_mat)) >= 0.8
ref_mat <- ref_mat[keep, ]
siz_mat <- siz_mat[keep, ]
cat('[', format(Sys.time()), '] Markers after 80% obs filter (210 samples):', nrow(ref_mat), '\n')

# ---- BATCHED multidog(), same resumable design as before ----
n_cores <- parallel::detectCores()
cat('Detected cores:', n_cores, '\n')

batch_size <- 5000
n_markers  <- nrow(ref_mat)
n_batches  <- ceiling(n_markers / batch_size)
batch_dir  <- file.path(work_dir, "updog_batches_v2_excl20")
dir.create(batch_dir, showWarnings = FALSE)

cat('Total markers:', n_markers, '| Batches:', n_batches, '| Cores:', n_cores, '\n\n')

for (b in seq_len(n_batches)) {
  batch_file <- file.path(batch_dir, paste0("batch_", b, ".rds"))
  if (file.exists(batch_file)) {
    cat('[', format(Sys.time()), '] Batch', b, 'of', n_batches, 'already done -- skipping\n')
    next
  }
  i_start <- (b - 1) * batch_size + 1
  i_end   <- min(b * batch_size, n_markers)
  cat('[', format(Sys.time()), '] Starting batch', b, 'of', n_batches, '(markers', i_start, '-', i_end, ')\n')

  t0 <- Sys.time()
  batch_result <- multidog(
    refmat  = ref_mat[i_start:i_end, , drop = FALSE],
    sizemat = siz_mat[i_start:i_end, , drop = FALSE],
    ploidy  = 4,
    model   = 'norm',
    nc      = n_cores
  )
  t1 <- Sys.time()
  saveRDS(batch_result, batch_file)
  cat('[', format(Sys.time()), '] Batch', b, 'done in', round(difftime(t1, t0, units = "mins"), 1), 'min. Saved.\n')
}

cat('\n[ALL BATCHES DONE -- ', format(Sys.time()), ']\n')

# ============================================================
# COMBINE BATCHES (numeric sort -- NOT alphabetical) + finalize
# ============================================================
batch_files <- list.files(batch_dir, pattern = "batch_.*\\.rds", full.names = TRUE)
batch_files_sorted <- batch_files[order(as.numeric(gsub("[^0-9]", "", basename(batch_files))))]

all_snp_stats <- list(); all_postmean <- list()
for (bf in batch_files_sorted) {
  bres <- readRDS(bf)
  all_snp_stats[[bf]] <- bres$snpdf
  all_postmean[[bf]]  <- format_multidog(bres, varname = 'postmean')
}
snp_stats    <- dplyr::bind_rows(all_snp_stats)
postmean_mat <- do.call(rbind, all_postmean)
cat('Combined:', nrow(postmean_mat), 'markers | rownames intact:', !is.null(rownames(postmean_mat)), '\n')

# ---- Quality filter: prop_mis <= 0.05 (correct column name) ----
keep_snps <- snp_stats$snp[!is.na(snp_stats$prop_mis) & snp_stats$prop_mis <= 0.05]
postmean_mat <- postmean_mat[keep_snps, , drop = FALSE]
cat('After prop_mis <= 0.05 filter:', nrow(postmean_mat), 'markers\n')

# ---- Rebuild full marker IDs with REF/ALT alleles ----
ref_alt <- read.table(ref_alt_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                       col.names = c("CHROM", "POS", "REF", "ALT"))
ref_alt$key <- paste0(ref_alt$CHROM, ':', ref_alt$POS)
ref_alt_dedup <- ref_alt[!duplicated(ref_alt$key), ]

marker_keys2 <- rownames(postmean_mat)
split_keys   <- strsplit(marker_keys2, ':', fixed = TRUE)
chrom_vec    <- sapply(split_keys, `[`, 1)
pos_vec      <- as.integer(sapply(split_keys, `[`, 2))

marker_info <- data.frame(key = marker_keys2, Chrom = chrom_vec, Position = pos_vec, stringsAsFactors = FALSE)
marker_info <- dplyr::left_join(marker_info, ref_alt_dedup[, c("key", "REF", "ALT")], by = "key")
marker_info$Marker <- paste0(marker_info$Chrom, ':', marker_info$Position, '_', marker_info$REF, '/', marker_info$ALT)
marker_info <- marker_info[, c("Marker", "Chrom", "Position")]

gwaspy_df <- cbind(marker_info, as.data.frame(postmean_mat))

dosage_out <- file.path(work_dir, "FlexSeq_2026_ARK133905_TetraploidPostmean_GWASready_v2_excl20.csv")
stats_out  <- file.path(work_dir, "FlexSeq_2026_ARK133905_UpdogQualityStats_v2_excl20.csv")

data.table::fwrite(gwaspy_df, dosage_out, row.names = FALSE)
data.table::fwrite(snp_stats, stats_out, row.names = FALSE)

cat('\n[FINAL DONE] Markers:', nrow(gwaspy_df), '| Samples:', ncol(gwaspy_df) - 3, '(expect 210)\n')
cat('Marker ID example:', gwaspy_df$Marker[1], '\n')
cat('Saved:', dosage_out, '\n')
