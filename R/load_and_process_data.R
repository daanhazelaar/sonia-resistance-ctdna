# =============================================================================
# Load and process data
# Loads ichorCNA (sWGS), TSO500 copy number, and TSO500 variant data.
# Produces: copy_number_ichorcna, copy_number_genes, copy_number_genes_wide,
#           somatic_tso500_ALL_trace_variants_paired,
#           somatic_tso500_PD_trace_variants_paired,
#           fastqseq_zscores
# =============================================================================

library(dplyr)
library(ggplot2)

source("./R/functions/import_ichorcna_data.R")

# --- Metadata ----------------------------------------------------------------

meta_data <- readxl::read_xlsx(
  "data/meta_data_SONIA_v6.xlsx",
  col_types = rep("text", 6)
)


n_patients_cohort <- length(unique(meta_data$patnr))
n_patients_arm_A  <- meta_data |> filter(arm == "arm_A") |> pull(patnr) |> unique() |> length()
n_patients_arm_B  <- meta_data |> filter(arm == "arm_B") |> pull(patnr) |> unique() |> length()

# --- ichorCNA sWGS copy number -----------------------------------------------

ichorcna_dirs <- list.files("data/ichorcna/", full.names = TRUE)

estimatesIchorCNA_sonia <- lapply(ichorcna_dirs, getEstimatesFromIchorCnaRdataFiles) |>
  bind_rows() |>
  mutate(bins = paste0(seqnames, "_", start)) |>
  dplyr::rename(sampleName_shWGS = sampleName) |>
  left_join(meta_data, by = join_by(sampleName_shWGS)) |>
  relocate(patnr)

# Tumor fraction estimates per sample, with manual entries for patients without sWGS
tf_estimates_tso500_sampleNames <- estimatesIchorCNA_sonia |>
  group_by(sampleName_shWGS) |>
  summarise(tumorFractionEstimate = unique(tumorFractionEstimate), .groups = "drop") |>
  right_join(meta_data, by = join_by(sampleName_shWGS)) |>
  mutate(
    sampleName_shWGS = case_when(
      sampleName_fastseq == "L4275"  ~ "L4275",
      sampleName_fastseq == "L10846" ~ "L10846",
      .default = sampleName_shWGS
    ),
    tumorFractionEstimate = case_when(
      sampleName_fastseq == "L4275"  ~ 0.141,
      sampleName_fastseq == "L10846" ~ 0.433,
      sampleName_fastseq == "L10176" ~  0.130,
      .default = tumorFractionEstimate
    )
  )

# --- Build 1 Mb genomic bins and fill ichorCNA CN across genome ---------------

chr_lengths <- c(
  249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
  159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
  115169878, 107349540, 102531392,  90354753,  81195210,  78077248,
   59128983,  63025520,  48129895,  51304566, 155270560,  59373566
)

create_bins <- function(seqnames, max_val, interval = 1e6) {
  breaks <- seq(1, max_val + interval, by = interval)
  tibble(
    seqnames = seqnames,
    start = breaks[-length(breaks)],
    end   = breaks[-1] - 1
  )
}

all_ichorCNA_bins <- mapply(
  create_bins,
  seqnames = paste0("chr", c(1:22, "X", "Y")),
  max_val  = chr_lengths,
  SIMPLIFY = FALSE
) |>
  bind_rows()

copy_number_ichorcna <- estimatesIchorCNA_sonia |>
  distinct(sampleName_shWGS) |>
  pull(sampleName_shWGS) |>
  lapply(\(sample) {
    sample_cn <- estimatesIchorCNA_sonia |>
      filter(sampleName_shWGS == .env$sample) |>
      select(patnr, sampleName_shWGS, seqnames, start, end, Corrected_Copy_Number)

    sample_cn |>
      right_join(all_ichorCNA_bins, by = join_by(seqnames, start, end)) |>
      arrange(seqnames, start) |>
      tidyr::fill(Corrected_Copy_Number, .direction = "downup") |>
      mutate(
        patnr = unique(sample_cn$patnr),
        sampleName_shWGS = .env$sample
      )
  }) |>
  bind_rows() |>
  select(sampleName_shWGS, seqnames, start, end, ichorCNA_cn = Corrected_Copy_Number)

# --- TSO500 gene-level copy number -------------------------------------------

estimate_gene_copy_number <- function(fold_change, tumor_fraction, normal_cn = 2) {
  if (tumor_fraction == 0) tumor_fraction <- 0.03
  (fold_change * normal_cn - (1 - tumor_fraction) * normal_cn) / tumor_fraction
}

copy_number_genes <- list.files("data/tso500_copy_number_files/", full.names = TRUE) |>
  lapply(\(path) {
    VariantAnnotation::readVcfAsVRanges(path) |>
      as_tibble() |>
      filter(!is.na(alt)) |>
      rename(sampleName_tso500 = sampleNames, genename = SEGID) |>
      mutate(SVLEN = unlist(SVLEN))
  }) |>
  bind_rows() |>
  left_join(tf_estimates_tso500_sampleNames, by = join_by(sampleName_tso500)) |>
  mutate(
    tumorFractionEstimate = round(tumorFractionEstimate, 4),
    positionComplete_copy_gene = paste0(seqnames, "_", start, "_", end)
  ) |>
  rowwise() |>
  mutate(
    cn_tumor_gene = estimate_gene_copy_number(SM, tumorFractionEstimate),
    cn_tumor_gene_int = case_when(
      cn_tumor_gene < 2 ~ pmax(1, round(cn_tumor_gene)),
      cn_tumor_gene > 2 ~ pmax(3, round(cn_tumor_gene))
    ),
    start_ichorcna_bin = floor(start / 1e6) * 1e6 + 1,
    end_ichorcna_bin   = start_ichorcna_bin + 1e6 - 1,
    alt = ifelse(alt == "<DUP>", "duplication", "deletion")
  ) |>
  ungroup() |>
  left_join(
    copy_number_ichorcna,
    by = join_by(sampleName_shWGS, seqnames, start_ichorcna_bin == start, end_ichorcna_bin == end)
  ) |>
  relocate(patnr)

# Classify focal CNA changes between baseline and progression
diff_variant_focal_cna <- copy_number_genes |>
  select(patnr, timepoint, genename, alt) |>
  tidyr::pivot_wider(names_from = timepoint, values_from = alt) |>
  mutate(
    diff_variant_focal_cna = case_when(
      BL == PD              ~ "both",
      !is.na(BL) & is.na(PD) ~ "BL",
      is.na(BL) & !is.na(PD) ~ "PD"
    )
  ) |>
  select(patnr, genename, diff_variant_focal_cna)

copy_number_genes <- copy_number_genes |>
  left_join(diff_variant_focal_cna, by = join_by(patnr, genename))

copy_number_genes_wide <- copy_number_genes |>
  select(patnr, arm, timepoint, genename, alt, diff_variant_focal_cna, cn_tumor_gene_int) |>
  tidyr::pivot_wider(names_from = timepoint, values_from = cn_tumor_gene_int, names_prefix = "cn_") |>
  tidyr::replace_na(list(cn_BL = 2, cn_PD = 2)) |>
  mutate(
    cn_category = case_when(
      cn_BL < 2 & cn_PD < 2                       ~ "stable_loss",
      cn_BL > 2 & cn_PD > 2                       ~ "stable_gain",
      cn_BL > 2 & cn_PD == 2                      ~ "gain_normalized",
      cn_BL < 2 & cn_PD == 2                      ~ "loss_normalized",
      cn_BL == 2 & cn_PD > 2                      ~ "acquired_gain",
      cn_BL == 2 & cn_PD < 2                      ~ "acquired_loss",
      cn_BL < 2 & cn_PD > 2                       ~ "loss_gain",
      cn_BL > 2 & cn_PD < 2                       ~ "gain_loss"
    )
  )

# --- TSO500 variant calling (trace files) ------------------------------------

trace_variants <- list.files("data/tso500_variant_trace_files/", full.names = TRUE) |>
  lapply(\(path) {
    readr::read_delim(path, delim = "\t", show_col_types = FALSE) |>
      rename_all(tolower) |>
      mutate(
        sampleName_tso500 = stringr::str_extract(path, "L[0-9]{4,5}"),
        positionComplete  = paste0(chromosome, "_", position)
      ) |>
      left_join(meta_data, by = join_by(sampleName_tso500)) |>
      mutate(
        sampleName_shWGS = case_when(
          sampleName_fastseq == "L4275"  ~ "L4275",
          sampleName_fastseq == "L10846" ~ "L10846",
          .default = sampleName_shWGS
        )
      ) |>
      relocate(patnr, sampleName_tso500, timepoint, positionComplete)
  }) |>
  bind_rows()

# --- Pair baseline/progression variants per patient --------------------------

trace_variants_paired <- unique(trace_variants$patnr) |>
  lapply(\(pat) {
    diff_variants <- trace_variants |>
      filter(patnr == .env$pat) |>
      select(patnr, arm, timepoint, positionComplete, genename, vaf) |>
      tidyr::pivot_wider(names_from = timepoint, values_from = vaf, names_prefix = "vaf_") |>
      left_join(copy_number_genes_wide, by = join_by(patnr, arm, genename)) |>
      tidyr::replace_na(list(cn_BL = 2, cn_PD = 2))

    pd_only <- diff_variants |> filter(is.na(vaf_BL)) |> pull(positionComplete)
    bl_only <- diff_variants |> filter(is.na(vaf_PD)) |> pull(positionComplete)

    trace_variants |>
      filter(patnr == .env$pat) |>
      mutate(
        diff_variant = case_when(
          positionComplete %in% pd_only ~ "PD",
          positionComplete %in% bl_only ~ "BL",
          .default = "both"
        )
      ) |>
      left_join(diff_variants, by = join_by(patnr, arm, positionComplete, genename)) |>
      filter(!is.na(genename))
  }) |>
  bind_rows()

# --- Resolve multi-gene annotations ------------------------------------------

all_genenames <- unique(trace_variants_paired$genename)

genename_lookup <- trace_variants_paired |>
  filter(stringr::str_detect(genename, ";")) |>
  distinct(genename) |>
  rename(genename_double = genename) |>
  tidyr::separate(genename_double, into = c("gene1", "gene2", "gene3"), sep = ";|-",
                  remove = FALSE) |>
  tidyr::pivot_longer(!genename_double, names_to = "slot", values_to = "genename_replace") |>
  filter(!is.na(genename_replace)) |>
  mutate(present = genename_replace %in% all_genenames) |>
  group_by(genename_double) |>
  filter(any(present) & present | (!any(present) & slot == "gene1")) |>
  ungroup() |>
  distinct(genename_double, genename_replace)

trace_variants_paired <- trace_variants_paired |>
  left_join(genename_lookup, by = join_by(genename == genename_double)) |>
  mutate(genename = coalesce(genename_replace, genename)) |>
  select(!genename_replace)

# --- Cancer cell fraction ----------------------------------------------------

estimate_ccf <- function(vaf, cn, tf) {
  if (tf <= 0 || tf > 1) stop("Tumor fraction must be between 0 and 1.")
  if (cn <= 0) stop("Copy number must be > 0.")

  denom <- (1 - tf) * 2 + tf * cn

  if (cn <= 2) {
    ccf <- (vaf * denom) / tf
  } else {
    ccf <- mean(sapply(1:cn, \(m) (vaf * denom) / (m * tf)))
  }

  pmin(pmax(ccf, 0), 1)
}

trace_variants_paired <- trace_variants_paired |>
  left_join(
    tf_estimates_tso500_sampleNames,
    by = join_by(patnr, sampleName_tso500, timepoint, arm, sampleName_shWGS, sampleName_fastseq)
  ) |>
  rowwise() |>
  mutate(
    ccf_BL = estimate_ccf(vaf_BL, cn_BL, tumorFractionEstimate),
    ccf_PD = estimate_ccf(vaf_PD, cn_PD, tumorFractionEstimate)
  ) |>
  ungroup()

# --- Filter somatic, protein-altering variants -------------------------------

protein_coding_consequences <- c(
  "missense_variant", "frameshift_variant", "stop_gained", "stop_lost",
  "start_lost", "inframe_deletion", "inframe_insertion", "protein_altering_variant",
  "splice_acceptor_variant", "splice_donor_variant", "stop_retained_variant",
  "inframe_deletion;splice_region_variant", "frameshift_variant;stop_gained",
  "missense_variant;splice_region_variant", "frameshift_variant;splice_region_variant",
  "inframe_insertion;stop_gained", "splice_region_variant;stop_gained",
  "coding_sequence_variant;splice_donor_variant",
  "coding_sequence_variant;intron_variant;splice_acceptor_variant",
  "coding_sequence_variant;intron_variant;splice_donor_variant",
  "5_prime_UTR_variant;coding_sequence_variant",
  "intron_variant;splice_acceptor_variant",
  "intron_variant;splice_donor_variant"
)

select_somatic_variants <- function(data, diff_filter, ccf_cols) {
  data |>
    filter(diff_variant %in% diff_filter, status == "Somatic",
           consequence %in% protein_coding_consequences) |>
    select(patnr, sampleName_tso500, timepoint, arm, genename, diff_variant,
           status, positionComplete, position, proteinchange, varianttype,
           refcall, altcall, consequence, vaf_BL, vaf_PD, cn_BL, cn_PD,
           ccf_BL, ccf_PD, depth, tumorFractionEstimate) |>
    filter(if_all(all_of(ccf_cols), \(x) x > 0.03))
}

somatic_BL   <- select_somatic_variants(trace_variants_paired, "BL", "ccf_BL")
somatic_PD   <- select_somatic_variants(trace_variants_paired, "PD", "ccf_PD")
somatic_both <- trace_variants_paired |>
  filter(diff_variant == "both", status == "Somatic",
         consequence %in% protein_coding_consequences, timepoint == "PD") |>
  select(patnr, sampleName_tso500, timepoint, arm, genename, diff_variant,
         status, positionComplete, position, proteinchange, varianttype,
         refcall, altcall, consequence, vaf_BL, vaf_PD, cn_BL, cn_PD,
         ccf_BL, ccf_PD, depth, tumorFractionEstimate) |>
  filter(ccf_BL > 0.03, ccf_PD > 0.03)

somatic_tso500_ALL_trace_variants_paired <- bind_rows(somatic_BL, somatic_PD, somatic_both) |>
  arrange(patnr, genename, positionComplete) |>
  mutate(mutation_id = paste0(genename, "_", positionComplete))

somatic_tso500_PD_trace_variants_paired <- bind_rows(somatic_PD, somatic_both) |>
  arrange(patnr, genename, positionComplete) |>
  mutate(mutation_id = paste0(genename, "_", positionComplete))


# --- Data for Griffin Analysis -------------------------------

# samples with ER- prediction
outlier_samples <- c("L8788", "L9747")

# Load GC bias profiles
gc_bias_files <- list.files("data/griffin_gc_bias/", pattern = "\\.GC_bias\\.txt$", full.names = TRUE)

gc_bias_all <- lapply(gc_bias_files, function(f) {
    sample_id <- stringr::str_extract(basename(f), "^[^.]+")
    readr::read_tsv(f, show_col_types = FALSE) |>
        mutate(sample = sample_id)
}) |>
    bind_rows()


# Filter to the fragment size range used in nucleosome profiling (100-200 bp)
gc_bias_relevant <- gc_bias_all |>
    filter(length >= 100, length <= 200, number_of_fragments > 0) |>
    mutate(is_outlier = sample %in% outlier_samples)



# Load Griffin predictions for annotation
griffin_predictions <- readr::read_tsv(
    "data/er_status_predictions.tsv",
    show_col_types = FALSE
) |>
    mutate(
        ER_prediction_label = ifelse(ER_prediction == 1, "ER+", "ER-"),
        patnr = as.character(patnr),
        is_outlier = sample %in% outlier_samples
    )


# Load uncorrected and GC-corrected coverage profiles (step 2)
coverage_dir <- "data/griffin_nucleosome_profiling/"
sample_dirs <- list.files(coverage_dir)

load_coverage <- function(sample_dirs, correction_type) {
    suffix <- ifelse(correction_type == "GC_corrected",
                     ".GC_corrected.coverage.tsv",
                     ".uncorrected.coverage.tsv")
    lapply(sample_dirs, function(s) {
        fpath <- file.path(coverage_dir, s, paste0(s, suffix))
        if (!file.exists(fpath)) return(NULL)
        df <- readr::read_tsv(fpath, show_col_types = FALSE)
        pos_cols <- grep("^-?[0-9]+$", colnames(df), value = TRUE)
        df |>
            select(site_name, all_of(pos_cols)) |>
            tidyr::pivot_longer(cols = all_of(pos_cols), names_to = "position", values_to = "coverage") |>
            mutate(position = as.integer(position), sample = s, correction = correction_type)
    }) |>
        bind_rows()
}

coverage_uncorrected <- load_coverage(sample_dirs, "uncorrected")
coverage_corrected <- load_coverage(sample_dirs, "GC_corrected")
coverage_both <- bind_rows(coverage_uncorrected, coverage_corrected)


# Load per-sample summary features from both correction types
load_summary_features <- function(sample_dirs, correction_type) {
    suffix <- ifelse(correction_type == "GC_corrected",
                     ".GC_corrected.coverage.tsv",
                     ".uncorrected.coverage.tsv")
    lapply(sample_dirs, function(s) {
        fpath <- file.path(coverage_dir, s, paste0(s, suffix))
        if (!file.exists(fpath)) return(NULL)
        df <- readr::read_tsv(fpath, show_col_types = FALSE)
        df |>
            select(mean_coverage, central_coverage, amplitude,
                   mean_reads_per_bp_in_normalization_window,
                   mean_reads_per_bp_in_saved_window,
                   site_name, sample) |>
            mutate(correction = correction_type)
    }) |>
        bind_rows()
}

features_uncorrected <- load_summary_features(sample_dirs, "uncorrected")
features_corrected <- load_summary_features(sample_dirs, "GC_corrected")
features_both <- bind_rows(features_uncorrected, features_corrected) |>
    mutate(
        site_name_clean = gsub("\\.5e-4_qval", "", site_name),
        is_outlier = sample %in% outlier_samples
    )

# --- mFast-SeqS Z-scores ------------------------------------------------------

fastqseq_zscores <- readxl::read_excel(
    "data/fastqseq_zscores.xlsx",
    col_types = c("text", "text", "text", "text", "numeric")
)

# --- Emerged chromosomal CNAs (manually curated) ------------------------------

emerged_chromosomal_cnas <- readxl::read_excel(
    "data/emerged_chromosomal_cnas.xlsx",
    col_types = c(rep("text", 5), "logical")
) |>
    filter(!is.na(cna_type))

# --- ESR1 mutation annotations ------------------------------------------------

esr1_annotations <- readr::read_csv("data/patients_with_esr1_mutations_annotated.csv",
                                    col_types = "cccclc") |>
    mutate(patnr = as.character(patnr)) |>
    select(patnr, timepoint, esr1_postive) |>
    rename(esr1_positive = esr1_postive) |>
    left_join(meta_data, by = join_by(patnr, timepoint))

# --- Tumor fraction estimates (simplified) ------------------------------------

tf_estimates <- estimatesIchorCNA_sonia |>
    group_by(sampleName_shWGS) |>
    summarise(tumor_fraction = unique(tumorFractionEstimate), .groups = "drop")

# --- Merged Griffin data frame ------------------------------------------------

griffin_data <- griffin_predictions |>
    rename(sampleName_shWGS = sample) |>
    mutate(patnr = as.character(patnr)) |>
    left_join(
        esr1_annotations |> select(patnr, timepoint, esr1_positive),
        by = c("patnr", "timepoint")
    ) |>
    left_join(tf_estimates, by = "sampleName_shWGS") |>
    mutate(
        esr1_positive = ifelse(is.na(esr1_positive), FALSE, esr1_positive),
        is_outlier = sampleName_shWGS %in% outlier_samples
    )

# --- Annotated coverage profiles ----------------------------------------------

coverage_annotated <- coverage_corrected |>
    left_join(
        griffin_data |> select(sampleName_shWGS, patnr, arm, timepoint, ER_prediction_label, is_outlier),
        by = c("sample" = "sampleName_shWGS")
    ) |>
    mutate(site_name_clean = gsub("\\.5e-4_qval", "", site_name))

# --- Griffin feature definitions ----------------------------------------------

feature_cols <- c(
    "central_coverage_ER_neg_heme.5e-4_qval",
    "mean_coverage_ER_neg_heme.5e-4_qval",
    "amplitude_ER_neg_heme.5e-4_qval",
    "central_coverage_ER_pos_heme.5e-4_qval",
    "mean_coverage_ER_pos_heme.5e-4_qval",
    "amplitude_ER_pos_heme.5e-4_qval",
    "central_coverage_ER_pos_specific.5e-4_qval",
    "mean_coverage_ER_pos_specific.5e-4_qval",
    "amplitude_ER_pos_specific.5e-4_qval",
    "central_coverage_ER_neg_specific.5e-4_qval",
    "mean_coverage_ER_neg_specific.5e-4_qval",
    "amplitude_ER_neg_specific.5e-4_qval"
)

feature_labels <- c(
    "Central cov (ER- heme)", "Mean cov (ER- heme)", "Amplitude (ER- heme)",
    "Central cov (ER+ heme)", "Mean cov (ER+ heme)", "Amplitude (ER+ heme)",
    "Central cov (ER+ specific)", "Mean cov (ER+ specific)", "Amplitude (ER+ specific)",
    "Central cov (ER- specific)", "Mean cov (ER- specific)", "Amplitude (ER- specific)"
)
names(feature_labels) <- feature_cols

# --- Model coefficients -------------------------------------------------------

model_coefs <- readr::read_tsv("data/model_coefficients.tsv",
                               show_col_types = FALSE)

# --- Shared analysis objects (used by compare_arms and compare_within_arm) ----

fisher_gene <- function(x, nA_total = 21, nB_total = 14) {
    n_mut_A <- x$arm_A
    n_mut_B <- x$arm_B
    cont_table <- matrix(
        c(n_mut_A, nA_total - n_mut_A, n_mut_B, nB_total - n_mut_B),
        nrow = 2, byrow = TRUE
    )
    test <- fisher.test(cont_table)
    tibble::tibble(
        p_value     = test$p.value,
        odds_ratio  = unname(test$estimate),
        conf_low    = test$conf.int[1],
        conf_high   = test$conf.int[2],
        n_mut_A     = n_mut_A,
        n_not_mut_A = nA_total - n_mut_A,
        n_mut_B     = n_mut_B,
        n_not_mut_B = nB_total - n_mut_B
    )
}

top_genes <- c("ESR1", "KMT2A", "CDH1", "RB1", "SPTA1", "PIK3CA", "PRKDC", "GATA3", "TP53", "LRP1B", "TBX3")

per_gene_tests <- somatic_tso500_PD_trace_variants_paired |>
    select(patnr, arm, genename) |>
    distinct() |>
    select(!patnr) |>
    group_by(arm, genename) |>
    summarize(n = n(), .groups = "drop") |>
    tidyr::pivot_wider(names_from = arm, values_from = n, values_fill = 0) |>
    group_by(genename) |>
    group_modify(~fisher_gene(.x)) |>
    ungroup() |>
    arrange(p_value) |>
    mutate(p_adj = p.adjust(p_value, method = "BH"))

top_genes_alteration <- per_gene_tests |>
    mutate(total_mut = n_mut_A + n_mut_B) |>
    filter(genename %in% top_genes) |>
    arrange(desc(total_mut)) |>
    tibble::rowid_to_column(var = "gene_order") |>
    mutate(p_adj = p.adjust(p_value, method = "fdr"))

patients_with_mut_in_top <- somatic_tso500_PD_trace_variants_paired |>
    dplyr::right_join(top_genes_alteration, by = join_by(genename)) |>
    pull(patnr) |>
    unique()

missing_patients <- meta_data |>
    filter(!patnr %in% patients_with_mut_in_top, timepoint == "PD")

mutation_labels <- somatic_tso500_PD_trace_variants_paired |>
    right_join(top_genes_alteration, by = join_by(genename)) |>
    bind_rows(missing_patients) |>
    group_by(patnr, arm, genename, diff_variant) |>
    summarize(n = n(), .groups = "drop") |>
    tidyr::pivot_wider(names_from = diff_variant, values_from = n, values_fill = 0) |>
    mutate(
        mutation_label = case_when(
            both == 1 & PD == 0 ~ "maintained",
            both > 1 & PD == 0 ~ "maintained_polyclonal",
            both > 0 & PD > 0  ~ "acquired_maintained",
            both == 0 & PD == 1 ~ "acquired",
            both == 0 & PD > 1 ~ "acquired_polyclonal",
            TRUE ~ NA
        ),
        mutation_label_id = case_when(
            mutation_label == "maintained"            ~ 2,
            mutation_label == "maintained_polyclonal" ~ 1,
            mutation_label == "acquired_maintained"   ~ 3,
            mutation_label == "acquired"              ~ 4,
            mutation_label == "acquired_polyclonal"   ~ 5,
            TRUE ~ NA
        )
    )

esr1_burden_per_patient <- somatic_tso500_PD_trace_variants_paired |>
    filter(genename == "ESR1") |>
    group_by(patnr, arm) |>
    summarize(esr1_burden = n(), .groups = "drop")

patient_order_somatic_variants <- mutation_labels |>
    select(patnr, arm, genename, mutation_label_id) |>
    tidyr::pivot_wider(names_from = genename, values_from = mutation_label_id, values_fill = 0) |>
    select(c(patnr, arm, all_of(top_genes_alteration$genename)), everything()) |>
    mutate(n_total = rowSums(across(all_of(top_genes_alteration$genename)) != 0)) |>
    left_join(esr1_burden_per_patient, by = join_by(patnr, arm)) |>
    mutate(esr1_burden = tidyr::replace_na(esr1_burden, 0)) |>
    arrange(esr1_burden, n_total, across(all_of(top_genes_alteration$genename), ~ .x)) |>
    mutate(patient_ordering = row_number()) |>
    relocate(patnr, patient_ordering)

gene_ordering_focal_cna <- copy_number_genes |>
    filter(timepoint == "PD", arm == "arm_A") |>
    group_by(genename) |>
    summarize(n = n(), .groups = "drop") |>
    mutate(PC_perc_focal_cna = n / 24) |>
    filter(PC_perc_focal_cna > 0.2) |>
    arrange(PC_perc_focal_cna) |>
    mutate(gene_ordering = row_number()) |>
    select(genename, gene_ordering)
