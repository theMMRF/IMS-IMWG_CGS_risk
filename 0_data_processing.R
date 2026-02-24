# Step 0: Data input and processing
# Steven Foltz

# This script reads in data,
# organizes gene expression and clinical data frames,
# processes copy number segment data,
# and saves output files for downstream use

# Libraries ----
library(tidyverse)

# Directories ----
proj_dir <- here::here()
data_dir <- here::here(proj_dir, "data")
tgen_dir <- here::here(data_dir, "CGS_from_TGen")
processed_data_dir <- here::here(proj_dir, "processed_data")
outputs_dir <- here::here(proj_dir, "outputs")
scripts_dir <- here::here(proj_dir, "scripts")

IA24_TGen_clinical_dir <- here::here(data_dir, "TGen_IA24_clinical")

source(here::here(scripts_dir, "utils.R"))

# Filenames ----
clinical_filename <- here::here(data_dir, "CoMMpass_clinical_n1143.updated_survival.txt")
clinical2_filename <- here::here(data_dir, "CoMMpass_NDMM_baseline_clinical.txt")
genome_cna_filename <- here::here(data_dir, "MMRF_CoMMpass_IA22_genome_gatk_cna.seg")
exome_cna_filename <- here::here(data_dir, "MMRF_CoMMpass_IA22_exome_gatk_cna.seg")
admixture_filename <- here::here(data_dir, "MMRF_exome_summary[41].xlsx")
tgen_df_filename <- here::here(tgen_dir, "MMRF_CoMMpass_IA24_Preliminary_imsIMWG_HRMM,tsv")

tgen_ia24_survival_filename <- here::here(IA24_TGen_clinical_dir, "MMRF_CoMMpass_IA24_STAND_ALONE_SURVIVAL.tsv")

# source of survival: IA24 harmonized
os_pfs_filename <- here::here(data_dir, "survival_pfs_deid.csv")
subject_filename <- here::here(data_dir, "subject_deid.csv")

# first line therapy: based on curated time to second line therapy
first_line_therapy_filename <- here::here(data_dir, "MMRF_CoMMpass_first_line_therapy.20260128.tsv")

# Define 1q and 1q21 region ----
chr1q <- list("chr" = "chr1", "start" = 135000000, "end" = 249000000)
chr1q21 <- list("chr" = "chr1", "start" = 143200001, "end" = 155100000)
chr1p32 <- list("chr" = "chr1", "start" = 50200001, "end" = 60800000)
chr17p13 <- list("chr" = "chr17", "start" = 0, "end" = 10800000)

# Read in data

clinical_df <- read_tsv(clinical_filename)
clinical2_df <- read_tsv(clinical2_filename) |>
  separate(Tumor_Sample_Barcode,
           into = c("public", "id", "visit", "source", "sort"),
           sep = "_") |>
  mutate("public_id" = str_c(public, id, sep = "_"))
genome_cna_df <- read_tsv(genome_cna_filename)
exome_cna_df <- read_tsv(exome_cna_filename)
admixture_df <- readxl::read_excel(admixture_filename)
tgen_df <- read_tsv(tgen_df_filename) |>
  filter(Reason_For_Collection == "Baseline",
         Source == "BM",
         IMS_IMWG_Genomic_Evaluable == 1)

os_pfs_df <- read_csv(os_pfs_filename) |>
  mutate(public_id = str_to_upper(public_id))

first_line_therapy_df <- read_tsv(first_line_therapy_filename)

tgen_ia24_survival_df <- read_tsv(tgen_ia24_survival_filename) |>
  select(public_id = PUBLIC_ID,
         tgen_ia24_censos = censos,
         tgen_ia24_ttcos = ttcos,
         tgen_ia24_censpfs = censpfs,
         tgen_ia24_ttcpfs = ttcpfs)

subject_df <- read_csv(subject_filename) |>
  mutate(public_id = str_to_upper(public_id)) |>
  select(public_id, subcohort)

# Process CNA data ----

genome_cna_processed_df <- genome_cna_df |>
  process_segment_data()

exome_cna_processed_df <- exome_cna_df |>
  process_segment_data()

genome_cna_chr1q_df <- genome_cna_processed_df |>
  process_segment_data(region = chr1q)

exome_cna_chr1q_df <- exome_cna_processed_df |>
  process_segment_data(region = chr1q)

genome_cna_chr1q21_df <- genome_cna_processed_df |>
  process_segment_data(region = chr1q21)

exome_cna_chr1q21_df <- exome_cna_processed_df |>
  process_segment_data(region = chr1q21)

genome_cna_chr1p32_df <- genome_cna_processed_df |>
  process_segment_data(region = chr1p32)

exome_cna_chr1p32_df <- exome_cna_processed_df |>
  process_segment_data(region = chr1p32)

genome_cna_chr17p13_df <- genome_cna_processed_df |>
  process_segment_data(region = chr17p13)

exome_cna_chr17p13_df <- exome_cna_processed_df |>
  process_segment_data(region = chr17p13)

# Weighted segment mean ----

genome_cna_chr1q_seg_mean_df <- genome_cna_chr1q_df |>
  weighted_segment_mean()

exome_cna_chr1q_seg_mean_df <- exome_cna_chr1q_df |>
  weighted_segment_mean()

genome_cna_chr1q21_seg_mean_df <- genome_cna_chr1q21_df |>
  weighted_segment_mean()

exome_cna_chr1q21_seg_mean_df <- exome_cna_chr1q21_df |>
  weighted_segment_mean()

genome_cna_chr1p32_seg_mean_df <- genome_cna_chr1p32_df |>
  weighted_segment_mean()

exome_cna_chr1p32_seg_mean_df <- exome_cna_chr1p32_df |>
  weighted_segment_mean()

genome_cna_chr17p13_seg_mean_df <- genome_cna_chr17p13_df |>
  weighted_segment_mean()

exome_cna_chr17p13_seg_mean_df <- exome_cna_chr17p13_df |>
  weighted_segment_mean()

# Combine data ----

# Combine genome 1q and 1p32 dfs with clinical data
# Limit to primary BM

genome_cna_1q_1p32_17p13_clinical_df <- genome_cna_chr1q_seg_mean_df |>
  filter(timepoint == 1,
         tissue_source == "BM") |>
  select(public_id,
         mean_log2_1q = mean_log2,
         mean_rescaled_1q = mean_rescaled,
         mean_winsorized_1q = mean_winsorized,
         mean_rounded_1q = mean_rounded,
         mean_rounded_cat_1q = mean_rounded_cat) |>
  left_join(genome_cna_chr1p32_seg_mean_df |>
              filter(timepoint == 1,
                     tissue_source == "BM") |>
              select(public_id,
                     mean_log2_1p32 = mean_log2,
                     mean_rescaled_1p32 = mean_rescaled,
                     mean_winsorized_1p32 = mean_winsorized,
                     mean_rounded_1p32 = mean_rounded,
                     mean_rounded_cat_1p32 = mean_rounded_cat),
            by = "public_id") |>
  left_join(genome_cna_chr17p13_seg_mean_df |>
              filter(timepoint == 1,
                     tissue_source == "BM") |>
              select(public_id,
                     mean_log2_17p13 = mean_log2,
                     mean_rescaled_17p13 = mean_rescaled,
                     mean_winsorized_17p13 = mean_winsorized,
                     mean_rounded_17p13 = mean_rounded,
                     mean_rounded_cat_17p13 = mean_rounded_cat),
            by = "public_id") |>
  left_join(clinical_df |>
              select(public_id, Gender, Ethnicity, Age, Serum_B2M, Creatinine, BMI,
                     ASCT_First, Trip_First,
                     censpfs, ttcpfs, censos, ttcos,
                     ISS_Stage, IMWG_Risk_Class, Cytogenetic_High_Risk,
                     TP53_Funct_Copies, TP53_NS_Mut_Count) |>
              mutate(reclassify_ethnicity = ifelse(Ethnicity == "Caucasian", "White", Ethnicity)),
            by = "public_id") |>
  left_join(clinical2_df |>
              select(public_id, "hispanic_latino" = Ethnicity, `1q21_amp`, `1q21_gain`,
                     `17p13_del`, `TP53 inactivation`,
                     `MAF/MAFB`, `t(11;14)`, `t(4;14)`),
            by = "public_id") |>
  left_join(tgen_df |>
              select(public_id = Patient_ID,
                     TP53_Mutation_Count,
                     CDKN2C_Mutation_Count,
                     RB1_Mutation_Count,
                     NSD2_t_4_14_pos,
                     MAF_t_14_16_pos,
                     MAFB_t_14_20_pos,
                     MAFA_t_8_14_pos),
            by = "public_id") |>
  left_join(os_pfs_df |>
              select(public_id,
                     final_survival_time_os = surivival_time_os,
                     final_censor_os = censor_os,
                     final_survival_time_pfs = survival_time_pfs,
                     final_censor_pfs = censor_pfs),
            by = "public_id") |>
  left_join(first_line_therapy_df |>
              select(public_id, received_asct, received_triplet),
            by = "public_id") |>
  left_join(tgen_ia24_survival_df,
            by = "public_id") |>
  left_join(subject_df,
            by = "public_id")

# Save processed CNV data ----

write_tsv(x = genome_cna_1q_1p32_17p13_clinical_df,
          file = here::here(processed_data_dir,
                            "genome_cna_1q_1p32_17p13_clinical_df.tsv"))

# Run ancestry
source(here::here(proj_dir, "scripts", "ethnicity_admixture.R"))
