# Comparison of Ethnicity and Admixture

# Input data
genome_cna_1q_1p32_17p13_clinical_df |>
  select(public_id, reclassify_ethnicity) |>
  left_join(admixture_df,
            by = c("public_id" = "ID"))

# Proportions

genome_cna_1q_1p32_17p13_clinical_df |>
  select(public_id, reclassify_ethnicity) |>
  left_join(admixture_df,
            by = c("public_id" = "ID")) |>
  arrange(`1KG-EUR-like`) |>
  pivot_longer(cols = starts_with("1KG"),
               names_to = "1KG",
               values_to = "1KG_proportion") |>
  ggplot(aes(x = `1KG_proportion`,
             y = public_id,
             fill = `1KG`)) +
  geom_col(width = 1) +
  facet_grid(rows = vars(reclassify_ethnicity),
             scales = "free",
             space = "free") +
  labs(x = "1KG proportion",
       y = "Public ID",
       fill = "1KG population",
       title = "CoMMpass admixture estimates by ethnicity") +
  scale_fill_viridis_d() +
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 20) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0),
        strip.clip = "off",
        strip.background = element_blank(),
        legend.position = "right",
        legend.direction = "vertical")

# AFR High

genome_cna_1q_1p32_17p13_clinical_df |>
  select(public_id, reclassify_ethnicity) |>
  left_join(admixture_df,
            by = c("public_id" = "ID")) |>
  ggplot(aes(x = `1KG-AFR-like`,
             fill = reclassify_ethnicity)) +
  geom_histogram() +
  geom_vline(xintercept = 0.5,
             lty = 2) +
  annotate(geom = "label", x = 0.6, y = 50, label = "AFR High", hjust = 0) +
  annotate(geom = "label", x = 0.4, y = 50, label = "AFR Low", hjust = 1) +
  labs(y = "N patients",
       fill = "Ethnicity\n(self-reported)",
       title = "Genomic similarity to 1KG-AFR") +
  theme_bw(base_size = 20) +
  theme(legend.position = "bottom")

# PCA

pca_input_df <- genome_cna_1q_1p32_17p13_clinical_df |>
  select(public_id, reclassify_ethnicity, hispanic_latino) |>
  left_join(admixture_df,
            by = c("public_id" = "ID")) |>
  filter(!is.na(`1KG-AFR-like`))

set.seed(1)
pca_output_df <- pca_input_df |>
  select(-public_id, -reclassify_ethnicity, -hispanic_latino) |>
  t() |>
  prcomp()

pca_df <- pca_input_df |>
  bind_cols(pca_output_df$rotation) |>
  replace_na(list(reclassify_ethnicity = "NA", hispanic_latino = "NA"))

pca_df |>
  mutate(reclassify_ethnicity2 = reclassify_ethnicity) |>
  bind_rows(pca_df |>
              mutate(reclassify_ethnicity2 = "All together")) |>
  ggplot(aes(x = PC1,
             y = PC2)) +
  geom_point(aes(color = reclassify_ethnicity), shape = 16, alpha = 0.5) +
  #geom_point(aes(color = hispanic_latino), shape = 16, alpha = 0.5) +
  geom_segment(x = -0.004, y = -0.045, yend = Inf, lty = 3) +
  geom_segment(x = -Inf, y = -0.045, xend = Inf, lty = 3) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~ reclassify_ethnicity2,
             nrow = 2) +
  labs(title = "Admixture principal components by ethnicity",
       color = "Self-reported\nethnicity") +
  theme_bw(base_size = 20) +
  theme(strip.clip = "off",
        strip.background = element_blank()) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

# reclassify labels

reclassify_max_pop_df <- admixture_df |>
  pivot_longer(cols = starts_with("1KG")) |>
  group_by(ID) |>
  summarize(reclassify_max_pop = c("AFR", "EUR", "EAS", "PEL")[which.max(value)],
            .groups = "drop") |>
  rename("public_id" = "ID")

reclassify_pca_df <- pca_df |>
  mutate(reclassify_pca = case_when(PC1 > -0.004 & PC2 > -0.045 ~ "White",
                                  PC1 < -0.004 & PC2 > -0.045 ~ "Asian",
                                  PC2 < -0.045 ~ "Black",
                                  .default = reclassify_ethnicity)) |>
  select(public_id, reclassify_pca)

afr_high_df <- admixture_df |>
  mutate(afr_high = ifelse(`1KG-AFR-like` > 0.5, "AFR High", "AFR Low")) |>
  select(public_id = ID,
         afr_high)

pca_df |>
  left_join(reclassify_pca_df, by = "public_id") |>
  left_join(reclassify_max_pop_df, by = "public_id") |>
  left_join(afr_high_df, by = "public_id") |>
  count(reclassify_ethnicity, reclassify_pca, reclassify_max_pop, afr_high) |>
  arrange(desc(n)) |>
  mutate(prop = n/nrow(pca_df))

ancestry_reclassification_df <- pca_df |>
  left_join(reclassify_pca_df, by = "public_id") |>
  left_join(reclassify_max_pop_df, by = "public_id") |>
  left_join(afr_high_df, by = "public_id") |>
  select(public_id, reclassify_ethnicity, reclassify_pca, reclassify_max_pop, afr_high)

write_tsv(x = ancestry_reclassification_df,
          file = here::here(processed_data_dir, "ancestry_df.tsv"))
