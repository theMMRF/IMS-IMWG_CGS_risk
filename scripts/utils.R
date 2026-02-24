process_segment_data <- function(seg_df,
                                 region = NULL) {

  if (!is.null(region)) {

    seg_df <- seg_df |>
      overlap_region(region = region) |>
      cap_region(region = region)

  }

  seg_df |>
  mutate(public_id = str_sub(SAMPLE, 1, 9),
         timepoint = str_sub(SAMPLE, 11, 11),
         tissue_source = str_sub(SAMPLE, 13, 14),
         Chromosome = factor(Chromosome,
                             levels = c(str_c("chr", seq(1:22)), "chrX", "chrY"),
                             ordered = TRUE),
         Start_Mb = Start/1e6,
         End_Mb = End/1e6,
         length = End - Start,
         length_Mb = length/1e6,
         rescaled = 2*(2^Segment_Mean),
         winsorized = case_when(rescaled >= 4 ~ 4,
                                TRUE ~ rescaled),
         rounded = round(winsorized, 0),
         rounded_cat = case_when(rounded == 0 ~ "0 copies",
                                 rounded == 1 ~ "1 copy",
                                 rounded == 2 ~ "2 copies",
                                 rounded == 3 ~ "3 copies",
                                 rounded >= 4 ~ "4+ copies")
  )


}

overlap_region <- function(seg_df, region) {

  seg_df |>
    filter(Chromosome == region$chr &
             !(End < region$start | Start > region$end))

}

cap_region <- function(seg_df, region) {

  seg_df |>
    mutate(Start = ifelse(Start < region$start, region$start, Start),
           End = ifelse(End > region$end, region$end, End))

}

weighted_segment_mean <- function(seg_df) {

  seg_df |>
    group_by(SAMPLE, public_id, timepoint, tissue_source) |>
    summarize(n_segments = n(),
              total_length = sum(length),
              mean_rescaled = sum(length*rescaled)/total_length,
              .groups = "drop") |>
    mutate(mean_log2 = log2(mean_rescaled/2),
           mean_winsorized = ifelse(mean_rescaled > 4, 4, mean_rescaled),
           mean_rounded = round(mean_winsorized, 0),
           mean_rounded_cat = case_when(mean_rounded == 0 ~ "0 copies",
                                        mean_rounded == 1 ~ "1 copy",
                                        mean_rounded == 2 ~ "2 copies",
                                        mean_rounded == 3 ~ "3 copies",
                                        mean_rounded >= 4 ~ "4+ copies"))

}
