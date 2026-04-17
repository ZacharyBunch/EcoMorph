  #### Libraries ####
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tidyr)
  
  #### Read in ####
  calibration <- read_excel("EcoMorph_Callibration_Data.xlsx")
  batch_results_new <- read.csv("batch_results_4_10_2026.csv", check.names = FALSE)
  
  #### Standardize Sample IDs ####
  batch_results_new <- batch_results_new %>%
    mutate(
      grid_pos = specimen_id,
      base = str_remove(image, "\\.png$"),
      base = str_remove(base, "_\\d{8}_\\d{6}$"),
      base = str_remove(base, "_Bombus_Impatiens$"),
      base = str_remove(base, "_Xylocopa_virginica$"),
      base = str_remove(base, "_Dolichovespula_Maculata$"),
      Sample = paste0(base, "_", grid_pos)
    ) %>%
    select(-base, -grid_pos)
  
  #### Fix Spelling ####
  calibration <- calibration %>%
    mutate(Sample = str_replace(Sample, "Wentsvile", "Wentsville"))
  
  #### Verify Matches ####
  matched_new <- sum(batch_results_new$Sample %in% calibration$Sample)
  cat("Batch samples matched in calibration:", matched_new, "/", nrow(batch_results_new), "\n")
  
  unmatched <- batch_results_new$Sample[!batch_results_new$Sample %in% calibration$Sample]
  if (length(unmatched) > 0) {
    cat("Unmatched batch samples:\n")
    print(head(unique(unmatched), 20))
  }
  
  #### Identify Metric Columns ####
  exclude_patterns <- c("rgb", "polygon", "pixels_per", "total_specimens",
                        "rotation_angle", "coverage_percent", "ITD")
  
  metric_cols <- batch_results_new %>%
    select(where(is.numeric)) %>%
    names()
  metric_cols <- metric_cols[!grepl(paste(exclude_patterns, collapse = "|"), metric_cols)]
  cat("Metrics to correlate:\n")
  print(metric_cols)
  
  #### Join Calibration and Batch ####
  itd_col <- ifelse("ITD_mm" %in% names(calibration), "ITD_mm", "ITD (mm)")
  cat("Using calibration ITD column:", itd_col, "\n")
  
  itd_compare <- inner_join(
    calibration %>% select(Sample, ITD_hand = all_of(itd_col), Species),
    batch_results_new %>% select(Sample, all_of(metric_cols)),
    by = "Sample"
  ) %>%
    mutate(box = str_remove(Sample, "_[A-Z]\\d+$"))
  
  cat("Specimens available for correlation:", nrow(itd_compare), "\n")
  
  #### Global Correlation Summary ####
  cor_results <- data.frame(
    metric = metric_cols,
    n = NA_integer_,
    pearson_r = NA_real_,
    pearson_p = NA_real_,
    spearman_rho = NA_real_,
    spearman_p = NA_real_
  )
  
  for (i in seq_along(metric_cols)) {
    complete <- itd_compare %>% filter(!is.na(ITD_hand), !is.na(.data[[metric_cols[i]]]))
    if (nrow(complete) < 3) next
    rp <- cor.test(complete$ITD_hand, complete[[metric_cols[i]]], method = "pearson")
    rs <- cor.test(complete$ITD_hand, complete[[metric_cols[i]]], method = "spearman")
    cor_results$n[i] <- nrow(complete)
    cor_results$pearson_r[i] <- round(rp$estimate, 4)
    cor_results$pearson_p[i] <- rp$p.value
    cor_results$spearman_rho[i] <- round(rs$estimate, 4)
    cor_results$spearman_p[i] <- rs$p.value
  }
  
  cor_results <- cor_results %>%
    filter(!is.na(pearson_r)) %>%
    arrange(desc(abs(pearson_r)))
  
  cat("\n===== Correlation Summary (sorted by |r|) =====\n")
  print(cor_results, row.names = FALSE)
  
  #### Per-Box Correlation Plots ####
  annot_box <- itd_compare %>%
    pivot_longer(cols = all_of(cor_results$metric), names_to = "metric", values_to = "value") %>%
    filter(!is.na(ITD_hand), !is.na(value)) %>%
    group_by(box, metric) %>%
    filter(n() >= 3) %>%
    summarise(
      r = round(cor(ITD_hand, value, use = "complete.obs"), 3),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(metric = factor(metric, levels = cor_results$metric))
  
  plot_data_box <- itd_compare %>%
    pivot_longer(cols = all_of(cor_results$metric), names_to = "metric", values_to = "value") %>%
    filter(!is.na(ITD_hand), !is.na(value)) %>%
    mutate(metric = factor(metric, levels = cor_results$metric))
  
  for (b in unique(plot_data_box$box)) {
    box_data <- plot_data_box %>% filter(box == b)
    box_annot <- annot_box %>% filter(box == b)
    
    if (nrow(box_data) == 0 || nrow(box_annot) == 0) next
    
    sp_label <- itd_compare %>%
      filter(box == b) %>%
      pull(Species) %>%
      unique() %>%
      na.omit() %>%
      paste(collapse = ", ")
    
    p <- ggplot(box_data, aes(x = ITD_hand, y = value)) +
      geom_point(alpha = 0.4, size = 1) +
      geom_smooth(method = "lm", se = TRUE, color = "steelblue", linewidth = 0.7) +
      geom_text(data = box_annot,
                aes(label = paste0("r = ", r, "\nn = ", n)),
                x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 2.5,
                inherit.aes = FALSE) +
      facet_wrap(~ metric, scales = "free_y", ncol = 4) +
      labs(x = "Hand-measured ITD (mm)",
           y = "EcoMorph Measurement (px)",
           title = paste0("Box: ", b, "  |  Species: ", sp_label)) +
      theme_classic() +
      theme(strip.text = element_text(face = "bold", size = 7))
    print(p)
  }
  
  #### Summary Table: Per-Box Per-Metric Pearson r and n ####
  summary_r <- annot_box %>%
    select(box, metric, r) %>%
    pivot_wider(names_from = metric, values_from = r, names_prefix = "r_")
  
  summary_n <- annot_box %>%
    select(box, metric, n) %>%
    pivot_wider(names_from = metric, values_from = n, names_prefix = "n_")
  
  summary_table <- left_join(summary_r, summary_n, by = "box") %>%
    arrange(box)
  
  print(summary_table, n = Inf)
  
  # Write to CSV
  write.csv(summary_table, "EcoMorph_correlation_summary_4_10_2026.csv", row.names = FALSE)
  
  
  #### Error Magnitude Table: RMSE and MAE per Box per Metric ####
  #### Error Magnitude Table: RMSE and MAE per Box per Metric ####
  error_table <- itd_compare %>%
    pivot_longer(cols = all_of(cor_results$metric), names_to = "metric", values_to = "value") %>%
    filter(!is.na(ITD_hand), !is.na(value)) %>%
    group_by(box, metric) %>%
    filter(n() >= 3) %>%
    summarise(
      n = n(),
      slope     = round(coef(lm(ITD_hand ~ value))[2], 6),
      intercept = round(coef(lm(ITD_hand ~ value))[1], 4),
      rmse      = round(sqrt(mean(resid(lm(ITD_hand ~ value))^2)), 4),
      mae       = round(mean(abs(resid(lm(ITD_hand ~ value)))), 4),
      .groups = "drop"
    ) %>%
    arrange(box, rmse)
  
  cat("\n===== Error Magnitude Summary =====\n")
  print(error_table, n = Inf)
  
  write.csv(error_table, "EcoMorph_error_magnitude_4_10_2026.csv", row.names = FALSE)
  
  
  #### ITD All ####
  # Compute per-metric annotations pooled across all boxes
  annot_pooled <- itd_compare %>%
    pivot_longer(cols = all_of(cor_results$metric), names_to = "metric", values_to = "value") %>%
    filter(!is.na(ITD_hand), !is.na(value)) %>%
    group_by(metric) %>%
    filter(n() >= 3) %>%
    summarise(
      r = round(cor(ITD_hand, value, use = "complete.obs"), 3),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(metric = factor(metric, levels = cor_results$metric))
  
  plot_data_pooled <- itd_compare %>%
    pivot_longer(cols = all_of(cor_results$metric), names_to = "metric", values_to = "value") %>%
    filter(!is.na(ITD_hand), !is.na(value)) %>%
    mutate(metric = factor(metric, levels = cor_results$metric))
  
  ggplot(plot_data_pooled, aes(x = ITD_hand, y = value, color = Species)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
    geom_text(data = annot_pooled,
              aes(label = paste0("r = ", r, "\nn = ", n)),
              x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 2.5,
              inherit.aes = FALSE) +
    facet_wrap(~ metric, scales = "free_y", ncol = 4) +
    labs(x = "Hand-measured ITD (mm)",
         y = "EcoMorph Measurement",
         title = "All Boxes Pooled: Hand-measured ITD vs EcoMorph Metrics",
         color = "Species") +
    theme_classic() +
    theme(strip.text = element_text(face = "bold", size = 7),
          legend.position = "bottom")
  
  
  #### Standalone Plot: Body Area (cm²) — Publication Ready ####
  
  library(ggplot2)
  
  body_area_col <- names(batch_results_new)[
    grepl("body_area", names(batch_results_new), ignore.case = TRUE) &
      grepl("cm", names(batch_results_new), ignore.case = TRUE)
  ]
  
  itd_body_area <- inner_join(
    calibration %>% select(Sample, ITD_hand = all_of(itd_col), Species),
    batch_results_new %>% select(Sample, all_of(body_area_col)),
    by = "Sample"
  ) %>%
    mutate(ITD_hand_cm = ITD_hand) %>%
    pivot_longer(cols = all_of(body_area_col), names_to = "metric", values_to = "value") %>%
    filter(!is.na(ITD_hand_cm), !is.na(value))
  
  annot_body_area <- itd_body_area %>%
    summarise(
      r = round(cor(ITD_hand_cm, value, use = "complete.obs"), 3),
      n = n()
    )
  
  species_palette <- c(
    "Bombus impatiens"        = "#E69F00",
    "Dolichovespula Maculata" = "#009E73",
    "Xylocopa virginica"      = "#0072B2"
  )
  
  ggplot(itd_body_area, aes(x = ITD_hand_cm, y = value, color = Species, fill = Species)) +
    geom_point(alpha = 0.55, size = 2, shape = 21, stroke = 0.3, color = "white") +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "black", fill = "grey75", linewidth = 0.9, alpha = 0.25) +
    annotate("text",
             x     = -Inf, y = Inf,
             label = paste0("r = ", annot_body_area$r, "\nn = ", annot_body_area$n),
             hjust = -0.12, vjust = 1.4,
             size  = 4.55, fontface = "italic", color = "grey25") +
    scale_color_manual(
      values = species_palette,
      labels = expression(italic("Bombus impatiens"),
                          italic("Dolichovespula maculata"),
                          italic("Xylocopa virginica"))
    ) +
    scale_fill_manual(
      values = species_palette,
      labels = expression(italic("Bombus impatiens"),
                          italic("Dolichovespula maculata"),
                          italic("Xylocopa virginica"))
    ) +
    labs(
      x     = "Intertegular Distance (mm)",
      y     = expression("Body Area (mm"^2*")"),
      color = NULL,
      fill  = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.title        = element_text(size = 16.8, color = "grey15"),
      axis.text.x       = element_text(size = 13.2, color = "grey25", margin = margin(t = 6)),
      axis.text.y       = element_text(size = 13.2, color = "grey25", margin = margin(r = 6)),
      axis.line         = element_line(color = "grey40", linewidth = 0.4),
      axis.ticks        = element_line(color = "grey40", linewidth = 0.4),
      legend.position   = "bottom",
      legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
      legend.text       = element_text(size = 10),
      legend.key.size   = unit(1.2, "cm"),
      plot.margin       = margin(12, 16, 8, 8)
    ) +
    guides(
      fill  = guide_legend(override.aes = list(size = 5, shape = 21, stroke = 0.3)),
      color = guide_legend(override.aes = list(size = 5, shape = 21, stroke = 0.3))
    )
  
  ggsave("body_area_ITD.png", width = 10, height = 8, dpi = 300)
  

  #### Standalone Plot: Thorax Width — Publication Ready ####
  
  library(ggplot2)
  
  thorax_width_col <- names(batch_results_new)[
    grepl("thorax_width", names(batch_results_new), ignore.case = TRUE) &
      grepl("cm", names(batch_results_new), ignore.case = TRUE)
  ]
  
  itd_thorax_width <- inner_join(
    calibration %>% select(Sample, ITD_hand = all_of(itd_col), Species),
    batch_results_new %>% select(Sample, all_of(thorax_width_col)),
    by = "Sample"
  ) %>%
    mutate(ITD_hand_cm = ITD_hand ) %>%
    pivot_longer(cols = all_of(thorax_width_col), names_to = "metric", values_to = "value") %>%
    filter(!is.na(ITD_hand_cm), !is.na(value))
  
  annot_thorax_width <- itd_thorax_width %>%
    summarise(
      r = round(cor(ITD_hand_cm, value, use = "complete.obs"), 3),
      n = n()
    )
  
  species_palette <- c(
    "Bombus impatiens"        = "#E69F00",
    "Dolichovespula Maculata" = "#009E73",
    "Xylocopa virginica"      = "#0072B2"
  )
  
  ggplot(itd_thorax_width, aes(x = ITD_hand_cm, y = value, color = Species, fill = Species)) +
    geom_point(alpha = 0.55, size = 2, shape = 21, stroke = 0.3, color = "white") +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "black", fill = "grey75", linewidth = 0.9, alpha = 0.25) +
    annotate("text",
             x     = -Inf, y = Inf,
             label = paste0("r = ", annot_thorax_width$r, "\nn = ", annot_thorax_width$n),
             hjust = -0.12, vjust = 1.4,
             size  = 4.55, fontface = "italic", color = "grey25") +
    scale_color_manual(
      values = species_palette,
      labels = expression(italic("Bombus impatiens"),
                          italic("Dolichovespula maculata"),
                          italic("Xylocopa virginica"))
    ) +
    scale_fill_manual(
      values = species_palette,
      labels = expression(italic("Bombus impatiens"),
                          italic("Dolichovespula maculata"),
                          italic("Xylocopa virginica"))
    ) +
    labs(
      x     = "Intertegular Distance (mm)",
      y     = "Thorax Width (mm)",
      color = NULL,
      fill  = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.title        = element_text(size = 16.8, color = "grey15"),
      axis.text.x       = element_text(size = 13.2, color = "grey25", margin = margin(t = 6)),
      axis.text.y       = element_text(size = 13.2, color = "grey25", margin = margin(r = 6)),
      axis.line         = element_line(color = "grey40", linewidth = 0.4),
      axis.ticks        = element_line(color = "grey40", linewidth = 0.4),
      legend.position   = "bottom",
      legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
      legend.text       = element_text(size = 10),
      legend.key.size   = unit(1.2, "cm"),
      plot.margin       = margin(12, 16, 8, 8)
    ) +
    guides(
      fill  = guide_legend(override.aes = list(size = 5, shape = 21, stroke = 0.3)),
      color = guide_legend(override.aes = list(size = 5, shape = 21, stroke = 0.3))
    )
  
  ggsave("thorax_width_ITD.png", width = 10, height = 8, dpi = 300)
