# =============================================================================
# 01_wild_vs_lab.R - Figure 2: Wild vs Lab midgut bacterial communities
#
# Data: 570 biological samples (380 Wild, 190 Lab) x 328 MC + GG2 >=97% OTUs
# No rarefaction (see Methods 2.11). Wild host species not used here.
#
# Outputs
#   output/figures/Fig2_wild_vs_lab.pdf/.png
#   output/figures/FigS_wild_vs_lab_per_sample_composition.pdf/.png
#   output/figures/FigS_wild_vs_lab_top8_bars.pdf/.png   (old panel A)
#   output/stats/01_all_genera_by_origin.tsv              (heatmap values)
#   output/stats/01_*.tsv                 all test results
#   output/stats/01_results_summary.txt   numbers ready for Results 3.5
#   output/tables/01_per_sample_alpha_depth.tsv
# =============================================================================

source(here::here("R", "00_setup.R"))
check_inputs()

# ---- 1. Load ----------------------------------------------------------------
meta    <- read_metadata()
otu_map <- read_otu_map()
bio     <- read_otu_matrix("biological")

sd <- tibble(Sample_ID = rownames(bio)) |>
  left_join(meta |> select(Sample_ID, Origin, Lab_Taxon), by = "Sample_ID") |>
  mutate(
    Origin    = droplevels(Origin),
    depth     = as.numeric(rowSums(bio)),
    log_depth = log10(depth),
    observed  = as.numeric(rowSums(bio > 0)),
    shannon   = as.numeric(diversity(bio, index = "shannon")),
    pielou    = if_else(observed > 1, shannon / log(observed), NA_real_)
  )
stopifnot(identical(sd$Sample_ID, rownames(bio)),
          !anyNA(sd$Origin),
          identical(levels(sd$Origin), c("Wild", "Lab")))

save_table(sd, "01_per_sample_alpha_depth")

# x-axis labels with sample sizes, e.g. "Wild\n(n = 380)"
n_by_origin   <- count(sd, Origin)
origin_labels <- setNames(paste0(n_by_origin$Origin, "\n(n = ", n_by_origin$n, ")"),
                          as.character(n_by_origin$Origin))

fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, digits = 2, format = "g")))

# ---- 2. Wilcoxon helper -----------------------------------------------------
compare_origin <- function(df, metrics) {
  map_dfr(metrics, function(m) {
    x <- df[[m]][df$Origin == "Wild"]
    y <- df[[m]][df$Origin == "Lab"]
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    t <- wilcox.test(x, y, exact = FALSE)
    tibble(
      metric      = m,
      n_wild      = length(x),
      median_wild = median(x),
      q1_wild     = as.numeric(quantile(x, 0.25)),
      q3_wild     = as.numeric(quantile(x, 0.75)),
      n_lab       = length(y),
      median_lab  = median(y),
      q1_lab      = as.numeric(quantile(y, 0.25)),
      q3_lab      = as.numeric(quantile(y, 0.75)),
      W           = as.numeric(t$statistic),
      p           = t$p.value
    )
  })
}

# ---- 3. Sequencing depth (for the first sentence of Results 3.5) ------------
depth_tab <- compare_origin(sd, "depth")
save_table(depth_tab, "01_depth_wild_vs_lab", paths$stats)

# ---- 4. Genus composition (Fig. 2A + supplementary per-sample plot) ---------
rel <- bio / rowSums(bio)

genus_long <- as_tibble(rel, rownames = "Sample_ID") |>
  pivot_longer(-Sample_ID, names_to = "OTU", values_to = "rel") |>
  left_join(select(otu_map, OTU, Genus), by = "OTU") |>
  mutate(Genus_plot = lump_genus(Genus)) |>
  group_by(Sample_ID, Genus_plot) |>
  summarise(rel = sum(rel), .groups = "drop") |>
  left_join(select(sd, Sample_ID, Origin), by = "Sample_ID")

comp_mean <- genus_long |>
  group_by(Origin, Genus_plot) |>
  summarise(mean_rel = mean(rel), .groups = "drop")
save_table(comp_mean |> mutate(mean_percent = round(100 * mean_rel, 2)),
           "01_genus_mean_relative_abundance", paths$stats)

# Reproducibility check against script 51 (Lab Elizabethkingia = 45.91%)
eliz_lab <- comp_mean |> filter(Origin == "Lab", Genus_plot == "Elizabethkingia") |> pull(mean_rel)
if (abs(100 * eliz_lab - 45.91) > 0.01) warning("Lab Elizabethkingia mean differs from script 51!")

# Legend labels: genus names italic, "Unclassified at genus" / "Other" plain
genus_label_expr <- as.expression(lapply(names(genus_colours), function(g) {
  lab <- pretty_genus(g)
  if (g %in% focal_genera) bquote(italic(.(lab))) else lab
}))
fill_genus <- scale_fill_manual(values = genus_colours, breaks = names(genus_colours),
                                labels = genus_label_expr, drop = FALSE, name = "Genus")

p_comp <- ggplot(comp_mean, aes(x = Origin, y = mean_rel, fill = Genus_plot)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  fill_genus +
  scale_x_discrete(labels = origin_labels) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 1.0001)) +
  labs(x = NULL, y = "Mean relative abundance") +
  theme(legend.key.size = unit(3.5, "mm"), legend.text = element_text(size = 8))

# Fig. 2A = every genus as a heatmap (all genera, per feedback). The top-8
# stacked bar is kept as a supplementary/slide figure.
heat_tab <- genus_mean_table(bio, sd$Origin)
save_table(heat_tab, "01_all_genera_by_origin", paths$stats)
heat_x <- setNames(paste0(n_by_origin$Origin, "\n(", n_by_origin$n, ")"), as.character(n_by_origin$Origin))
p_heat <- genus_heatmap(heat_tab, x_labels = heat_x,
                        breaks = c(1e-3, 0.1, 10), legend_width = 26)
save_fig(p_comp, "FigS_wild_vs_lab_top8_bars", width = 110, height = 110)

# Supplementary: every sample, ordered by Elizabethkingia then Wolbachia
sample_order <- genus_long |>
  select(Sample_ID, Origin, Genus_plot, rel) |>
  pivot_wider(names_from = Genus_plot, values_from = rel) |>
  arrange(Origin, desc(Elizabethkingia), desc(Wolbachia)) |>
  pull(Sample_ID)

p_per_sample <- genus_long |>
  mutate(Sample_ID = factor(Sample_ID, levels = sample_order)) |>
  ggplot(aes(x = Sample_ID, y = rel, fill = Genus_plot)) +
  geom_col(width = 1) +
  fill_genus +
  facet_grid(~ Origin, scales = "free_x", space = "free_x",
             labeller = as_labeller(origin_labels)) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "Individual mosquito midguts", y = "Relative abundance") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_fig(p_per_sample, "FigS_wild_vs_lab_per_sample_composition", width = 250, height = 110)

# ---- 5. Alpha diversity (Fig. 2B) -------------------------------------------
alpha_metrics <- c("observed", "shannon", "pielou")
alpha_names   <- c(observed = "Observed\nOTUs", shannon = "Shannon\ndiversity",
                   pielou = "Pielou's\nevenness")

alpha_tab <- compare_origin(sd, alpha_metrics) |>
  mutate(p_BH = p.adjust(p, method = "BH"))
save_table(alpha_tab, "01_alpha_wild_vs_lab", paths$stats)

alpha_long <- sd |>
  select(Sample_ID, Origin, all_of(alpha_metrics)) |>
  pivot_longer(all_of(alpha_metrics), names_to = "metric", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(metric = factor(alpha_names[metric], levels = alpha_names))

alpha_lab <- alpha_tab |>
  mutate(metric = factor(alpha_names[metric], levels = alpha_names),
         label  = fmt_p(p_BH))

p_alpha <- ggplot(alpha_long, aes(x = Origin, y = value)) +
  geom_jitter(aes(colour = Origin), width = 0.2, height = 0, size = 0.7, alpha = 0.35) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5, linewidth = 0.4) +
  geom_text(data = alpha_lab, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.3, size = 3) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_colour_manual(values = origin_colours, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = NULL) +
  theme(strip.text = element_text(size = 9, lineheight = 0.9),
        strip.clip = "off")   # never cut titles off at the panel edge

# ---- 6. Beta diversity: Jaccard PCoA, PERMANOVA, PERMDISP (Fig. 2C) ---------
sd_df <- as.data.frame(sd)
d_jac <- vegdist(bio, method = "jaccard", binary = TRUE)

perm_origin <- adonis2(d_jac ~ Origin, data = sd_df, permutations = 999, by = "terms")
perm_depth  <- adonis2(d_jac ~ log_depth + Origin, data = sd_df,
                       permutations = 999, by = "margin")
disp        <- betadisper(d_jac, sd_df$Origin)
disp_test   <- permutest(disp, permutations = 999)

tidy_adonis <- function(a, model) {
  as.data.frame(a) |> rownames_to_column("term") |> mutate(model = model, .before = 1)
}
save_table(bind_rows(tidy_adonis(perm_origin, "Origin only"),
                     tidy_adonis(perm_depth,  "log10(depth) + Origin, marginal")),
           "01_permanova_wild_vs_lab", paths$stats)

disp_out <- tibble(
  Origin = levels(sd$Origin),
  mean_distance_to_centroid = as.numeric(tapply(disp$distances, disp$group, mean)),
  F = disp_test$tab[1, "F"],
  p = disp_test$tab[1, "Pr(>F)"]
)
save_table(disp_out, "01_permdisp_wild_vs_lab", paths$stats)

pc      <- cmdscale(d_jac, k = 2, eig = TRUE)
var_exp <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
pcoa_df <- sd |>
  mutate(PCo1 = pc$points[Sample_ID, 1], PCo2 = pc$points[Sample_ID, 2])

R2_o <- perm_depth["Origin", "R2"]
F_o  <- perm_depth["Origin", "F"]
p_o  <- perm_depth["Origin", "Pr(>F)"]
p_d  <- disp_test$tab[1, "Pr(>F)"]

p_pcoa <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, colour = Origin)) +
  geom_point(size = 1.1, alpha = 0.55) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 0.8) +
  scale_colour_manual(values = origin_colours, labels = origin_labels, name = NULL) +
  labs(x = sprintf("PCoA 1 (%.1f%%)", var_exp[1]),
       y = sprintf("PCoA 2 (%.1f%%)", var_exp[2]),
       subtitle = sprintf("Jaccard | PERMANOVA (controlling for depth): R² = %.3f, F = %.1f, %s | PERMDISP: %s",
                          R2_o, F_o, fmt_p(p_o), fmt_p(p_d))) +
  theme(plot.subtitle = element_text(size = 8), legend.position = "bottom") +
  coord_fixed()   # equal axis scaling: distances on the plot match Jaccard distances

# ---- 7. Assemble Figure 2 ---------------------------------------------------
fig2 <- (free(p_heat) | (p_alpha / p_pcoa + plot_layout(heights = c(1, 1.25)))) +
  plot_layout(widths = c(1.2, 2.6)) +
  plot_annotation(tag_levels = "A")
save_fig(fig2, "Fig2_wild_vs_lab", width = 215, height = 250)

# ---- 8. Sensitivity: drop samples with < 1,000 classified reads -------------
keep    <- sd$depth >= 1000
sd_s    <- sd[keep, ]
sd_s_df <- as.data.frame(sd_s)
d_s     <- vegdist(bio[keep, ], method = "jaccard", binary = TRUE)

perm_s  <- adonis2(d_s ~ log_depth + Origin, data = sd_s_df, permutations = 999, by = "margin")
disp_s  <- permutest(betadisper(d_s, sd_s_df$Origin), permutations = 999)
alpha_s <- compare_origin(sd_s, alpha_metrics) |> mutate(p_BH = p.adjust(p, method = "BH"))

save_table(bind_rows(
  tidy_adonis(perm_s, "SENSITIVITY >=1000 reads: log10(depth) + Origin, marginal"),
  tibble(model = "SENSITIVITY >=1000 reads: PERMDISP", term = "Origin",
         F = disp_s$tab[1, "F"], `Pr(>F)` = disp_s$tab[1, "Pr(>F)"])
), "01_sensitivity_min1000_beta", paths$stats)
save_table(alpha_s, "01_sensitivity_min1000_alpha", paths$stats)

# ---- 9. Plain-text summary for writing Results ------------------------------
row_of <- function(tab, m) tab[tab$metric == m, ]
desc_line <- function(r, nm) sprintf(
  "%-18s Wild median %.2f (IQR %.2f-%.2f, n=%d) vs Lab median %.2f (IQR %.2f-%.2f, n=%d); W = %.0f, %s%s",
  nm, r$median_wild, r$q1_wild, r$q3_wild, r$n_wild,
  r$median_lab, r$q1_lab, r$q3_lab, r$n_lab, r$W, fmt_p(r$p),
  if ("p_BH" %in% names(r)) paste0(" (BH-adjusted ", fmt_p(r$p_BH), ")") else "")

summary_lines <- c(
  "01_wild_vs_lab.R - results summary",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "SEQUENCING DEPTH (classified reads per sample)",
  desc_line(depth_tab, "Depth"),
  "",
  "ALPHA DIVERSITY (Wilcoxon rank-sum, BH across 3 metrics)",
  desc_line(row_of(alpha_tab, "observed"), "Observed OTUs"),
  desc_line(row_of(alpha_tab, "shannon"),  "Shannon"),
  desc_line(row_of(alpha_tab, "pielou"),   "Pielou"),
  "",
  "BETA DIVERSITY (Jaccard, presence/absence, 999 permutations)",
  sprintf("PERMANOVA Origin only:              R2 = %.3f, F = %.2f, %s",
          perm_origin["Origin", "R2"], perm_origin["Origin", "F"], fmt_p(perm_origin["Origin", "Pr(>F)"])),
  sprintf("PERMANOVA Origin | log10(depth):    R2 = %.3f, F = %.2f, %s", R2_o, F_o, fmt_p(p_o)),
  sprintf("PERMANOVA log10(depth) | Origin:    R2 = %.3f, F = %.2f, %s",
          perm_depth["log_depth", "R2"], perm_depth["log_depth", "F"], fmt_p(perm_depth["log_depth", "Pr(>F)"])),
  sprintf("PERMDISP: F = %.2f, %s; mean distance to centroid Wild %.3f, Lab %.3f",
          disp_test$tab[1, "F"], fmt_p(p_d),
          disp_out$mean_distance_to_centroid[1], disp_out$mean_distance_to_centroid[2]),
  sprintf("PCoA axes: PCo1 %.1f%%, PCo2 %.1f%% (of positive eigenvalues)", var_exp[1], var_exp[2]),
  "",
  sprintf("SENSITIVITY (samples >= 1,000 reads: %d Wild, %d Lab)",
          sum(sd_s$Origin == "Wild"), sum(sd_s$Origin == "Lab")),
  sprintf("PERMANOVA Origin | log10(depth):    R2 = %.3f, F = %.2f, %s",
          perm_s["Origin", "R2"], perm_s["Origin", "F"], fmt_p(perm_s["Origin", "Pr(>F)"])),
  sprintf("PERMDISP: F = %.2f, %s", disp_s$tab[1, "F"], fmt_p(disp_s$tab[1, "Pr(>F)"])),
  desc_line(row_of(alpha_s, "observed"), "Observed OTUs"),
  desc_line(row_of(alpha_s, "shannon"),  "Shannon"),
  desc_line(row_of(alpha_s, "pielou"),   "Pielou"),
  "",
  "TOP GENERA, mean relative abundance (%) Wild / Lab",
  comp_mean |>
    mutate(pct = 100 * mean_rel) |>
    select(Origin, Genus_plot, pct) |>
    pivot_wider(names_from = Origin, values_from = pct) |>
    mutate(line = sprintf("  %-28s %6.2f / %6.2f", pretty_genus(Genus_plot), Wild, Lab)) |>
    pull(line)
)
dir.create(paths$stats, showWarnings = FALSE, recursive = TRUE)
writeLines(summary_lines, file.path(paths$stats, "01_results_summary.txt"))
cat(summary_lines, sep = "\n")

# ---- 10. Figure S3: (A) mean top-8 bars + (B) per-mosquito profiles ----------
figS3 <- (p_comp + theme(legend.position = "none") | p_per_sample) +
  plot_layout(widths = c(1, 3.4)) +
  plot_annotation(tag_levels = "A")
save_fig(figS3, "FigS3_genus_composition", width = 260, height = 110)
