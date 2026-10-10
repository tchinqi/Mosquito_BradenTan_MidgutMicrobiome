# =============================================================================
# 06_lab_species.R - Lab Ae. albopictus vs Lab Culex
#
# Data: 190 laboratory-reared samples (165 Ae. albopictus, 25 Culex) x 328 OTUs
# Same methods as 01_wild_vs_lab.R (no rarefaction; Jaccard presence/absence;
# PERMANOVA with log10 depth, marginal; PERMDISP; BH within each set of tests).
# Lab samples only; wild host species are analysed in 11_host_species.R.
#
# Outputs
#   output/figures/Fig_lab_species.pdf/.png
#   output/figures/FigS_lab_species_top8_bars.pdf/.png   (old panel A)
#   output/stats/06_all_genera_by_lab_species.tsv         (heatmap values)
#   output/stats/06_*.tsv
#   output/stats/06_results_summary.txt   numbers ready for Results
# =============================================================================

source(here::here("R", "00_setup.R"))
check_inputs()

# ---- 1. Load ----------------------------------------------------------------
meta    <- read_metadata()
otu_map <- read_otu_map()
lab     <- read_otu_matrix("lab")        # 190 Lab samples

taxon_levels  <- c("Ae. albopictus", "Culex")
taxon_colours <- c("Ae. albopictus" = "#4a3aa7", "Culex" = "#eb6834")  # as in FigS_wolbachia

sd <- tibble(Sample_ID = rownames(lab)) |>
  left_join(select(meta, Sample_ID, Origin, Lab_Taxon), by = "Sample_ID") |>
  mutate(
    Taxon     = factor(Lab_Taxon, levels = taxon_levels),
    depth     = as.numeric(rowSums(lab)),
    log_depth = log10(depth),
    observed  = as.numeric(rowSums(lab > 0)),
    shannon   = as.numeric(diversity(lab, index = "shannon")),
    pielou    = if_else(observed > 1, shannon / log(observed), NA_real_)
  )
stopifnot(identical(sd$Sample_ID, rownames(lab)),
          all(sd$Origin == "Lab"), !anyNA(sd$Taxon),
          sum(sd$Taxon == "Ae. albopictus") == 165, sum(sd$Taxon == "Culex") == 25)
save_table(sd, "06_per_sample_alpha_depth_lab_species")

n_tax <- table(sd$Taxon)
# Axis / legend labels with italic species names and sample sizes
tax_label_expr <- as.expression(list(
  bquote(atop(italic("Ae. albopictus"), "(n = " * .(as.integer(n_tax[1])) * ")")),
  bquote(atop(italic("Culex") ~ "spp.", "(n = " * .(as.integer(n_tax[2])) * ")"))))
tax_legend_expr <- as.expression(list(
  bquote(italic("Ae. albopictus") ~ "(n = " * .(as.integer(n_tax[1])) * ")"),
  bquote(italic("Culex") ~ "spp. (n = " * .(as.integer(n_tax[2])) * ")")))

fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, digits = 2, format = "g")))

compare_taxa <- function(df, metrics) {
  map_dfr(metrics, function(m) {
    x <- df[[m]][df$Taxon == "Ae. albopictus"]; y <- df[[m]][df$Taxon == "Culex"]
    x <- x[!is.na(x)]; y <- y[!is.na(y)]
    t <- wilcox.test(x, y, exact = FALSE)
    tibble(metric = m,
           n_alb = length(x), median_alb = median(x),
           q1_alb = as.numeric(quantile(x, 0.25)), q3_alb = as.numeric(quantile(x, 0.75)),
           n_clx = length(y), median_clx = median(y),
           q1_clx = as.numeric(quantile(y, 0.25)), q3_clx = as.numeric(quantile(y, 0.75)),
           W = as.numeric(t$statistic), p = t$p.value)
  })
}

# ---- 2. Sequencing depth ----------------------------------------------------
depth_tab <- compare_taxa(sd, "depth")
save_table(depth_tab, "06_depth_lab_species", paths$stats)

# ---- 3. Genus composition (Fig. A) ------------------------------------------
# Phytobacter is ~33% of Culex reads but is not one of the 8 house genera, so it
# gets its own colour HERE (brown, unused elsewhere); all other colours unchanged.
genus_colours_lab <- c(genus_colours[focal_genera], "Phytobacter" = "#8c5a2b",
                       genus_colours[c("Unclassified_at_genus", "Other")])
focal_lab <- c(focal_genera, "Phytobacter")
lump_genus_lab <- function(genus) {
  case_when(genus %in% focal_lab              ~ genus,
            genus == "Unclassified_at_genus" ~ "Unclassified_at_genus",
            TRUE                             ~ "Other") |>
    factor(levels = names(genus_colours_lab))
}

genus_long <- as_tibble(lab / rowSums(lab), rownames = "Sample_ID") |>
  pivot_longer(-Sample_ID, names_to = "OTU", values_to = "rel") |>
  left_join(select(otu_map, OTU, Genus), by = "OTU") |>
  mutate(Genus_plot = lump_genus_lab(Genus)) |>
  group_by(Sample_ID, Genus_plot) |>
  summarise(rel = sum(rel), .groups = "drop") |>
  left_join(select(sd, Sample_ID, Taxon), by = "Sample_ID")

comp_mean <- genus_long |>
  group_by(Taxon, Genus_plot) |>
  summarise(mean_rel = mean(rel), .groups = "drop")
save_table(comp_mean |> mutate(mean_percent = round(100 * mean_rel, 2)),
           "06_genus_mean_relative_abundance_lab_species", paths$stats)

genus_label_expr <- as.expression(lapply(names(genus_colours_lab), function(g) {
  lab_g <- pretty_genus(g)
  if (g %in% focal_lab) bquote(italic(.(lab_g))) else lab_g
}))

p_comp <- ggplot(comp_mean, aes(x = Taxon, y = mean_rel, fill = Genus_plot)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = genus_colours_lab, breaks = names(genus_colours_lab),
                    labels = genus_label_expr, drop = FALSE, name = "Genus") +
  scale_x_discrete(labels = tax_label_expr) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 1.0001)) +
  labs(x = NULL, y = "Mean relative abundance") +
  theme(legend.key.size = unit(3.5, "mm"), legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 8))

# Panel A = every genus detected in Lab mosquitoes, as a heatmap. The stacked
# bar (top 8 + Phytobacter) is kept as a supplementary/slide figure.
heat_tab <- genus_mean_table(lab, sd$Taxon)
save_table(heat_tab, "06_all_genera_by_lab_species", paths$stats)
heat_x_expr <- setNames(as.expression(list(
  bquote(atop(italic("Ae. alb."), "(" * .(as.integer(n_tax[1])) * ")")),
  bquote(atop(italic("Culex"), "(" * .(as.integer(n_tax[2])) * ")")))), taxon_levels)
p_heat <- genus_heatmap(heat_tab, x_labels = function(b) heat_x_expr[b],
                        breaks = c(1e-3, 0.1, 10), legend_width = 26)
save_fig(p_comp, "FigS_lab_species_top8_bars", width = 120, height = 110)

# ---- 4. Alpha diversity (Fig. B) --------------------------------------------
alpha_metrics <- c("observed", "shannon", "pielou")
alpha_names   <- c(observed = "Observed\nOTUs", shannon = "Shannon\ndiversity",
                   pielou = "Pielou's\nevenness")

alpha_tab <- compare_taxa(sd, alpha_metrics) |> mutate(p_BH = p.adjust(p, method = "BH"))
save_table(alpha_tab, "06_alpha_lab_species", paths$stats)

alpha_long <- sd |>
  select(Sample_ID, Taxon, all_of(alpha_metrics)) |>
  pivot_longer(all_of(alpha_metrics), names_to = "metric", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(metric = factor(alpha_names[metric], levels = alpha_names))
alpha_lab <- alpha_tab |>
  mutate(metric = factor(alpha_names[metric], levels = alpha_names), label = fmt_p(p_BH))

p_alpha <- ggplot(alpha_long, aes(x = Taxon, y = value)) +
  geom_jitter(aes(colour = Taxon), width = 0.2, height = 0, size = 0.7, alpha = 0.4) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5, linewidth = 0.4) +
  geom_text(data = alpha_lab, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.3, size = 2.7) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_colour_manual(values = taxon_colours, guide = "none") +
  scale_x_discrete(labels = c("Ae. albopictus" = "Ae. alb.", "Culex" = "Culex")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = NULL) +
  theme(strip.text = element_text(size = 9, lineheight = 0.9), strip.clip = "off",
        axis.text.x = element_text(face = "italic", size = 8, angle = 45, hjust = 1)) +
  coord_cartesian(clip = "off")

# ---- 5. Beta diversity (Fig. C) ---------------------------------------------
sd_df <- as.data.frame(sd)
d_jac <- vegdist(lab, method = "jaccard", binary = TRUE)

perm_tax   <- adonis2(d_jac ~ Taxon, data = sd_df, permutations = 999, by = "terms")
perm_depth <- adonis2(d_jac ~ log_depth + Taxon, data = sd_df, permutations = 999, by = "margin")
disp       <- betadisper(d_jac, sd_df$Taxon)
disp_test  <- permutest(disp, permutations = 999)

tidy_adonis <- function(a, model) {
  as.data.frame(a) |> rownames_to_column("term") |> mutate(model = model, .before = 1)
}
save_table(bind_rows(tidy_adonis(perm_tax, "Taxon only"),
                     tidy_adonis(perm_depth, "log10(depth) + Taxon, marginal")),
           "06_permanova_lab_species", paths$stats)

disp_out <- tibble(Taxon = levels(sd$Taxon),
                   mean_distance_to_centroid = as.numeric(tapply(disp$distances, disp$group, mean)),
                   F = disp_test$tab[1, "F"], p = disp_test$tab[1, "Pr(>F)"])
save_table(disp_out, "06_permdisp_lab_species", paths$stats)

pc      <- cmdscale(d_jac, k = 2, eig = TRUE)
var_exp <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
pcoa_df <- sd |> mutate(PCo1 = pc$points[Sample_ID, 1], PCo2 = pc$points[Sample_ID, 2])

R2_t <- perm_depth["Taxon", "R2"]; F_t <- perm_depth["Taxon", "F"]
p_t  <- perm_depth["Taxon", "Pr(>F)"]; p_d <- disp_test$tab[1, "Pr(>F)"]

p_pcoa <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, colour = Taxon)) +
  geom_point(size = 1.2, alpha = 0.6) +
  stat_ellipse(type = "t", level = 0.95, linewidth = 0.8) +
  scale_colour_manual(values = taxon_colours, labels = tax_legend_expr, name = NULL) +
  labs(x = sprintf("PCoA 1 (%.1f%%)", var_exp[1]), y = sprintf("PCoA 2 (%.1f%%)", var_exp[2]),
       subtitle = sprintf("Jaccard | PERMANOVA (controlling for depth): R² = %.3f, F = %.1f, %s | PERMDISP: %s",
                          R2_t, F_t, fmt_p(p_t), fmt_p(p_d))) +
  theme(plot.subtitle = element_text(size = 8), legend.position = "bottom") +
  coord_fixed()   # equal axis scaling: distances on the plot match Jaccard distances

fig <- (free(p_heat) | (p_alpha / p_pcoa + plot_layout(heights = c(1, 1.25)))) +
  plot_layout(widths = c(1.2, 2.6)) +
  plot_annotation(tag_levels = "A")
save_fig(fig, "Fig_lab_species", width = 215, height = 230)

# ---- 6. Sensitivity: samples with >= 1,000 classified reads -----------------
keep    <- sd$depth >= 1000
sd_s    <- sd[keep, ] |> mutate(Taxon = droplevels(Taxon))
d_s     <- vegdist(lab[keep, ], method = "jaccard", binary = TRUE)
perm_s  <- adonis2(d_s ~ log_depth + Taxon, data = as.data.frame(sd_s), permutations = 999, by = "margin")
disp_s  <- permutest(betadisper(d_s, sd_s$Taxon), permutations = 999)
alpha_s <- compare_taxa(sd_s, alpha_metrics) |> mutate(p_BH = p.adjust(p, method = "BH"))
save_table(bind_rows(
  tidy_adonis(perm_s, "SENSITIVITY >=1000 reads: log10(depth) + Taxon, marginal"),
  tibble(model = "SENSITIVITY >=1000 reads: PERMDISP", term = "Taxon",
         F = disp_s$tab[1, "F"], `Pr(>F)` = disp_s$tab[1, "Pr(>F)"])
), "06_sensitivity_min1000_beta", paths$stats)
save_table(alpha_s, "06_sensitivity_min1000_alpha", paths$stats)

# ---- 7. Genus tests (genera in >= 10% of the 190 Lab samples) ---------------
genus_of   <- otu_map$Genus[match(colnames(lab), otu_map$OTU)]
gen_counts <- t(rowsum(t(lab), group = genus_of))
gen_rel    <- gen_counts / rowSums(gen_counts)
is_a <- sd$Taxon == "Ae. albopictus"; is_c <- sd$Taxon == "Culex"
prev_all <- colMeans(gen_counts > 0)
tested   <- setdiff(names(prev_all)[prev_all >= 0.10], "Unclassified_at_genus")

gen_res <- map_dfr(tested, function(g) {
  x <- gen_rel[is_a, g]; y <- gen_rel[is_c, g]
  ft <- fisher.test(matrix(c(sum(gen_counts[is_a, g] > 0), sum(gen_counts[is_a, g] == 0),
                             sum(gen_counts[is_c, g] > 0), sum(gen_counts[is_c, g] == 0)), nrow = 2))
  wt <- wilcox.test(x, y, exact = FALSE)
  tibble(genus = g,
         mean_pct_alb = 100 * mean(x), mean_pct_clx = 100 * mean(y),
         prev_pct_alb = 100 * mean(gen_counts[is_a, g] > 0),
         prev_pct_clx = 100 * mean(gen_counts[is_c, g] > 0),
         W = as.numeric(wt$statistic), p_abundance = wt$p.value,
         fisher_OR = as.numeric(ft$estimate), p_prevalence = ft$p.value)
}) |>
  mutate(p_abundance_BH  = p.adjust(p_abundance,  method = "BH"),
         p_prevalence_BH = p.adjust(p_prevalence, method = "BH"),
         rank_higher = if_else(W > sum(is_a) * sum(is_c) / 2, "Ae. albopictus", "Culex"),
         higher_in   = if_else(p_abundance_BH >= 0.05, "n.s.", rank_higher)) |>
  arrange(p_abundance_BH)
save_table(gen_res, "06_genus_tests_lab_species", paths$stats)

# ---- 8. Plain-text summary --------------------------------------------------
row_of <- function(tab, m) tab[tab$metric == m, ]
desc_line <- function(r, nm) sprintf(
  "%-15s Ae. alb median %.2f (IQR %.2f-%.2f, n=%d) vs Culex median %.2f (IQR %.2f-%.2f, n=%d); W = %.1f, %s%s",
  nm, r$median_alb, r$q1_alb, r$q3_alb, r$n_alb, r$median_clx, r$q1_clx, r$q3_clx, r$n_clx,
  r$W, fmt_p(r$p), if ("p_BH" %in% names(r)) paste0(" (BH ", fmt_p(r$p_BH), ")") else "")

summary_lines <- c(
  "06_lab_species.R - results summary (Lab Ae. albopictus vs Lab Culex)",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "SEQUENCING DEPTH", desc_line(depth_tab, "Depth"),
  "",
  "ALPHA DIVERSITY (Wilcoxon, BH across 3 metrics)",
  desc_line(row_of(alpha_tab, "observed"), "Observed OTUs"),
  desc_line(row_of(alpha_tab, "shannon"),  "Shannon"),
  desc_line(row_of(alpha_tab, "pielou"),   "Pielou"),
  "",
  "BETA DIVERSITY (Jaccard presence/absence, 999 permutations)",
  sprintf("PERMANOVA Taxon only:             R2 = %.3f, F = %.2f, %s",
          perm_tax["Taxon", "R2"], perm_tax["Taxon", "F"], fmt_p(perm_tax["Taxon", "Pr(>F)"])),
  sprintf("PERMANOVA Taxon | log10(depth):   R2 = %.3f, F = %.2f, %s", R2_t, F_t, fmt_p(p_t)),
  sprintf("PERMANOVA log10(depth) | Taxon:   R2 = %.3f, F = %.2f, %s",
          perm_depth["log_depth", "R2"], perm_depth["log_depth", "F"], fmt_p(perm_depth["log_depth", "Pr(>F)"])),
  sprintf("PERMDISP: F = %.2f, %s; distance to centroid Ae. alb %.3f, Culex %.3f",
          disp_test$tab[1, "F"], fmt_p(p_d), disp_out$mean_distance_to_centroid[1],
          disp_out$mean_distance_to_centroid[2]),
  sprintf("PCoA axes: PCo1 %.1f%%, PCo2 %.1f%%", var_exp[1], var_exp[2]),
  "",
  sprintf("SENSITIVITY (>= 1,000 reads: %d Ae. alb, %d Culex)",
          sum(sd_s$Taxon == "Ae. albopictus"), sum(sd_s$Taxon == "Culex")),
  sprintf("PERMANOVA Taxon | log10(depth):   R2 = %.3f, F = %.2f, %s",
          perm_s["Taxon", "R2"], perm_s["Taxon", "F"], fmt_p(perm_s["Taxon", "Pr(>F)"])),
  sprintf("PERMDISP: F = %.2f, %s", disp_s$tab[1, "F"], fmt_p(disp_s$tab[1, "Pr(>F)"])),
  "",
  sprintf("GENUS TESTS (%d genera in >= 10%% of Lab samples)", length(tested)),
  sprintf("  %-26s Ae.alb %5.1f%% vs Culex %5.1f%% (BH %s) | prev %3.0f%% vs %3.0f%% (BH %s) | %s",
          pretty_genus(gen_res$genus), gen_res$mean_pct_alb, gen_res$mean_pct_clx,
          fmt_p(gen_res$p_abundance_BH), gen_res$prev_pct_alb, gen_res$prev_pct_clx,
          fmt_p(gen_res$p_prevalence_BH), gen_res$higher_in)
)
writeLines(summary_lines, file.path(paths$stats, "06_results_summary.txt"))
cat(summary_lines, sep = "\n")