# =============================================================================
# 02_locations.R - Report Figure 4: Wild mosquitoes across collection locations
#
# Statistics: the 345 Wild samples with confirmed locations ONLY
#   Mandai Zoo (MZ) 194, Kranji Marshes (KM) 88, Bird Paradise (BP) 39,
#   River Safari (RS) 24. The 35 unmapped Wild samples are excluded.
# Lab (n = 190) is drawn in grey on the composition plot and PCoA as a
# visual reference only - it never enters a test.
# No rarefaction. Wild host species are analysed in 11_host_species.R.
#
# Outputs
#   output/figures/Fig5_locations.pdf/.png
#   output/figures/FigS_locations_wild_only_pcoa.pdf/.png
#   output/figures/FigS_locations_top8_bars.pdf/.png   (old panel A)
#   output/stats/02_all_genera_by_location.tsv          (heatmap values)
#   output/stats/02_*.tsv
#   output/stats/02_results_summary.txt
# =============================================================================

source(here::here("R", "00_setup.R"))
check_inputs()

# ---- 1. Load ----------------------------------------------------------------
meta    <- read_metadata()
otu_map <- read_otu_map()
wild    <- read_otu_matrix("location")   # 345 mapped Wild samples
lab     <- read_otu_matrix("lab")        # 190 Lab samples (reference only)

site_levels <- c("Mandai Zoo", "Kranji Marshes", "Bird Paradise", "River Safari")
site_abbr   <- c("Mandai Zoo" = "MZ", "Kranji Marshes" = "KM",
                 "Bird Paradise" = "BP", "River Safari" = "RS", "Lab" = "Lab")

make_sd <- function(m) {
  tibble(Sample_ID = rownames(m)) |>
    left_join(select(meta, Sample_ID, Origin, Location, Lab_Taxon), by = "Sample_ID") |>
    mutate(
      depth     = as.numeric(rowSums(m)),
      log_depth = log10(depth),
      observed  = as.numeric(rowSums(m > 0)),
      shannon   = as.numeric(diversity(m, index = "shannon")),
      pielou    = if_else(observed > 1, shannon / log(observed), NA_real_)
    )
}

sd_w <- make_sd(wild) |> mutate(Location = factor(as.character(Location), levels = site_levels))
sd_l <- make_sd(lab)

site_n <- table(sd_w$Location)
stopifnot(
  identical(sd_w$Sample_ID, rownames(wild)),
  identical(sd_l$Sample_ID, rownames(lab)),
  all(sd_w$Origin == "Wild"), all(sd_l$Origin == "Lab"),
  !anyNA(sd_w$Location),
  all(as.numeric(site_n) == c(194, 88, 39, 24))
)
save_table(sd_w, "02_per_sample_alpha_depth_wild_locations")

# Group = site for Wild, "Lab" for reference
group_levels <- c(site_levels, "Lab")
group_n      <- c(as.numeric(site_n), nrow(sd_l))
names(group_n) <- group_levels
group_labels_short <- setNames(paste0(site_abbr[group_levels], "\n(n = ", group_n, ")"), group_levels)
group_labels_full  <- setNames(paste0(group_levels, " (n = ", group_n, ")"), group_levels)
group_labels_full["Lab"] <- paste0("Lab, reference (n = ", group_n["Lab"], ")")

fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, digits = 2, format = "g")))

# ---- 2. Genus composition (Fig. 5A) -----------------------------------------
comb  <- rbind(wild, lab)
sd_c  <- bind_rows(sd_w |> mutate(Group = as.character(Location)),
                   sd_l |> mutate(Group = "Lab")) |>
  mutate(Group = factor(Group, levels = group_levels))
stopifnot(identical(sd_c$Sample_ID, rownames(comb)))

genus_long <- as_tibble(comb / rowSums(comb), rownames = "Sample_ID") |>
  pivot_longer(-Sample_ID, names_to = "OTU", values_to = "rel") |>
  left_join(select(otu_map, OTU, Genus), by = "OTU") |>
  mutate(Genus_plot = lump_genus(Genus)) |>
  group_by(Sample_ID, Genus_plot) |>
  summarise(rel = sum(rel), .groups = "drop") |>
  left_join(select(sd_c, Sample_ID, Group), by = "Sample_ID")

comp_mean <- genus_long |>
  group_by(Group, Genus_plot) |>
  summarise(mean_rel = mean(rel), .groups = "drop")
save_table(comp_mean |> mutate(mean_percent = round(100 * mean_rel, 2)),
           "02_genus_mean_relative_abundance_by_location", paths$stats)

genus_label_expr <- as.expression(lapply(names(genus_colours), function(g) {
  lab_g <- pretty_genus(g)
  if (g %in% focal_genera) bquote(italic(.(lab_g))) else lab_g
}))

p_comp <- ggplot(comp_mean, aes(x = Group, y = mean_rel, fill = Genus_plot)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = 4.5, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  scale_fill_manual(values = genus_colours, breaks = names(genus_colours),
                    labels = genus_label_expr, drop = FALSE, name = "Genus") +
  scale_x_discrete(labels = site_abbr) +   # n per group is in the PCoA legend
  scale_y_continuous(labels = scales::percent, expand = c(0, 0), limits = c(0, 1.0001)) +
  labs(x = NULL, y = "Mean relative abundance") +
  theme(legend.key.size = unit(3.5, "mm"), legend.text = element_text(size = 8),
        axis.text.x = element_text(size = 8))

# Panel A = every genus as a heatmap (sites + Lab reference column). The top-8
# stacked bar is kept as a supplementary/slide figure.
heat_tab <- genus_mean_table(comb, sd_c$Group,
                             facet = factor(if_else(sd_c$Group == "Lab", "Ref.", "Wild sites"),
                                            levels = c("Wild sites", "Ref.")))
save_table(heat_tab, "02_all_genera_by_location", paths$stats)
p_heat <- genus_heatmap(heat_tab, x_labels = site_abbr,
                        breaks = c(1e-3, 0.1, 10), legend_width = 30)
save_fig(p_comp, "FigS_locations_top8_bars", width = 120, height = 110)

# ---- 3. Alpha diversity (Fig. 5B) - Wild sites only -------------------------
alpha_metrics <- c("observed", "shannon", "pielou")
alpha_names   <- c(observed = "Observed\nOTUs", shannon = "Shannon\ndiversity",
                   pielou = "Pielou's\nevenness")

kw_tab <- map_dfr(alpha_metrics, function(m) {
  t <- kruskal.test(sd_w[[m]] ~ sd_w$Location)
  tibble(metric = m, H = as.numeric(t$statistic), df = as.numeric(t$parameter), p = t$p.value)
}) |> mutate(p_BH = p.adjust(p, method = "BH"))
save_table(kw_tab, "02_alpha_kruskal_wallis", paths$stats)

alpha_medians <- sd_w |>
  pivot_longer(all_of(alpha_metrics), names_to = "metric", values_to = "value") |>
  filter(!is.na(value)) |>
  group_by(metric, Location) |>
  summarise(n = n(), median = median(value),
            q1 = as.numeric(quantile(value, 0.25)), q3 = as.numeric(quantile(value, 0.75)),
            .groups = "drop")
save_table(alpha_medians, "02_alpha_medians_by_location", paths$stats)

# Pairwise Wilcoxon between sites, BH within each metric (6 comparisons)
alpha_pairwise <- map_dfr(alpha_metrics, function(m) {
  pt <- pairwise.wilcox.test(sd_w[[m]], sd_w$Location, p.adjust.method = "BH", exact = FALSE)
  as.data.frame(as.table(pt$p.value)) |>
    filter(!is.na(Freq)) |>
    transmute(metric = m, group1 = as.character(Var2), group2 = as.character(Var1), p_BH = Freq)
})
save_table(alpha_pairwise, "02_alpha_pairwise_wilcoxon", paths$stats)

alpha_long <- sd_w |>
  select(Sample_ID, Location, all_of(alpha_metrics)) |>
  pivot_longer(all_of(alpha_metrics), names_to = "metric", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(metric = factor(alpha_names[metric], levels = alpha_names))

# Kruskal-Wallis p goes into the panel title (no text squeezed into the plot)
strip_labels <- setNames(
  paste0(alpha_names[kw_tab$metric], "\nK-W ", fmt_p(kw_tab$p_BH)),
  alpha_names[kw_tab$metric])

p_alpha <- ggplot(alpha_long, aes(x = Location, y = value)) +
  geom_jitter(aes(colour = Location), width = 0.2, height = 0, size = 0.6, alpha = 0.4) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.55, linewidth = 0.4) +
  facet_wrap(~ metric, scales = "free_y", labeller = as_labeller(strip_labels)) +
  scale_colour_manual(values = location_colours, guide = "none") +
  scale_x_discrete(labels = site_abbr) +
  labs(x = NULL, y = NULL) +
  theme(strip.text  = element_text(size = 8.5, lineheight = 0.9),
        strip.clip  = "off",
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

# ---- 4. Beta diversity - statistics on Wild sites only ----------------------
sd_w_df <- as.data.frame(sd_w)
d_w     <- vegdist(wild, method = "jaccard", binary = TRUE)

perm_loc_only <- adonis2(d_w ~ Location, data = sd_w_df, permutations = 999, by = "terms")
perm_loc      <- adonis2(d_w ~ log_depth + Location, data = sd_w_df,
                         permutations = 999, by = "margin")
disp_w        <- betadisper(d_w, sd_w_df$Location)
disp_w_test   <- permutest(disp_w, permutations = 999)

tidy_adonis <- function(a, model) {
  as.data.frame(a) |> rownames_to_column("term") |> mutate(model = model, .before = 1)
}
save_table(bind_rows(tidy_adonis(perm_loc_only, "Location only"),
                     tidy_adonis(perm_loc, "log10(depth) + Location, marginal")),
           "02_permanova_locations", paths$stats)

disp_out <- tibble(
  Location = levels(sd_w$Location),
  n = as.numeric(site_n),
  mean_distance_to_centroid = as.numeric(tapply(disp_w$distances, disp_w$group, mean)),
  F = disp_w_test$tab[1, "F"],
  p = disp_w_test$tab[1, "Pr(>F)"]
)
save_table(disp_out, "02_permdisp_locations", paths$stats)

# Pairwise PERMANOVA + PERMDISP between sites (6 pairs, BH-adjusted)
dm    <- as.matrix(d_w)
pairs <- combn(site_levels, 2, simplify = FALSE)
pairwise_beta <- map_dfr(pairs, function(pr) {
  idx <- sd_w$Location %in% pr
  sub <- as.dist(dm[idx, idx])
  dat <- as.data.frame(sd_w[idx, ]) |> mutate(Location = droplevels(Location))
  a   <- adonis2(sub ~ log_depth + Location, data = dat, permutations = 999, by = "margin")
  dsp <- permutest(betadisper(sub, dat$Location), permutations = 999)
  tibble(group1 = pr[1], group2 = pr[2],
         n1 = sum(dat$Location == pr[1]), n2 = sum(dat$Location == pr[2]),
         R2 = a["Location", "R2"], F = a["Location", "F"], p = a["Location", "Pr(>F)"],
         permdisp_F = dsp$tab[1, "F"], permdisp_p = dsp$tab[1, "Pr(>F)"])
}) |>
  mutate(p_BH = p.adjust(p, method = "BH"),
         permdisp_p_BH = p.adjust(permdisp_p, method = "BH"))
save_table(pairwise_beta, "02_pairwise_permanova_locations", paths$stats)

# ---- 5. PCoA (Fig. 5C): Wild sites + Lab reference in one ordination --------
d_c  <- vegdist(comb, method = "jaccard", binary = TRUE)
pc   <- cmdscale(d_c, k = 2, eig = TRUE)
ve_c <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
pcoa_c <- sd_c |> mutate(PCo1 = pc$points[Sample_ID, 1], PCo2 = pc$points[Sample_ID, 2])

site_shapes <- c("Mandai Zoo" = 16, "Kranji Marshes" = 17, "Bird Paradise" = 15,
                 "River Safari" = 18, "Lab" = 1)

R2_l <- perm_loc["Location", "R2"]
F_l  <- perm_loc["Location", "F"]
p_l  <- perm_loc["Location", "Pr(>F)"]
p_dl <- disp_w_test$tab[1, "Pr(>F)"]

pcoa_plot <- function(df, ve, show_lab = TRUE, subtitle = NULL) {
  sites_df <- filter(df, Group != "Lab")
  p <- ggplot(mapping = aes(x = PCo1, y = PCo2, colour = Group, shape = Group))
  if (show_lab) {
    p <- p + geom_point(data = filter(df, Group == "Lab"), size = 1, alpha = 0.35)
  }
  p +
    geom_point(data = sites_df, size = 1.3, alpha = 0.7) +
    stat_ellipse(data = sites_df, aes(shape = NULL), type = "t", level = 0.95, linewidth = 0.8) +
    scale_colour_manual(values = location_colours, labels = group_labels_full, name = NULL) +
    scale_shape_manual(values = site_shapes, labels = group_labels_full, name = NULL) +
    labs(x = sprintf("PCoA 1 (%.1f%%)", ve[1]), y = sprintf("PCoA 2 (%.1f%%)", ve[2]),
         subtitle = subtitle) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(plot.subtitle = element_text(size = 8), legend.position = "bottom",
          legend.text = element_text(size = 8)) +
      coord_fixed()   # equal axis scaling: distances on the plot match Jaccard distances
}

stats_subtitle <- sprintf(
  "Jaccard | Wild sites only\nPERMANOVA (controlling for depth): R² = %.3f, F = %.1f, %s | PERMDISP: %s",
  R2_l, F_l, fmt_p(p_l), fmt_p(p_dl))

p_pcoa <- pcoa_plot(pcoa_c, ve_c, show_lab = TRUE, subtitle = stats_subtitle)

# Supplementary: Wild-only ordination (the space the statistics were run in)
pc_w   <- cmdscale(d_w, k = 2, eig = TRUE)
ve_w   <- 100 * pc_w$eig[1:2] / sum(pc_w$eig[pc_w$eig > 0])
pcoa_w <- sd_w |>
  mutate(Group = factor(as.character(Location), levels = group_levels),
         PCo1 = pc_w$points[Sample_ID, 1], PCo2 = pc_w$points[Sample_ID, 2])
p_pcoa_w <- pcoa_plot(pcoa_w, ve_w, show_lab = FALSE, subtitle = stats_subtitle)
save_fig(p_pcoa_w, "FigS_locations_wild_only_pcoa", width = 180, height = 120)

# ---- 6. Assemble Figure 5 ---------------------------------------------------
fig5 <- (free(p_heat) | (p_alpha / p_pcoa + plot_layout(heights = c(1, 1.25)))) +
  plot_layout(widths = c(1.35, 2.1)) +
  plot_annotation(tag_levels = "A")
save_fig(fig5, "Fig5_locations", width = 230, height = 250)

# ---- 7. Sensitivity: Wild location samples with >= 1,000 classified reads ---
keep   <- sd_w$depth >= 1000
sd_ws  <- sd_w[keep, ] |> mutate(Location = droplevels(Location))
d_ws   <- vegdist(wild[keep, ], method = "jaccard", binary = TRUE)
perm_s <- adonis2(d_ws ~ log_depth + Location, data = as.data.frame(sd_ws),
                  permutations = 999, by = "margin")
disp_s <- permutest(betadisper(d_ws, sd_ws$Location), permutations = 999)
kw_s   <- map_dfr(alpha_metrics, function(m) {
  t <- kruskal.test(sd_ws[[m]] ~ sd_ws$Location)
  tibble(metric = m, H = as.numeric(t$statistic), df = as.numeric(t$parameter), p = t$p.value)
}) |> mutate(p_BH = p.adjust(p, method = "BH"))

save_table(bind_rows(
  tidy_adonis(perm_s, "SENSITIVITY >=1000 reads: log10(depth) + Location, marginal"),
  tibble(model = "SENSITIVITY >=1000 reads: PERMDISP", term = "Location",
         F = disp_s$tab[1, "F"], `Pr(>F)` = disp_s$tab[1, "Pr(>F)"])
), "02_sensitivity_min1000_beta", paths$stats)
save_table(kw_s, "02_sensitivity_min1000_alpha", paths$stats)

# ---- 8. Plain-text summary --------------------------------------------------
med_line <- function(m, nm) {
  r <- filter(alpha_medians, metric == m)
  k <- filter(kw_tab, metric == m)
  paste0(sprintf("%-15s ", nm),
         paste(sprintf("%s %.2f (%.2f-%.2f)", site_abbr[as.character(r$Location)], r$median, r$q1, r$q3),
               collapse = "; "),
         sprintf(" | Kruskal-Wallis H = %.2f, df = %d, %s (BH %s)", k$H, as.integer(k$df), fmt_p(k$p), fmt_p(k$p_BH)))
}

sig_pairs <- function(df, pcol) {
  s <- df[df[[pcol]] < 0.05, ]
  if (nrow(s) == 0) return("  none (all BH-adjusted p >= 0.05)")
  sprintf("  %s vs %s: BH %s", site_abbr[s$group1], site_abbr[s$group2], fmt_p(s[[pcol]]))
}

summary_lines <- c(
  "02_locations.R - results summary (Wild samples with confirmed locations)",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  do.call(sprintf, c(list("n: MZ %d, KM %d, BP %d, RS %d (total %d)"), as.list(as.integer(c(site_n, sum(site_n)))))),
  "",
  "SEQUENCING DEPTH (classified reads), median by site",
  paste0("  ", paste(sprintf("%s %.0f", site_abbr[site_levels],
                             tapply(sd_w$depth, sd_w$Location, median)), collapse = "; ")),
  sprintf("  Kruskal-Wallis: %s", fmt_p(kruskal.test(depth ~ Location, data = sd_w)$p.value)),
  "",
  "ALPHA DIVERSITY, median (IQR) by site",
  med_line("observed", "Observed OTUs"),
  med_line("shannon",  "Shannon"),
  med_line("pielou",   "Pielou"),
  "Significant pairwise Wilcoxon (BH within metric):",
  unlist(lapply(alpha_metrics, function(m) {
    c(paste0(" ", m, ":"), sig_pairs(filter(alpha_pairwise, metric == m), "p_BH"))
  })),
  "",
  "BETA DIVERSITY (Jaccard presence/absence, Wild sites only, 999 permutations)",
  sprintf("PERMANOVA Location only:            R2 = %.3f, F = %.2f, %s",
          perm_loc_only["Location", "R2"], perm_loc_only["Location", "F"], fmt_p(perm_loc_only["Location", "Pr(>F)"])),
  sprintf("PERMANOVA Location | log10(depth):  R2 = %.3f, F = %.2f, %s", R2_l, F_l, fmt_p(p_l)),
  sprintf("PERMANOVA log10(depth) | Location:  R2 = %.3f, F = %.2f, %s",
          perm_loc["log_depth", "R2"], perm_loc["log_depth", "F"], fmt_p(perm_loc["log_depth", "Pr(>F)"])),
  sprintf("PERMDISP: F = %.2f, %s; distance to centroid %s",
          disp_w_test$tab[1, "F"], fmt_p(p_dl),
          paste(sprintf("%s %.3f", site_abbr[disp_out$Location], disp_out$mean_distance_to_centroid), collapse = ", ")),
  "Pairwise PERMANOVA (depth-adjusted, BH over 6 pairs):",
  sprintf("  %s vs %s: R2 = %.3f, F = %.2f, BH %s | PERMDISP BH %s",
          site_abbr[pairwise_beta$group1], site_abbr[pairwise_beta$group2],
          pairwise_beta$R2, pairwise_beta$F, fmt_p(pairwise_beta$p_BH), fmt_p(pairwise_beta$permdisp_p_BH)),
  sprintf("PCoA (Wild + Lab reference): PCo1 %.1f%%, PCo2 %.1f%%", ve_c[1], ve_c[2]),
  sprintf("PCoA (Wild only, supplementary): PCo1 %.1f%%, PCo2 %.1f%%", ve_w[1], ve_w[2]),
  "",
  sprintf("SENSITIVITY (>= 1,000 reads: %s)",
          paste(sprintf("%s %d", site_abbr[levels(sd_ws$Location)], as.numeric(table(sd_ws$Location))), collapse = ", ")),
  sprintf("PERMANOVA Location | log10(depth):  R2 = %.3f, F = %.2f, %s",
          perm_s["Location", "R2"], perm_s["Location", "F"], fmt_p(perm_s["Location", "Pr(>F)"])),
  sprintf("PERMDISP: F = %.2f, %s", disp_s$tab[1, "F"], fmt_p(disp_s$tab[1, "Pr(>F)"])),
  sprintf("Kruskal-Wallis BH: observed %s, Shannon %s, Pielou %s",
          fmt_p(kw_s$p_BH[1]), fmt_p(kw_s$p_BH[2]), fmt_p(kw_s$p_BH[3])),
  "",
  "TOP GENERA, mean relative abundance (%) by group (Lab = reference)",
  comp_mean |>
    mutate(pct = 100 * mean_rel) |>
    select(Group, Genus_plot, pct) |>
    pivot_wider(names_from = Group, values_from = pct) |>
    mutate(line = sprintf("  %-24s MZ %5.1f | KM %5.1f | BP %5.1f | RS %5.1f | Lab %5.1f",
                          pretty_genus(Genus_plot), `Mandai Zoo`, `Kranji Marshes`,
                          `Bird Paradise`, `River Safari`, Lab)) |>
    pull(line)
)
dir.create(paths$stats, showWarnings = FALSE, recursive = TRUE)
writeLines(summary_lines, file.path(paths$stats, "02_results_summary.txt"))
cat(summary_lines, sep = "\n")
