# =============================================================================
# 04_qc_supplementary.R - three supplementary figures
#
#   A. Wolbachia: relative abundance + prevalence in Wild, Lab Ae. albopictus
#      and Lab Culex (Wild pooled here; wild host species are analysed in 11)
#   B. Negative controls: genus composition of the 6 H12 controls next to the
#      mean Wild and Lab profiles, with read counts
#   C. Depth vs richness: classified reads vs observed OTUs (no-rarefaction
#      justification, Methods 2.11)
#
# Outputs
#   output/figures/FigS_wolbachia.pdf/.png
#   output/figures/FigS_negative_controls.pdf/.png
#   output/figures/FigS_depth_vs_richness.pdf/.png
#   output/stats/04_*.tsv, output/stats/04_results_summary.txt
# =============================================================================

source(here::here("R", "00_setup.R"))
check_inputs()

meta    <- read_metadata()
otu_map <- read_otu_map()
bio     <- read_otu_matrix("biological")
ctrl    <- read_otu_matrix("controls")

fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, digits = 2, format = "g")))

genus_of <- function(m) otu_map$Genus[match(colnames(m), otu_map$OTU)]

sd <- tibble(Sample_ID = rownames(bio)) |>
  left_join(select(meta, Sample_ID, Origin, Lab_Taxon), by = "Sample_ID") |>
  mutate(Origin   = droplevels(Origin),
         depth    = as.numeric(rowSums(bio)),
         observed = as.numeric(rowSums(bio > 0)),
         Group    = case_when(Origin == "Wild" ~ "Wild",
                              Lab_Taxon == "Ae. albopictus" ~ "Lab Ae. albopictus",
                              Lab_Taxon == "Culex" ~ "Lab Culex"),
         Group    = factor(Group, levels = c("Wild", "Lab Ae. albopictus", "Lab Culex")))
stopifnot(identical(sd$Sample_ID, rownames(bio)), !anyNA(sd$Group))

# =============================================================================
# A. Wolbachia
# =============================================================================
wol_otus <- otu_map$OTU[otu_map$Genus == "Wolbachia"]
stopifnot(length(wol_otus) > 0)
sd <- sd |>
  mutate(wol_reads = as.numeric(rowSums(bio[, wol_otus, drop = FALSE])),
         wol_pct   = 100 * wol_reads / depth,
         wol_pres  = wol_reads > 0)

wol_summary <- sd |>
  group_by(Group) |>
  summarise(n = n(),
            prevalence_pct = 100 * mean(wol_pres),
            mean_pct = mean(wol_pct), median_pct = median(wol_pct),
            q1 = as.numeric(quantile(wol_pct, 0.25)), q3 = as.numeric(quantile(wol_pct, 0.75)),
            .groups = "drop")
save_table(wol_summary, "04_wolbachia_by_group", paths$stats)

wol_kw <- kruskal.test(wol_pct ~ Group, data = sd)
wol_pw <- pairwise.wilcox.test(sd$wol_pct, sd$Group, p.adjust.method = "BH", exact = FALSE)
wol_pw_tab <- as.data.frame(as.table(wol_pw$p.value)) |>
  filter(!is.na(Freq)) |>
  transmute(group1 = as.character(Var2), group2 = as.character(Var1), p_BH = Freq)
save_table(wol_pw_tab, "04_wolbachia_pairwise", paths$stats)

group_cols <- c("Wild" = "#1baf7a", "Lab Ae. albopictus" = "#4a3aa7", "Lab Culex" = "#eb6834")
group_labs <- c("Wild" = "Wild",
                "Lab Ae. albopictus" = "Lab\nAe. albopictus",
                "Lab Culex" = "Lab\nCulex")
group_labs_n <- setNames(paste0(group_labs[as.character(wol_summary$Group)], "\n(n = ", wol_summary$n, ")"),
                         as.character(wol_summary$Group))

p_wol <- ggplot(sd, aes(x = Group, y = wol_pct)) +
  geom_jitter(aes(colour = Group), width = 0.22, height = 0, size = 0.8, alpha = 0.4) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5, linewidth = 0.4) +
  geom_text(data = wol_summary,
            aes(x = Group, y = 104, label = sprintf("%.0f%% carry", prevalence_pct)),
            size = 3, inherit.aes = FALSE) +
  scale_colour_manual(values = group_cols, guide = "none") +
  scale_x_discrete(labels = as.expression(list(
    bquote(atop("Wild", "(n = " * .(wol_summary$n[1]) * ")")),
    bquote(atop("Lab" ~ italic("Ae. albopictus"), "(n = " * .(wol_summary$n[2]) * ")")),
    bquote(atop("Lab" ~ italic("Culex"), "(n = " * .(wol_summary$n[3]) * ")"))))) +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25)) +
  labs(x = NULL, y = expression(italic(Wolbachia)~"relative abundance (%)"),
       subtitle = paste0("Kruskal-Wallis ", fmt_p(wol_kw$p.value))) +
  theme(plot.subtitle = element_text(size = 8))
save_fig(p_wol, "FigS_wolbachia", width = 110, height = 110)

# =============================================================================
# B. Negative controls
# =============================================================================
plate_of <- function(id) str_extract(id, "(AG|ZX|BT)_P0[12]")

lumped_profile <- function(m) {
  rel <- m / rowSums(m)
  as_tibble(rel, rownames = "Sample_ID") |>
    pivot_longer(-Sample_ID, names_to = "OTU", values_to = "rel") |>
    left_join(select(otu_map, OTU, Genus), by = "OTU") |>
    mutate(Genus_plot = lump_genus(Genus)) |>
    group_by(Sample_ID, Genus_plot) |>
    summarise(rel = sum(rel), .groups = "drop")
}

ctrl_prof <- lumped_profile(ctrl) |>
  mutate(Bar = paste0("NC ", plate_of(Sample_ID)))
bio_mean  <- lumped_profile(bio) |>
  left_join(select(sd, Sample_ID, Origin), by = "Sample_ID") |>
  group_by(Origin, Genus_plot) |>
  summarise(rel = mean(rel), .groups = "drop") |>
  mutate(Bar = paste0(Origin, " (mean)"))

ctrl_reads <- tibble(Sample_ID = rownames(ctrl),
                     Bar = paste0("NC ", plate_of(rownames(ctrl))),
                     reads = as.numeric(rowSums(ctrl)),
                     otus = as.numeric(rowSums(ctrl > 0)))
stopifnot(nrow(ctrl_reads) == 6, !anyNA(ctrl_reads$Bar))
save_table(ctrl_reads, "04_negative_control_reads", paths$stats)

bar_levels <- c(sort(unique(ctrl_prof$Bar)), "Wild (mean)", "Lab (mean)")
nc_df <- bind_rows(select(ctrl_prof, Bar, Genus_plot, rel),
                   select(bio_mean,  Bar, Genus_plot, rel)) |>
  mutate(Bar = factor(Bar, levels = bar_levels))

read_lab <- bind_rows(
  ctrl_reads |> transmute(Bar, label = paste0(scales::comma(reads), "\nreads")),
  tibble(Bar = c("Wild (mean)", "Lab (mean)"),
         label = c(paste0("median\n", scales::comma(median(sd$depth[sd$Origin == "Wild"]))),
                   paste0("median\n", scales::comma(median(sd$depth[sd$Origin == "Lab"])))))
) |> mutate(Bar = factor(Bar, levels = bar_levels))

genus_label_expr <- as.expression(lapply(names(genus_colours), function(g) {
  lab_g <- pretty_genus(g)
  if (g %in% focal_genera) bquote(italic(.(lab_g))) else lab_g
}))

p_nc <- ggplot(nc_df, aes(x = Bar, y = rel, fill = Genus_plot)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_vline(xintercept = 6.5, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_text(data = read_lab, aes(x = Bar, y = 1.08, label = label),
            inherit.aes = FALSE, size = 2.5, lineheight = 0.85) +
  scale_fill_manual(values = genus_colours, breaks = names(genus_colours),
                    labels = genus_label_expr, drop = FALSE, name = "Genus") +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.25),
                     limits = c(0, 1.15), expand = c(0, 0)) +
  labs(x = NULL, y = "Relative abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.key.size = unit(3.5, "mm"), legend.text = element_text(size = 8))
save_fig(p_nc, "FigS_negative_controls", width = 170, height = 110)

# Which genera appear in the controls at all?
ctrl_gen <- t(rowsum(t(ctrl), group = genus_of(ctrl)))   # controls x genera
ctrl_genera <- tibble(Genus = colnames(ctrl_gen),
                      reads = as.numeric(colSums(ctrl_gen)),
                      n_controls_present = as.numeric(colSums(ctrl_gen > 0))) |>
  filter(reads > 0) |>
  arrange(desc(reads))
save_table(ctrl_genera, "04_negative_control_genera", paths$stats)

# =============================================================================
# C. Depth vs richness
# =============================================================================
rho_all <- cor.test(sd$depth, sd$observed, method = "spearman", exact = FALSE)
rho_grp <- sd |>
  group_by(Origin) |>
  summarise(rho = cor(depth, observed, method = "spearman"),
            p = cor.test(depth, observed, method = "spearman", exact = FALSE)$p.value,
            n = n(), .groups = "drop")
save_table(bind_rows(tibble(Origin = "All", rho = unname(rho_all$estimate),
                            p = rho_all$p.value, n = nrow(sd)),
                     mutate(rho_grp, Origin = as.character(Origin))),
           "04_depth_vs_richness_spearman", paths$stats)

rho_txt <- sprintf("Spearman rho: all %.2f | Wild %.2f | Lab %.2f",
                   rho_all$estimate, rho_grp$rho[rho_grp$Origin == "Wild"],
                   rho_grp$rho[rho_grp$Origin == "Lab"])

p_depth <- ggplot(sd, aes(x = depth, y = observed, colour = Origin)) +
  geom_jitter(width = 0, height = 0.15, size = 1, alpha = 0.5) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.9) +
  geom_vline(xintercept = 1000, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  scale_x_log10(labels = scales::comma) +
  scale_colour_manual(values = origin_colours, name = NULL) +
  labs(x = "Classified reads per sample (log scale)", y = "Observed OTUs",
       subtitle = paste0(rho_txt, " | dashed line: 1,000-read sensitivity cut-off")) +
  theme(plot.subtitle = element_text(size = 8), legend.position = "top")
save_fig(p_depth, "FigS_depth_vs_richness", width = 140, height = 110)

# =============================================================================
# Summary
# =============================================================================
summary_lines <- c(
  "04_qc_supplementary.R - results summary",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "WOLBACHIA (relative abundance % of classified reads)",
  sprintf("  %-20s n = %3d | carry %5.1f%% | mean %5.2f%% | median %5.2f%% (IQR %.2f-%.2f)",
          as.character(wol_summary$Group), wol_summary$n, wol_summary$prevalence_pct,
          wol_summary$mean_pct, wol_summary$median_pct, wol_summary$q1, wol_summary$q3),
  sprintf("  Kruskal-Wallis H = %.2f, df = %d, %s", wol_kw$statistic, as.integer(wol_kw$parameter), fmt_p(wol_kw$p.value)),
  sprintf("  %s vs %s: BH %s", wol_pw_tab$group1, wol_pw_tab$group2, fmt_p(wol_pw_tab$p_BH)),
  "",
  "NEGATIVE CONTROLS",
  sprintf("  %-10s %6s reads, %2d OTUs", ctrl_reads$Bar, scales::comma(ctrl_reads$reads),
          as.integer(ctrl_reads$otus)),
  "  Genera detected in controls (reads, number of controls):",
  sprintf("    %-26s %6s reads in %d control(s)", pretty_genus(ctrl_genera$Genus),
          scales::comma(ctrl_genera$reads), as.integer(ctrl_genera$n_controls_present)),
  "",
  "DEPTH VS RICHNESS",
  sprintf("  All: rho = %.3f, %s (n = %d)", rho_all$estimate, fmt_p(rho_all$p.value), nrow(sd)),
  sprintf("  %s: rho = %.3f, %s (n = %d)", as.character(rho_grp$Origin), rho_grp$rho, fmt_p(rho_grp$p), rho_grp$n)
)
writeLines(summary_lines, file.path(paths$stats, "04_results_summary.txt"))
cat(summary_lines, sep = "\n")
