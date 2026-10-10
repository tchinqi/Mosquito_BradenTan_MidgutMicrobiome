# =============================================================================
# 11_host_species.R - host species (COI) and the midgut microbiome
#
# Needs data/raw/species/wild_host_species_COI.tsv  (run 11a_build_host_species.R)
# Lab species come from the sample IDs (BT_ALB / BT_CLX), wild species from COI.
#
# PART A (Results 3.5.2) - does rearing origin change the microbiome WITHIN a
#   species?  Species present in both origins only: Ae. albopictus and Culex spp.
#   PERMANOVA  Jaccard ~ log10(depth) + Origin + Species        (marginal; main effects)
#              Jaccard ~ log10(depth) + Origin * Species        (sequential; interaction)
#   Pairwise   wild vs lab within each species; Ae. alb vs Culex within each origin
#              (PERMANOVA with depth, PERMDISP, Wilcoxon alpha; BH over the 4 pairs)
#   Core (90%) per Origin x Species group
#
# PART B (Results 3.5.3) - do WILD host species differ?  Wild only:
#   Ae. albopictus vs Mansonia indiana vs Culex spp. (Anopheles, Aedeomyia,
#   unidentified and pending samples excluded: too few / unknown).
#   Kruskal-Wallis + pairwise Wilcoxon; PERMANOVA (depth) + PERMDISP + pairwise;
#   PERMANOVA with collection site as an extra covariate (site-adjusted);
#   Kranji Marshes only: Ae. albopictus vs Mansonia indiana (site held constant).
#   Lab mosquitoes appear in the Fig. 7 heatmap/PCoA as a grey reference only.
#
# Outputs
#   output/figures/Fig6_origin_x_species.pdf/.png
#   output/figures/Fig7_wild_host_species.pdf/.png
#   output/stats/11_*.tsv, output/stats/11_results_summary.txt
# Caveats: wild Culex is a mix of species (mostly genus-level COI IDs) and the lab
# Culex colony's species is unconfirmed; Culex groups are small (n ~ 25).
# =============================================================================

source(here::here("R", "00_setup.R"))
check_inputs()

host_file <- file.path(paths$species, "wild_host_species_COI.tsv")
stopifnot("Run 11a_build_host_species.R first" = file.exists(host_file))

fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, digits = 2, format = "g")))
sp_levels  <- c("Ae. albopictus", "Mansonia indiana", "Culex spp.")
sp_short   <- c("Ae. albopictus" = "Ae. alb.", "Mansonia indiana" = "Mansonia", "Culex spp." = "Culex")
sp_colours <- c("Ae. albopictus" = "#4a3aa7", "Mansonia indiana" = "#008300", "Culex spp." = "#eb6834")

# ---- 1. Load + per-sample metrics -------------------------------------------
meta <- read_metadata()
host <- read_tsv(host_file, show_col_types = FALSE) |> select(Sample_ID, Host_group)
bio  <- read_otu_matrix("biological")

sd <- tibble(Sample_ID = rownames(bio)) |>
  left_join(select(meta, Sample_ID, Origin, Location, Lab_Taxon), by = "Sample_ID") |>
  left_join(host, by = "Sample_ID") |>
  mutate(Origin   = droplevels(Origin),
         Species  = case_when(Origin == "Lab" & Lab_Taxon == "Ae. albopictus" ~ "Ae. albopictus",
                              Origin == "Lab" & Lab_Taxon == "Culex"          ~ "Culex spp.",
                              Origin == "Wild"                                ~ Host_group,
                              TRUE                                            ~ NA_character_),
         depth     = as.numeric(rowSums(bio)),
         log_depth = log10(depth),
         observed  = as.numeric(rowSums(bio > 0)),
         shannon   = as.numeric(diversity(bio, index = "shannon")),
         pielou    = if_else(observed > 1, shannon / log(observed), NA_real_))
stopifnot(identical(sd$Sample_ID, rownames(bio)),
          sum(!is.na(sd$Host_group)) == 380)        # every wild sample has a row
save_table(count(sd, Origin, Species), "11_species_counts", paths$stats)
save_table(sd |> filter(Origin == "Wild") |> count(Species, Location) |>
             pivot_wider(names_from = Location, values_from = n, values_fill = 0),
           "11_wild_species_by_site", paths$stats)

alpha_metrics <- c("observed", "shannon", "pielou")
alpha_names   <- c(observed = "Observed\nOTUs", shannon = "Shannon\ndiversity", pielou = "Pielou's\nevenness")

# ---- 2. Shared helpers ------------------------------------------------------
jac <- function(ids) vegdist(bio[ids, , drop = FALSE], method = "jaccard", binary = TRUE)

# PERMANOVA (depth + group, marginal) + PERMDISP between two levels of a column
pair_beta <- function(df, col, a, b) {
  d2 <- as.data.frame(df[df[[col]] %in% c(a, b), ])
  d2$G <- factor(as.character(d2[[col]]), levels = c(a, b))
  dd <- jac(d2$Sample_ID)
  ad <- adonis2(dd ~ log_depth + G, data = d2, permutations = 999, by = "margin")
  ds <- permutest(betadisper(dd, d2$G), permutations = 999)
  tibble(group1 = a, group2 = b, n1 = sum(d2$G == a), n2 = sum(d2$G == b),
         R2 = ad["G", "R2"], F = ad["G", "F"], p = ad["G", "Pr(>F)"],
         permdisp_F = ds$tab[1, "F"], permdisp_p = ds$tab[1, "Pr(>F)"])
}
# Wilcoxon on every alpha metric between two levels of a column
pair_alpha <- function(df, col, a, b) {
  map_dfr(alpha_metrics, function(m) {
    x <- df[[m]][df[[col]] == a]; y <- df[[m]][df[[col]] == b]
    x <- x[!is.na(x)]; y <- y[!is.na(y)]
    t <- wilcox.test(x, y, exact = FALSE)
    tibble(metric = m, group1 = a, group2 = b, n1 = length(x), n2 = length(y),
           median1 = median(x), median2 = median(y), W = as.numeric(t$statistic), p = t$p.value)
  })
}
# Prevalence per genus per group (Unclassified excluded) and core at 90%
core_by <- function(df, col) {
  g <- genus_counts(bio[df$Sample_ID, , drop = FALSE])
  g <- g[, colnames(g) != "Unclassified_at_genus", drop = FALSE]
  map_dfr(levels(df[[col]]), function(lv) {
    ids <- df$Sample_ID[df[[col]] == lv]
    tibble(Group = lv, n = length(ids), Genus = colnames(g),
           prevalence = unname(100 * colMeans(g[ids, , drop = FALSE] > 0)))
  }) |> mutate(core_90 = prevalence >= 90 - 1e-9)
}
core_lines <- function(ct) {
  ct |> filter(core_90) |> arrange(Group, desc(prevalence)) |>
    group_by(Group) |>
    summarise(n = first(n), k = n(),
              genera = paste(sprintf("%s %.0f%%", pretty_genus(Genus), prevalence), collapse = ", "),
              .groups = "drop") |>
    mutate(line = sprintf("  %-26s (n = %3d): %d core - %s", Group, n, k, genera)) |>
    pull(line)
}
tidy_adonis <- function(a, model) as.data.frame(a) |> rownames_to_column("term") |> mutate(model = model, .before = 1)

# =============================================================================
# PART A - Origin x Species (Ae. albopictus and Culex spp., wild + lab)
# =============================================================================
A_levels <- c("Wild Ae. albopictus", "Lab Ae. albopictus", "Wild Culex spp.", "Lab Culex spp.")
A <- sd |>
  filter(Species %in% c("Ae. albopictus", "Culex spp.")) |>
  mutate(Species = factor(Species, levels = c("Ae. albopictus", "Culex spp.")),
         Group   = factor(paste(Origin, Species), levels = A_levels))
A_n <- table(A$Group)
A_df <- as.data.frame(A)
dA <- jac(A$Sample_ID)

pA_main <- adonis2(dA ~ log_depth + Origin + Species, data = A_df, permutations = 999, by = "margin")
pA_int  <- adonis2(dA ~ log_depth + Origin * Species, data = A_df, permutations = 999, by = "terms")
dispA   <- permutest(betadisper(dA, A_df$Group), permutations = 999)
save_table(bind_rows(tidy_adonis(pA_main, "depth + Origin + Species, marginal"),
                     tidy_adonis(pA_int,  "depth + Origin * Species, sequential")),
           "11A_permanova_origin_x_species", paths$stats)

A_pairs <- list(c("Wild Ae. albopictus", "Lab Ae. albopictus"), c("Wild Culex spp.", "Lab Culex spp."),
                c("Wild Ae. albopictus", "Wild Culex spp."),    c("Lab Ae. albopictus", "Lab Culex spp."))
A_beta  <- map_dfr(A_pairs, ~ pair_beta(A, "Group", .x[1], .x[2])) |>
  mutate(p_BH = p.adjust(p, "BH"), permdisp_p_BH = p.adjust(permdisp_p, "BH"))
A_alpha <- map_dfr(A_pairs, ~ pair_alpha(A, "Group", .x[1], .x[2])) |>
  group_by(metric) |> mutate(p_BH = p.adjust(p, "BH")) |> ungroup()
save_table(A_beta,  "11A_pairwise_beta", paths$stats)
save_table(A_alpha, "11A_pairwise_alpha", paths$stats)

A_core <- core_by(A, "Group")
save_table(A_core |> filter(prevalence >= 50), "11A_core_prevalence", paths$stats)

# ---- Fig. 6 -----------------------------------------------------------------
A_lab <- setNames(as.expression(lapply(A_levels, function(g) {          # heatmap: origin + n
  bquote(atop(.(word(g, 1)), "(" * .(as.integer(A_n[[g]])) * ")"))
})), A_levels)
A_lab_box <- setNames(as.expression(lapply(A_levels, function(g) {      # boxplots: origin over species
  word(g, 1)
})), A_levels)
A_heat_tab <- genus_mean_table(bio[A$Sample_ID, ], A$Group,
                               facet = factor(sp_short[as.character(A$Species)], levels = c("Ae. alb.", "Culex")))
save_table(A_heat_tab, "11A_all_genera_by_origin_x_species", paths$stats)
p6_heat <- genus_heatmap(A_heat_tab, x_labels = function(b) A_lab[b],
                         breaks = c(1e-3, 0.1, 10), legend_width = 30) +
  theme(strip.text = element_text(face = "bold.italic", size = 8.5))

A_long <- A |> select(Sample_ID, Group, Species, all_of(alpha_metrics)) |>
  pivot_longer(all_of(alpha_metrics), names_to = "metric", values_to = "value") |>
  filter(!is.na(value)) |> mutate(metric = factor(alpha_names[metric], levels = alpha_names))
p6_alpha <- ggplot(A_long, aes(x = Group, y = value)) +
  geom_jitter(aes(colour = Species), width = 0.2, height = 0, size = 0.6, alpha = 0.4) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.55, linewidth = 0.4) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_colour_manual(values = sp_colours, guide = "none") +
  scale_x_discrete(labels = function(b) A_lab_box[b]) +
  labs(x = NULL, y = NULL) +
  theme(strip.text = element_text(size = 8.5, lineheight = 0.9), strip.clip = "off",
        axis.text.x = element_text(size = 7))

pcA <- cmdscale(dA, k = 2, eig = TRUE)
veA <- 100 * pcA$eig[1:2] / sum(pcA$eig[pcA$eig > 0])
A_pc <- A |> mutate(PCo1 = pcA$points[Sample_ID, 1], PCo2 = pcA$points[Sample_ID, 2])
p6_pcoa <- ggplot(A_pc, aes(x = PCo1, y = PCo2, colour = Species)) +
  geom_point(aes(shape = Origin), size = 1.2, alpha = 0.6) +
  stat_ellipse(aes(group = Group, linetype = Origin), type = "t", level = 0.95, linewidth = 0.7) +
  scale_colour_manual(values = sp_colours, name = NULL, breaks = c("Ae. albopictus", "Culex spp."),
                      labels = c(expression(italic("Ae. albopictus")), expression(italic("Culex") ~ "spp."))) +
  scale_shape_manual(values = c(Wild = 16, Lab = 1), name = NULL) +
  scale_linetype_manual(values = c(Wild = "solid", Lab = "dashed"), name = NULL) +
  guides(linetype = guide_legend(title = "Origin", override.aes = list(colour = "black"),
                                 keywidth = unit(10, "mm")),
         shape    = guide_legend(title = "Origin"),
         colour   = guide_legend(title = "Species", override.aes = list(linetype = "solid"))) +
  labs(x = sprintf("PCoA 1 (%.1f%%)", veA[1]), y = sprintf("PCoA 2 (%.1f%%)", veA[2]),
       subtitle = sprintf("Jaccard | PERMANOVA (depth-adjusted): Origin R² = %.3f, %s; Species R² = %.3f, %s;\nOrigin × Species R² = %.3f, %s | PERMDISP (4 groups): %s",
                          pA_main["Origin", "R2"], fmt_p(pA_main["Origin", "Pr(>F)"]),
                          pA_main["Species", "R2"], fmt_p(pA_main["Species", "Pr(>F)"]),
                          pA_int["Origin:Species", "R2"], fmt_p(pA_int["Origin:Species", "Pr(>F)"]),
                          fmt_p(dispA$tab[1, "Pr(>F)"]))) +
  theme(plot.subtitle = element_text(size = 7.5), legend.position = "bottom",
        legend.text = element_text(size = 8)) +
  coord_fixed()

fig6 <- (free(p6_heat) | (p6_alpha / p6_pcoa + plot_layout(heights = c(1, 1.3)))) +
  plot_layout(widths = c(1.3, 2.5)) + plot_annotation(tag_levels = "A")
save_fig(fig6, "Fig6_origin_x_species", width = 215, height = 250)
if (isTRUE(capabilities("cairo")))
  ggsave(file.path(paths$figures, "Fig6_origin_x_species.pdf"), fig6, width = 215, height = 250,
         units = "mm", device = cairo_pdf)

# =============================================================================
# PART B - wild host species (wild only)
# =============================================================================
B <- sd |>
  filter(Origin == "Wild", Species %in% sp_levels) |>
  mutate(Species = factor(Species, levels = sp_levels))
B_n  <- table(B$Species)
B_df <- as.data.frame(B)
dB   <- jac(B$Sample_ID)

B_kw <- map_dfr(alpha_metrics, function(m) {
  t <- kruskal.test(B[[m]] ~ B$Species)
  tibble(metric = m, H = as.numeric(t$statistic), df = as.numeric(t$parameter), p = t$p.value)
}) |> mutate(p_BH = p.adjust(p, "BH"))
B_pairs <- combn(sp_levels, 2, simplify = FALSE)
B_alpha <- map_dfr(B_pairs, ~ pair_alpha(B, "Species", .x[1], .x[2])) |>
  group_by(metric) |> mutate(p_BH = p.adjust(p, "BH")) |> ungroup()
pB    <- adonis2(dB ~ log_depth + Species, data = B_df, permutations = 999, by = "margin")
dispB <- permutest(betadisper(dB, B_df$Species), permutations = 999)
B_beta <- map_dfr(B_pairs, ~ pair_beta(B, "Species", .x[1], .x[2])) |>
  mutate(p_BH = p.adjust(p, "BH"), permdisp_p_BH = p.adjust(permdisp_p, "BH"))

# Site-adjusted: species effect after collection site (location-assigned only)
BL    <- B |> filter(!is.na(Location)) |> mutate(Location = droplevels(Location))
dBL     <- jac(BL$Sample_ID)
pB_site <- adonis2(dBL ~ log_depth + Location + Species, data = as.data.frame(BL),
                   permutations = 999, by = "margin")

# Kranji Marshes only: Ae. albopictus vs Mansonia indiana
KM <- B |> filter(Location == "Kranji Marshes", Species %in% c("Ae. albopictus", "Mansonia indiana")) |>
  mutate(Species = droplevels(Species))
km_ok <- all(table(KM$Species) >= 5) && nlevels(KM$Species) == 2
if (km_ok) {
  KM_beta  <- pair_beta(KM, "Species", "Ae. albopictus", "Mansonia indiana")
  KM_alpha <- pair_alpha(KM, "Species", "Ae. albopictus", "Mansonia indiana") |> mutate(p_BH = p.adjust(p, "BH"))
  KM_core  <- core_by(KM, "Species")
  save_table(KM_beta, "11B_KM_albopictus_vs_mansonia_beta", paths$stats)
  save_table(KM_alpha, "11B_KM_albopictus_vs_mansonia_alpha", paths$stats)
}

B_core <- core_by(B, "Species")
save_table(B_kw, "11B_alpha_kruskal", paths$stats)
save_table(B_alpha, "11B_pairwise_alpha", paths$stats)
save_table(bind_rows(tidy_adonis(pB, "depth + Species, marginal"),
                     tidy_adonis(pB_site, "depth + Location + Species, marginal (location-assigned)")),
           "11B_permanova_wild_species", paths$stats)
save_table(B_beta, "11B_pairwise_beta", paths$stats)
save_table(B_core |> filter(prevalence >= 50), "11B_core_prevalence", paths$stats)

# ---- Fig. 7 (stats wild only; lab = grey reference) -------------------------
lab_ids <- sd$Sample_ID[sd$Origin == "Lab"]
C7 <- bind_rows(B |> mutate(Group = as.character(Species)),
                sd |> filter(Origin == "Lab") |> mutate(Group = "Lab")) |>
  mutate(Group = factor(Group, levels = c(sp_levels, "Lab")))
B_lab <- setNames(as.expression(c(
  lapply(sp_levels, function(s) bquote(atop(italic(.(sp_short[[s]])), "(" * .(as.integer(B_n[[s]])) * ")"))),
  list(bquote(atop("Lab", "(" * .(length(lab_ids)) * ")"))))), c(sp_levels, "Lab"))
B_heat_tab <- genus_mean_table(bio[C7$Sample_ID, ], C7$Group,
                               facet = factor(if_else(C7$Group == "Lab", "Ref.", "Wild species"),
                                              levels = c("Wild species", "Ref.")))
save_table(B_heat_tab, "11B_all_genera_by_wild_species", paths$stats)
p7_heat <- genus_heatmap(B_heat_tab, x_labels = function(b) B_lab[b],
                         breaks = c(1e-3, 0.1, 10), legend_width = 30)

B_long <- B |> select(Sample_ID, Species, all_of(alpha_metrics)) |>
  pivot_longer(all_of(alpha_metrics), names_to = "metric", values_to = "value") |>
  filter(!is.na(value)) |> mutate(metric = factor(alpha_names[metric], levels = alpha_names))
B_strip <- setNames(paste0(alpha_names[B_kw$metric], "\nK-W ", fmt_p(B_kw$p_BH)), alpha_names[B_kw$metric])
p7_alpha <- ggplot(B_long, aes(x = Species, y = value)) +
  geom_jitter(aes(colour = Species), width = 0.2, height = 0, size = 0.6, alpha = 0.4) +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.55, linewidth = 0.4) +
  facet_wrap(~ metric, scales = "free_y", labeller = as_labeller(B_strip)) +
  scale_colour_manual(values = sp_colours, guide = "none") +
  scale_x_discrete(labels = sp_short) +
  labs(x = NULL, y = NULL) +
  theme(strip.text = element_text(size = 8.5, lineheight = 0.9), strip.clip = "off",
        axis.text.x = element_text(face = "italic", size = 8, angle = 45, hjust = 1))

d7  <- jac(C7$Sample_ID)
pc7 <- cmdscale(d7, k = 2, eig = TRUE)
ve7 <- 100 * pc7$eig[1:2] / sum(pc7$eig[pc7$eig > 0])
pc7_df <- C7 |> mutate(PCo1 = pc7$points[Sample_ID, 1], PCo2 = pc7$points[Sample_ID, 2])
leg7 <- as.expression(c(
  lapply(sp_levels, function(s) {
    if (s == "Culex spp.") bquote(italic("Culex") ~ "spp. (n =" ~ .(as.integer(B_n[[s]])) * ")")
    else bquote(italic(.(s)) ~ "(n =" ~ .(as.integer(B_n[[s]])) * ")")
  }),
  list(bquote("Lab, reference (n =" ~ .(length(lab_ids)) * ")"))))
brk7 <- c(sp_levels, "Lab")
p7_pcoa <- ggplot(mapping = aes(x = PCo1, y = PCo2, colour = Group, shape = Group)) +
  geom_point(data = filter(pc7_df, Group == "Lab"), size = 1, alpha = 0.35) +
  geom_point(data = filter(pc7_df, Group != "Lab"), size = 1.3, alpha = 0.7) +
  stat_ellipse(data = filter(pc7_df, Group != "Lab"), aes(shape = NULL),
               type = "t", level = 0.95, linewidth = 0.8) +
  scale_colour_manual(values = c(sp_colours, Lab = "#9b9a95"), breaks = brk7, labels = leg7, name = NULL) +
  scale_shape_manual(values = c("Ae. albopictus" = 16, "Mansonia indiana" = 17, "Culex spp." = 15, Lab = 1),
                     breaks = brk7, labels = leg7, name = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE), shape = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = sprintf("PCoA 1 (%.1f%%)", ve7[1]), y = sprintf("PCoA 2 (%.1f%%)", ve7[2]),
       subtitle = sprintf("Jaccard | Wild species only\nPERMANOVA (depth-adjusted): R² = %.3f, F = %.1f, %s | PERMDISP: %s",
                          pB["Species", "R2"], pB["Species", "F"], fmt_p(pB["Species", "Pr(>F)"]),
                          fmt_p(dispB$tab[1, "Pr(>F)"]))) +
  theme(plot.subtitle = element_text(size = 8), legend.position = "bottom",
        legend.text = element_text(size = 8)) +
  coord_fixed()

fig7 <- (free(p7_heat) | (p7_alpha / p7_pcoa + plot_layout(heights = c(1, 1.3)))) +
  plot_layout(widths = c(1.5, 2.3)) + plot_annotation(tag_levels = "A")
save_fig(fig7, "Fig7_wild_host_species", width = 225, height = 250)
if (isTRUE(capabilities("cairo")))
  ggsave(file.path(paths$figures, "Fig7_wild_host_species.pdf"), fig7, width = 225, height = 250,
         units = "mm", device = cairo_pdf)

# =============================================================================
# Summary
# =============================================================================
ad_line <- function(a, term, label) sprintf("  %-34s R2 = %.3f, F = %.2f, %s", label,
                                            a[term, "R2"], a[term, "F"], fmt_p(a[term, "Pr(>F)"]))
beta_lines <- function(t) sprintf("  %-20s vs %-20s (n %3d/%3d): R2 = %.3f, F = %.2f, BH %s | PERMDISP BH %s",
                                  t$group1, t$group2, t$n1, t$n2, t$R2, t$F, fmt_p(t$p_BH), fmt_p(t$permdisp_p_BH))
alpha_lines <- function(t) sprintf("  %-8s %-20s %.2f vs %-20s %.2f: W = %.1f, BH %s",
                                   t$metric, t$group1, t$median1, t$group2, t$median2, t$W, fmt_p(t$p_BH))
top_genera <- function(tab) tab |> filter(Genus != "Unclassified_at_genus") |>
  group_by(Group) |> slice_max(mean_rel, n = 6, with_ties = FALSE) |>
  summarise(top = paste(sprintf("%s %.1f%% (prev %.0f%%)", pretty_genus(Genus), mean_rel, prevalence), collapse = "; "),
            .groups = "drop") |>
  mutate(line = sprintf("  %-22s %s", as.character(Group), top)) |> pull(line)

summary_lines <- c(
  "11_host_species.R - results summary",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "", "SAMPLE COUNTS (wild species from COI)",
  with(count(sd, Origin, Species), sprintf("  %-6s %-40s %d", as.character(Origin), as.character(Species), n)),
  "",
  "=== PART A (3.5.2): Origin x Species, Ae. albopictus + Culex spp. ===",
  sprintf("  n: %s", paste(sprintf("%s %d", names(A_n), as.integer(A_n)), collapse = ", ")),
  "PERMANOVA (Jaccard, 999 perm.)",
  ad_line(pA_main, "Origin",    "Origin | depth + Species:"),
  ad_line(pA_main, "Species",   "Species | depth + Origin:"),
  ad_line(pA_main, "log_depth", "log10(depth) | Origin + Species:"),
  ad_line(pA_int,  "Origin:Species", "Origin x Species (after main):"),
  sprintf("  PERMDISP (4 groups): F = %.2f, %s; distance to centroid %s", dispA$tab[1, "F"], fmt_p(dispA$tab[1, "Pr(>F)"]),
          paste(sprintf("%s %.3f", names(tapply(betadisper(dA, A_df$Group)$distances, A_df$Group, mean)),
                        tapply(betadisper(dA, A_df$Group)$distances, A_df$Group, mean)), collapse = ", ")),
  "Pairwise PERMANOVA (depth-adjusted; BH over 4 pairs)", beta_lines(A_beta),
  "Pairwise alpha diversity (medians; Wilcoxon, BH over 4 pairs within metric)", alpha_lines(A_alpha),
  "Core genera (>= 90%)", core_lines(A_core),
  "Top genera (mean relative abundance)", top_genera(A_heat_tab),
  "",
  "=== PART B (3.5.3): wild host species (wild only) ===",
  sprintf("  n: %s", paste(sprintf("%s %d", names(B_n), as.integer(B_n)), collapse = ", ")),
  "Alpha diversity (Kruskal-Wallis, BH over 3 metrics)",
  sprintf("  %-8s H = %.2f, df = %d, BH %s", B_kw$metric, B_kw$H, as.integer(B_kw$df), fmt_p(B_kw$p_BH)),
  alpha_lines(B_alpha),
  "PERMANOVA",
  ad_line(pB, "Species", "Species | depth:"),
  sprintf("  PERMDISP: F = %.2f, %s", dispB$tab[1, "F"], fmt_p(dispB$tab[1, "Pr(>F)"])),
  ad_line(pB_site, "Species",  "Species | depth + site:"),
  ad_line(pB_site, "Location", "Site | depth + species:"),
  sprintf("  (site-adjusted model uses %d location-assigned wild mosquitoes)", nrow(BL)),
  "Pairwise PERMANOVA (depth-adjusted; BH over 3 pairs)", beta_lines(B_beta),
  "Core genera (>= 90%)", core_lines(B_core),
  "Top genera (mean relative abundance)", top_genera(B_heat_tab |> filter(Group != "Lab")),
  "",
  "=== KRANJI MARSHES ONLY: Ae. albopictus vs Mansonia indiana ===",
  if (km_ok) c(sprintf("  n: %s", paste(sprintf("%s %d", levels(KM$Species), as.integer(table(KM$Species))), collapse = ", ")),
               beta_lines(mutate(KM_beta, p_BH = p, permdisp_p_BH = permdisp_p)),
               alpha_lines(KM_alpha), core_lines(KM_core))
  else "  skipped: fewer than 5 mosquitoes of one species at Kranji Marshes",
  "",
  "Wild species x collection site:",
  capture.output(print(sd |> filter(Origin == "Wild") |> count(Species, Location) |>
                         pivot_wider(names_from = Location, values_from = n, values_fill = 0), n = Inf))
)
writeLines(summary_lines, file.path(paths$stats, "11_results_summary.txt"))
cat(summary_lines, sep = "\n")