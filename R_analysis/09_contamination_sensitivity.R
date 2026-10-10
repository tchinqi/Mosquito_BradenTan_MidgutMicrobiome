# =============================================================================
# 09_contamination_sensitivity.R - are the results robust to cross-contamination?
#
# The negative controls (H12 of every plate) contained the same dominant
# genera as mosquito samples, including Wolbachia (an obligate intracellular
# symbiont that cannot come from reagents) -> sample-to-sample
# cross-contamination (well-to-well transfer and/or tag misassignment).
# Contamination of this kind mostly adds LOW-abundance reads, which inflate
# "presence". This script repeats the key analyses after requiring a minimum
# within-sample relative abundance before a taxon counts as present.
#
#   1. Core genera (as script 07) with presence = >= 1 read, >= 0.5%, >= 1%,
#      >= 2%, >= 5% of a mosquito's classified reads
#   2. Wild vs Lab and among-site tests (as scripts 01/02) after removing
#      OTU detections below 1% within each mosquito
#   3. Composition of each negative control vs the mean of the samples on the
#      same plate (descriptive: do blanks look like their plate's samples?)
#
# Outputs
#   output/figures/FigS_core_contamination_sensitivity.pdf/.png
#   output/stats/09_*.tsv
#   output/stats/09_results_summary.txt
# =============================================================================

source(here::here("R", "00_setup.R"))
check_inputs()

ra_thresholds <- c(0, 0.005, 0.01, 0.02, 0.05)   # 0 = >= 1 read (main analysis)
otu_filter    <- 0.01                            # for section 2
is_core <- function(prev, t) prev >= t - 1e-9
fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, digits = 2, format = "g")))
plate_of <- function(id) str_extract(id, "(AG|ZX|BT)_P0[12]")

# ---- 1. Load ----------------------------------------------------------------
meta <- read_metadata()
bio  <- read_otu_matrix("biological")
ctrl <- read_otu_matrix("controls")
gen     <- genus_counts(bio)
gen_rel <- gen / rowSums(gen)

sd <- tibble(Sample_ID = rownames(bio)) |>
  left_join(select(meta, Sample_ID, Origin, Location, Lab_Taxon), by = "Sample_ID") |>
  mutate(Origin = droplevels(Origin), Plate = plate_of(Sample_ID))
stopifnot(identical(sd$Sample_ID, rownames(gen)), !anyNA(sd$Plate))

# ---- 2. Groups (identical to script 07) -------------------------------------
group_levels <- c("Wild", "Lab", "Mandai Zoo", "Kranji Marshes", "Bird Paradise",
                  "River Safari", "Ae. albopictus", "Culex")
membership <- bind_rows(
  sd |> filter(Origin == "Wild")  |> transmute(Sample_ID, Group = "Wild"),
  sd |> filter(Origin == "Lab")   |> transmute(Sample_ID, Group = "Lab"),
  sd |> filter(!is.na(Location))  |> transmute(Sample_ID, Group = as.character(Location)),
  sd |> filter(!is.na(Lab_Taxon)) |> transmute(Sample_ID, Group = Lab_Taxon)
) |> mutate(Group = factor(Group, levels = group_levels))
stopifnot(identical(as.integer(table(membership$Group)), c(380L, 190L, 194L, 88L, 39L, 24L, 165L, 25L)))
n_groups <- length(group_levels)
genera   <- setdiff(colnames(gen), "Unclassified_at_genus")

# ---- 3. Core genera under stricter presence definitions ---------------------
prev_tab <- map_dfr(ra_thresholds, function(t) {
  pres <- if (t == 0) gen > 0 else gen_rel >= t
  membership |>
    group_by(Group) |>
    group_modify(function(d, k) {
      tibble(Genus = genera,
             prevalence = unname(100 * colMeans(pres[d$Sample_ID, genera, drop = FALSE])))
    }) |>
    ungroup() |>
    mutate(min_ra = t)
})
save_table(prev_tab, "09_genus_prevalence_by_presence_threshold", paths$stats)

core_summary <- prev_tab |>
  group_by(min_ra, Group) |>
  summarise(core90 = sum(is_core(prevalence, 90)), core70 = sum(is_core(prevalence, 70)),
            core90_genera = paste(pretty_genus(Genus[is_core(prevalence, 90)]), collapse = ", "),
            .groups = "drop")
shared <- prev_tab |>
  group_by(min_ra, Genus) |>
  summarise(core_in = sum(is_core(prevalence, 90)), min_prev = min(prevalence), .groups = "drop") |>
  group_by(min_ra) |>
  summarise(shared_core_90 = paste(pretty_genus(Genus[core_in == n_groups]), collapse = ", "),
            n_shared = sum(core_in == n_groups), .groups = "drop")
save_table(core_summary, "09_core_by_presence_threshold", paths$stats)
save_table(shared, "09_shared_core_by_presence_threshold", paths$stats)

# Figure: prevalence of the three dominant genera vs presence threshold
focus <- c("Stenotrophomonas_A_615274", "Elizabethkingia", "Wolbachia")
grp_cols <- c("Wild" = "#1baf7a", "Lab" = "#4a3aa7", "Mandai Zoo" = "#2a78d6",
              "Kranji Marshes" = "#eb6834", "Bird Paradise" = "#008300",
              "River Safari" = "#e87ba4", "Ae. albopictus" = "#8c5a2b", "Culex" = "#9b9a95")
thr_lab <- function(t) ifelse(t == 0, "≥1 read", paste0("≥", 100 * t, "%"))
p_sens <- prev_tab |>
  filter(Genus %in% focus) |>
  mutate(Genus = factor(pretty_genus(Genus), levels = pretty_genus(focus)),
         thr = factor(thr_lab(min_ra), levels = thr_lab(ra_thresholds))) |>
  ggplot(aes(x = thr, y = prevalence, colour = Group, group = Group)) +
  geom_hline(yintercept = 90, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.4) +
  facet_wrap(~ Genus, nrow = 1) +
  scale_colour_manual(values = grp_cols, name = NULL, breaks = names(grp_cols),
                      labels = c(as.list(names(grp_cols)[1:6]),
                                 expression(italic("Ae. albopictus")),
                                 expression(italic("Culex") ~ "spp."))) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  labs(x = "Minimum share of a mosquito's reads for the genus to count as present",
       y = "Prevalence (% of mosquitoes)") +
  theme(strip.text = element_text(face = "bold.italic"), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
save_fig(p_sens, "FigS_core_contamination_sensitivity", width = 200, height = 105)
if (isTRUE(capabilities("cairo"))) {
  ggsave(file.path(paths$figures, "FigS_core_contamination_sensitivity.pdf"), p_sens,
         width = 200, height = 105, units = "mm", device = cairo_pdf)
}

# ---- 4. Wild vs Lab and sites after removing OTU detections < 1% ------------
filt <- bio
filt[(bio / rowSums(bio)) < otu_filter] <- 0
keep_s <- rowSums(filt) > 0
stopifnot(all(keep_s))
sd_f <- sd |>
  mutate(depth     = as.numeric(rowSums(filt)),
         log_depth = log10(as.numeric(rowSums(bio))),     # original depth as covariate
         observed  = as.numeric(rowSums(filt > 0)),
         shannon   = as.numeric(diversity(filt, index = "shannon")),
         pielou    = if_else(observed > 1, shannon / log(observed), NA_real_))

alpha_f <- map_dfr(c("observed", "shannon", "pielou"), function(m) {
  x <- sd_f[[m]][sd_f$Origin == "Wild"]; y <- sd_f[[m]][sd_f$Origin == "Lab"]
  t <- wilcox.test(x[!is.na(x)], y[!is.na(y)], exact = FALSE)
  tibble(metric = m, median_wild = median(x, na.rm = TRUE), median_lab = median(y, na.rm = TRUE),
         W = as.numeric(t$statistic), p = t$p.value)
}) |> mutate(p_BH = p.adjust(p, method = "BH"))
save_table(alpha_f, "09_alpha_wild_vs_lab_otu_filtered", paths$stats)

sd_f_df <- as.data.frame(sd_f)
d_f   <- vegdist(filt, method = "jaccard", binary = TRUE)
perm_f <- adonis2(d_f ~ log_depth + Origin, data = sd_f_df, permutations = 999, by = "margin")
disp_f <- permutest(betadisper(d_f, sd_f_df$Origin), permutations = 999)

loc <- !is.na(sd_f$Location)
sd_l_df <- as.data.frame(sd_f[loc, ]) |> mutate(Location = droplevels(Location))
d_l    <- vegdist(filt[loc, ], method = "jaccard", binary = TRUE)
perm_l <- adonis2(d_l ~ log_depth + Location, data = sd_l_df, permutations = 999, by = "margin")
disp_l <- permutest(betadisper(d_l, sd_l_df$Location), permutations = 999)

beta_f <- tibble(
  comparison = c("Wild vs Lab", "Wild sites"),
  R2 = c(perm_f["Origin", "R2"], perm_l["Location", "R2"]),
  F  = c(perm_f["Origin", "F"],  perm_l["Location", "F"]),
  p  = c(perm_f["Origin", "Pr(>F)"], perm_l["Location", "Pr(>F)"]),
  permdisp_p = c(disp_f$tab[1, "Pr(>F)"], disp_l$tab[1, "Pr(>F)"]))
save_table(beta_f, "09_beta_otu_filtered", paths$stats)

# ---- 5. Do blanks resemble the samples on their own plate? ------------------
ctrl_gen <- genus_counts(ctrl)
ctrl_rel <- ctrl_gen / rowSums(ctrl_gen)
blank_vs_plate <- map_dfr(rownames(ctrl_rel), function(cid) {
  pl <- plate_of(cid)
  plate_mean <- colMeans(gen_rel[sd$Plate == pl, , drop = FALSE])
  b <- ctrl_rel[cid, names(plate_mean)]
  top <- names(sort(b, decreasing = TRUE))[1:3]
  tibble(Control = cid, Plate = pl, reads = sum(ctrl_gen[cid, ]),
         bray_curtis_to_plate_mean = as.numeric(vegdist(rbind(b, plate_mean), method = "bray")),
         top3_genera = paste(sprintf("%s %.0f%%", pretty_genus(top), 100 * b[top]), collapse = "; "))
})
save_table(blank_vs_plate, "09_blanks_vs_plate_samples", paths$stats)

# ---- 6. Summary -------------------------------------------------------------
# Main-analysis values for comparison, read from the outputs of scripts 01 and 02
main_line <- tryCatch({
  rd <- function(f) read_tsv(file.path(paths$stats, f), show_col_types = FALSE)
  o  <- rd("01_permanova_wild_vs_lab.tsv") |> filter(model == "log10(depth) + Origin, marginal", term == "Origin")
  l  <- rd("02_permanova_locations.tsv")   |> filter(model == "log10(depth) + Location, marginal", term == "Location")
  sprintf("  (main analysis: Wild vs Lab R2 = %.3f, PERMDISP %s; sites R2 = %.3f, PERMDISP %s)",
          o$R2, fmt_p(rd("01_permdisp_wild_vs_lab.tsv")$p[1]),
          l$R2, fmt_p(rd("02_permdisp_locations.tsv")$p[1]))
}, error = function(e) "  (main analysis: run scripts 01 and 02 first to print their values here)")

foc_lines <- prev_tab |>
  filter(Genus %in% focus) |>
  group_by(min_ra, Genus) |>
  summarise(range = sprintf("%.0f-%.0f%%", min(prevalence), max(prevalence)),
            lowest = as.character(Group[which.min(prevalence)]), .groups = "drop") |>
  mutate(line = sprintf("  %-10s %-20s prevalence across groups %s (lowest: %s)",
                        thr_lab(min_ra), pretty_genus(Genus), range, lowest)) |>
  pull(line)

summary_lines <- c(
  "09_contamination_sensitivity.R - results summary",
  paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  "1. SHARED CORE (>= 90% in all 8 groups) BY PRESENCE DEFINITION",
  sprintf("  %-10s %d genera: %s", thr_lab(shared$min_ra), shared$n_shared, shared$shared_core_90),
  "",
  "   Prevalence of the three dominant genera (Stenotrophomonas A, Elizabethkingia, Wolbachia):",
  foc_lines,
  "",
  "   Number of core genera (90%) per group:",
  core_summary |>
    select(min_ra, Group, core90) |>
    pivot_wider(names_from = Group, values_from = core90) |>
    mutate(line = paste0(sprintf("  %-10s", thr_lab(min_ra)),
                         apply(across(-min_ra), 1, function(v) paste(sprintf("%4d", v), collapse = "")))) |>
    pull(line),
  paste0("  ", strrep(" ", 10), paste(sprintf("%4s", c("Wild", "Lab", "MZ", "KM", "BP", "RS", "Alb", "Clx")), collapse = "")),
  "",
  sprintf("2. COMPARISONS AFTER REMOVING OTU DETECTIONS < %.0f%% WITHIN EACH MOSQUITO", 100 * otu_filter),
  sprintf("  %-9s Wild median %.2f vs Lab %.2f; W = %.0f, BH %s",
          alpha_f$metric, alpha_f$median_wild, alpha_f$median_lab, alpha_f$W, fmt_p(alpha_f$p_BH)),
  sprintf("  %-12s PERMANOVA (depth-adjusted) R2 = %.3f, F = %.2f, %s | PERMDISP %s",
          beta_f$comparison, beta_f$R2, beta_f$F, fmt_p(beta_f$p), fmt_p(beta_f$permdisp_p)),
  main_line,
  "",
  "3. NEGATIVE CONTROLS VS THE MEAN PROFILE OF SAMPLES ON THE SAME PLATE",
  sprintf("  %-36s %-7s %5d reads | Bray-Curtis to plate mean %.2f | top: %s",
          blank_vs_plate$Control, blank_vs_plate$Plate, as.integer(blank_vs_plate$reads),
          blank_vs_plate$bray_curtis_to_plate_mean, blank_vs_plate$top3_genera),
  "  (Bray-Curtis 0 = identical profile, 1 = no shared genera)"
)
writeLines(summary_lines, file.path(paths$stats, "09_results_summary.txt"))
cat(summary_lines, sep = "\n")
