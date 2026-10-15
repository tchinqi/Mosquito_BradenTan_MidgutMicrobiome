# =============================================================================
# 11a_build_host_species.R - tag COI host-species IDs onto the wild 16S samples
#
# Wild 16S plates: AG_P01, AG_P02, ZX_P01, ZX_P02 (380 wild midguts).
# COI identifications (one per specimen):
#   Agabus plate 1 : CO1_P1_Gbif-Results.xlsx (or .csv)  - GBIF/BOLD BLAST
#   Agabus plate 2 : CO1_P2_Gbif-Results.xlsx (or .csv)  - OPTIONAL, used if present
#                    otherwise TYTMZ2_COI_morpho_id.xlsx  - MORPHOLOGICAL IDs (Oct 2026)
#   Zi Xin         : ZiXin_Taxa_classified_COI.xlsx, Sheet1 (Sheet2 is a subset)
# 16S sample <-> specimen: demux_metadata_audit.csv (from Atlas9 01_demux_metadata/)
#
# Match key = specimen name in upper case with all non-alphanumerics removed
# (COI files write e.g. "COI07km3-13"; 16S writes "07KM3-13").
#
# Put the input files here:
#   data/raw/species/CO1_P1_Gbif-Results.xlsx
#   data/raw/species/ZiXin_Taxa_classified_COI.xlsx
#   data/raw/metadata/demux_metadata_audit.csv
# Output:
#   data/raw/species/wild_host_species_COI.tsv   (one row per wild 16S sample)
# Plate 2: Aggy's morpho sheet data/raw/species/TYTMZ2_COI_morpho_id.xlsx is used
# (columns Specimen_code_number, Specimen_morpho_ID). If a GBIF/BLAST file for
# plate 2 (CO1_P2_Gbif-Results) is added later, it takes precedence.
# =============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(readxl))

sp_dir  <- paths$species
f_audit <- c(here("data", "raw", "metadata", "demux_metadata_audit.csv"),
             file.path(paths$species, "demux_metadata_audit.csv"))
f_audit <- f_audit[file.exists(f_audit)][1]
if (is.na(f_audit)) f_audit <- "missing"
find_file <- function(base) {
  f <- file.path(sp_dir, paste0(base, c(".xlsx", ".csv")))
  f[file.exists(f)][1]
}
read_any <- function(f) if (grepl("\\.csv$", f)) read_csv(f, show_col_types = FALSE) else read_excel(f)
norm <- function(x) gsub("[^A-Z0-9]", "", toupper(as.character(x)))

f_ag1 <- find_file("CO1_P1_Gbif-Results")
f_ag2 <- find_file("CO1_P2_Gbif-Results")          # GBIF/BLAST version (not available yet)
f_ag2m <- file.path(sp_dir, "TYTMZ2_COI_morpho_id.xlsx")   # morphological IDs
if (!file.exists(f_ag2m)) f_ag2m <- NA_character_
f_zx  <- file.path(sp_dir, "ZiXin_Taxa_classified_COI.xlsx")
stopifnot("demux_metadata_audit.csv not found in data/raw/metadata/" = file.exists(f_audit),
          "CO1_P1_Gbif-Results not found in data/raw/species/"       = !is.na(f_ag1),
          "ZiXin_Taxa_classified_COI.xlsx not found in data/raw/species/" = file.exists(f_zx))

# ---- 1. The 380 wild 16S samples --------------------------------------------
wild <- read_csv(f_audit, show_col_types = FALSE) |>
  filter(Plate_ID %in% c("AG_P01_16S", "AG_P02_16S", "ZX_P01_16S", "ZX_P02_16S"),
         !str_detect(Specimen, regex("neg", ignore_case = TRUE))) |>
  transmute(Sample_ID = paste0("sample_", SampleID), Specimen,
            Plate = str_remove(Plate_ID, "_16S"), Well, key = norm(Specimen))
stopifnot(nrow(wild) == 380, !anyDuplicated(wild$key))

# ---- 2. COI identifications -------------------------------------------------
last_n <- function(x, n) vapply(str_split(x, "_"), function(p) p[length(p) - n], character(1))
read_gbif <- function(f, label) {
  d <- read_any(f)
  tibble(key          = norm(sub("^COI", "", sub("_all.*$", "", d$occurrenceId))),
         COI_source   = label,
         COI_call_raw = last_n(d$classification, 0),
         COI_rank     = "SPECIES",
         COI_support  = sprintf("%.2f%% identity", as.numeric(d$identity)),
         Genus_raw    = last_n(d$classification, 1))
}
# Aggy plate 2 morpho sheet: abbreviated names ("Ae. albopictus", "Ma. uniformis")
# and long specimen codes that the 16S plate writes in short form:
#   240926MZFF_97 -> 0926FF_97   (Mandai: drop year and "MZ")
#   250122MT1_1   -> 0122MT1_1   (MT: drop year)
#   230907KM1_9   -> 07KM1_09    (Kranji: day only, number zero-padded)
to_16s_name <- function(code) {
  code <- str_trim(as.character(code))
  m <- str_match(code, "^(\\d{2})(\\d{2})(\\d{2})(?:MZ)?([A-Z]+\\d*)_(\\d+)$")
  out <- if_else(str_starts(m[, 5], "KM"),
                 sprintf("%s%s_%02d", m[, 4], m[, 5], as.integer(m[, 6])),
                 sprintf("%s%s%s_%s", m[, 3], m[, 4], m[, 5], m[, 6]))
  if_else(is.na(m[, 1]), code, out)
}
read_morpho <- function(f, label) {
  d  <- read_excel(f)
  id <- str_squish(as.character(d$Specimen_morpho_ID))
  genus_map <- c("Ae." = "Aedes", "Ma." = "Mansonia", "Cx." = "Culex",
                 "An." = "Anopheles", "Ar." = "Armigeres", "Ur." = "Uranotaenia")
  ab    <- str_extract(id, "^[A-Z][a-z]+\\.")
  genus <- coalesce(unname(genus_map[ab]), word(id, 1))
  full  <- if_else(!is.na(ab), paste(genus, word(id, 2)), id)   # "Ae. albopictus" -> "Aedes albopictus"
  tibble(key          = norm(to_16s_name(d$Specimen_code_number)),
         COI_source   = label,
         COI_call_raw = full,
         COI_rank     = "MORPHOLOGY",
         COI_support  = "",
         Genus_raw    = genus)
}

zx <- read_excel(f_zx, sheet = "Sheet1")
zx_tab <- tibble(key          = norm(zx$SpecimenID),
                 COI_source   = "Zi Xin (tree placement)",
                 COI_call_raw = coalesce(as.character(zx$Species), as.character(zx$Genus)),
                 COI_rank     = as.character(zx$`Tax Rank`),
                 COI_support  = ifelse(is.na(zx$Confidence), "",
                                       sprintf("%.0f%% confidence", as.numeric(zx$Confidence))),
                 Genus_raw    = as.character(zx$Genus))

coi <- bind_rows(read_gbif(f_ag1, "Agabus P1 (GBIF/BOLD BLAST)"),
                 if (!is.na(f_ag2)) read_gbif(f_ag2, "Agabus P2 (GBIF/BOLD BLAST)")
                 else if (!is.na(f_ag2m)) read_morpho(f_ag2m, "Agabus P2 (morphological ID)"),
                 zx_tab)
stopifnot("a specimen appears in more than one COI file" = !anyDuplicated(coi$key))

# ---- 3. Join + harmonise names ----------------------------------------------
host <- wild |>
  left_join(coi, by = "key") |>
  mutate(Host_group = case_when(
    Plate == "AG_P02" & is.na(COI_call_raw) & is.na(f_ag2) & is.na(f_ag2m) ~ "Pending (AG_P02 COI not yet available)",
    is.na(COI_call_raw)                                      ~ "Unidentified",
    COI_call_raw == "Aedes albopictus"                       ~ "Ae. albopictus",
    COI_call_raw == "Aedes aegypti"                          ~ "Ae. aegypti",
    Genus_raw == "Culex"                                     ~ "Culex spp.",
    Genus_raw == "Mansonia" & str_detect(COI_call_raw, "indiana") ~ "Mansonia indiana",
    Genus_raw == "Mansonia"                                  ~ "Mansonia spp.",
    Genus_raw == "Anopheles"                                 ~ "Anopheles spp.",
    TRUE                                                     ~ paste("Other:", COI_call_raw))) |>
  select(Sample_ID, Specimen, Plate, Well, Host_group, COI_call_raw, COI_rank, COI_support, COI_source) |>
  arrange(Plate, Sample_ID)

out <- file.path(sp_dir, "wild_host_species_COI.tsv")
write_tsv(host, out)

# ---- 4. Report --------------------------------------------------------------
cat("\nHost species of the 380 wild 16S samples:\n")
print(count(host, Host_group, sort = TRUE), n = Inf)
cat("\nBy plate:\n")
print(count(host, Plate, Host_group), n = Inf)
cat("\nCulex detail:\n")
print(count(filter(host, Host_group == "Culex spp."), COI_call_raw, sort = TRUE))
unmatched <- setdiff(coi$key, wild$key)
cat("\nCOI records with no wild 16S sample (other projects / not sequenced):", length(unmatched), "\n")
if (!is.na(f_ag2m) || !is.na(f_ag2)) {
  p2 <- filter(host, Plate == "AG_P02")
  cat("\nAG_P02: ", sum(!is.na(p2$COI_call_raw)), " of ", nrow(p2),
      " wild 16S samples matched a plate-2 ID\n", sep = "")
  miss <- filter(p2, is.na(COI_call_raw))
  if (nrow(miss)) { cat("AG_P02 samples with NO plate-2 ID (check specimen names):\n"); print(select(miss, Sample_ID, Specimen), n = Inf) }
}
if (is.na(f_ag2) && is.na(f_ag2m)) {
  expected <- c("Ae. albopictus" = 219, "Mansonia indiana" = 36, "Culex spp." = 23,
                "Pending (AG_P02 COI not yet available)" = 95, "Unidentified" = 4)
  got <- table(host$Host_group)[names(expected)]
  if (!isTRUE(all(got == expected))) warning("Counts differ from the 2026-10-08 build - check the input files")
}
message("Saved ", out)