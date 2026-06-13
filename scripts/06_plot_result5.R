#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."), mustWork = TRUE)
input <- file.path(root, "results", "sample_abundance_matrix.tsv")
figdir <- file.path(root, "results", "figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input)) {
  stop("Missing results/sample_abundance_matrix.tsv. Expected columns: sample_id, latitude, longitude, soil_pH, Candidatus_Solibacter_abundance, pmtA_abundance, pmtA_per_Solibacter_ratio")
}

df <- read_tsv(input, show_col_types = FALSE)
required <- c("sample_id", "latitude", "longitude", "soil_pH",
              "Candidatus_Solibacter_abundance", "pmtA_abundance",
              "pmtA_per_Solibacter_ratio")
missing <- setdiff(required, names(df))
if (length(missing) > 0) {
  stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
}

cor_out <- tibble(
  comparison = c("Solibacter_vs_pH", "pmtA_vs_pH", "pmtA_per_Solibacter_vs_pH", "Solibacter_vs_pmtA"),
  spearman_rho = c(
    cor(df$Candidatus_Solibacter_abundance, df$soil_pH, method = "spearman", use = "complete.obs"),
    cor(df$pmtA_abundance, df$soil_pH, method = "spearman", use = "complete.obs"),
    cor(df$pmtA_per_Solibacter_ratio, df$soil_pH, method = "spearman", use = "complete.obs"),
    cor(df$Candidatus_Solibacter_abundance, df$pmtA_abundance, method = "spearman", use = "complete.obs")
  )
)
write_tsv(cor_out, file.path(root, "results", "result5_correlations.tsv"))

p_map <- ggplot(df, aes(x = longitude, y = latitude)) +
  geom_point(aes(size = Candidatus_Solibacter_abundance, color = pmtA_abundance), alpha = 0.75) +
  scale_color_viridis_c(option = "C", trans = "log10") +
  scale_size_continuous(range = c(1, 7)) +
  coord_quickmap(xlim = c(-180, 180), ylim = c(-60, 85)) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", color = "pmtA abundance", size = "Candidatus Solibacter abundance")
ggsave(file.path(figdir, "Fig5d_global_map.pdf"), p_map, width = 8, height = 4.8)

plot_cor <- function(x, y, xlab, ylab, file) {
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(alpha = 0.75) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    theme_classic() +
    labs(x = xlab, y = ylab)
  ggsave(file.path(figdir, file), p, width = 4.8, height = 4.2)
}

plot_cor("soil_pH", "Candidatus_Solibacter_abundance", "Soil pH", "Candidatus Solibacter abundance", "Fig5e_solibacter_vs_pH.pdf")
plot_cor("soil_pH", "pmtA_abundance", "Soil pH", "pmtA abundance", "Fig5f_pmta_vs_pH.pdf")
plot_cor("soil_pH", "pmtA_per_Solibacter_ratio", "Soil pH", "pmtA / Candidatus Solibacter abundance", "Fig5f_pmta_ratio_vs_pH.pdf")
plot_cor("Candidatus_Solibacter_abundance", "pmtA_abundance", "Candidatus Solibacter abundance", "pmtA abundance", "Fig5g_solibacter_vs_pmta.pdf")
