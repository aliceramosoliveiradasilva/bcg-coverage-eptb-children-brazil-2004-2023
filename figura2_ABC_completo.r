# ============================================================
# FIGURE 2 - Brazil EPTB
# (A) Descriptive map: cumulative notification rate quartiles
# (B) LISA map: High-High / Low-Low / Not significant
# (C) SaTScan space-time clusters (permutation), default = HIGH only
#
# Outputs:
#   figura2_ABC.png
#   figura2_ABC.tiff
#   figura2_ABC.pdf
# ============================================================

library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(spdep)
library(patchwork)

# -------------------------
# 1) INPUTS - EDIT HERE
# -------------------------

# municipality-year panel with cases and population
arq_painel <- "/home/ramos/Documentos/R/TUBERCULOSE/SATSCAN/painel_municipio_ano_ibge7_2004_2023_ok.csv"

# municipal shapefile
arq_brasil <- "/home/ramos/Documentos/R/TUBERCULOSE/IBGE/BR_Municipios_2025/BR_Municipios_2025.shp"

# SaTScan directory for panel C
# choose one:
# "/home/ramos/Documentos/R/TUBERCULOSE/SATSCAN/resultado da permutacao 10"
# "/home/ramos/Documentos/R/TUBERCULOSE/SATSCAN/resultado permutacao 20"
# "/home/ramos/Documentos/R/TUBERCULOSE/SATSCAN/resultado permutacao 50"
dir_cluster <- "/home/ramos/Documentos/R/TUBERCULOSE/SATSCAN/resultado permutacao 50"

# period for panels A and B
ano_ini <- 2004
ano_fim <- 2023

# show low clusters in panel C?
mostrar_low <- FALSE

# column names in municipality-year panel
col_mun   <- "mun_code"
col_ano   <- "year"
col_cases <- "cases"
col_pop   <- "population"

# rate multiplier
mult_taxa <- 100000

# -------------------------
# 2) HELPER FUNCTIONS
# -------------------------

achar_arquivo <- function(dir, padrao) {
  arqs <- list.files(dir, pattern = padrao, full.names = TRUE, ignore.case = TRUE)
  if (length(arqs) == 0) stop(paste("Arquivo não encontrado:", dir, "padrão:", padrao))
  arqs[1]
}

ler_painel <- function(arquivo) {
  ext <- tolower(tools::file_ext(arquivo))
  if (ext == "csv") {
    readr::read_csv(arquivo, show_col_types = FALSE)
  } else if (ext == "rds") {
    readRDS(arquivo)
  } else if (ext == "fst") {
    if (!requireNamespace("fst", quietly = TRUE)) {
      stop("Pacote fst não instalado. Instale fst ou salve o painel em CSV/RDS.")
    }
    fst::read_fst(arquivo, as.data.table = FALSE)
  } else {
    stop("Formato não suportado. Use CSV, RDS ou FST.")
  }
}

ler_col_txt_permutacao <- function(arquivo) {
  suppressMessages({
    readr::read_table2(
      file = arquivo,
      col_names = c(
        "cluster","loc_id","lat","lon","radius_km","span_km",
        "start_date","end_date","n_locations","test_stat","p_value",
        "observed","expected","ode","gini"
      ),
      locale = readr::locale(decimal_mark = ",", grouping_mark = ".")
    )
  }) %>%
    mutate(
      cluster_type = if_else(ode > 1, "High", "Low"),
      period = paste(start_date, "to", end_date)
    )
}

ler_shp_col <- function(dir_cluster) {
  shp <- achar_arquivo(dir_cluster, "\\.col\\.shp$")
  x <- st_read(shp, quiet = TRUE)

  nm <- names(x)
  id_col <- nm[str_detect(tolower(nm), "^cluster$|clust")]
  if (length(id_col) == 0) stop("Não encontrei a coluna de cluster no shapefile do SaTScan.")

  x %>%
    mutate(cluster = as.integer(.data[[id_col[1]]])) %>%
    select(cluster, geometry)
}

# -------------------------
# 3) READ DATA
# -------------------------

painel <- ler_painel(arq_painel)
brasil <- st_read(arq_brasil, quiet = TRUE)

# municipality code in shapefile
cd_col <- names(brasil)[str_detect(tolower(names(brasil)), "^cd_mun$|cod.*mun|geocod")]
if (length(cd_col) == 0) stop("Não encontrei coluna do código municipal no shapefile.")
cd_col <- cd_col[1]

brasil <- brasil %>%
  mutate(mun_code = str_pad(as.character(.data[[cd_col]]), width = 7, side = "left", pad = "0"))

# standardize panel
painel2 <- painel %>%
  mutate(
    mun_code   = str_pad(as.character(.data[[col_mun]]), width = 7, side = "left", pad = "0"),
    year       = as.integer(.data[[col_ano]]),
    cases      = as.numeric(.data[[col_cases]]),
    population = as.numeric(.data[[col_pop]])
  ) %>%
  filter(year >= ano_ini, year <= ano_fim)

# -------------------------
# 4) MUNICIPAL AGGREGATION FOR A AND B
# -------------------------

mun_agg <- painel2 %>%
  group_by(mun_code) %>%
  summarise(
    cases_sum = sum(cases, na.rm = TRUE),
    pop_sum   = sum(population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    taxa_acum = if_else(pop_sum > 0, mult_taxa * cases_sum / pop_sum, NA_real_)
  )

map_ab <- brasil %>%
  left_join(mun_agg, by = "mun_code")

# -------------------------
# 5) PANEL A - DESCRIPTIVE QUARTILES
# -------------------------

q_vals <- quantile(map_ab$taxa_acum[map_ab$cases_sum > 0], probs = c(0.25, 0.50, 0.75), na.rm = TRUE)

map_a <- map_ab %>%
  mutate(
    classe_A = case_when(
      is.na(taxa_acum) | is.na(cases_sum) ~ "No reported EPTB cases",
      cases_sum == 0 ~ "No reported EPTB cases",
      taxa_acum <= q_vals[1] ~ "1st quartile",
      taxa_acum <= q_vals[2] ~ "2nd quartile",
      taxa_acum <= q_vals[3] ~ "3rd quartile",
      taxa_acum >  q_vals[3] ~ "4th quartile",
      TRUE ~ "No reported EPTB cases"
    )
  )

map_a$classe_A <- factor(
  map_a$classe_A,
  levels = c("1st quartile", "2nd quartile", "3rd quartile", "4th quartile", "No reported EPTB cases")
)

pA <- ggplot() +
  geom_sf(data = map_a, aes(fill = classe_A), color = NA) +
  scale_fill_manual(
    values = c(
      "1st quartile" = "#f6d37b",
      "2nd quartile" = "#fbb040",
      "3rd quartile" = "#ff6b00",
      "4th quartile" = "#9e1b0c",
      "No reported EPTB cases" = "#6b6b6b"
    ),
    drop = FALSE
  ) +
  labs(fill = NULL, title = "(A)") +
  coord_sf() +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
    axis.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0)
  )

# -------------------------
# 6) PANEL B - LISA
# -------------------------

nb <- spdep::poly2nb(map_ab, queen = TRUE)
lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

x <- map_ab$taxa_acum
x[is.na(x)] <- 0

lisa <- spdep::localmoran(x, lw, zero.policy = TRUE)
p_col <- grep("^Pr", colnames(lisa), value = TRUE)[1]
lag_x <- spdep::lag.listw(lw, x, zero.policy = TRUE)

x_center <- x - mean(x, na.rm = TRUE)
lag_center <- lag_x - mean(lag_x, na.rm = TRUE)

map_b <- map_ab %>%
  mutate(
    lisa_p = lisa[, p_col],
    lisa_class = case_when(
      lisa_p <= 0.05 & x_center > 0 & lag_center > 0 ~ "High-High",
      lisa_p <= 0.05 & x_center < 0 & lag_center < 0 ~ "Low-Low",
      TRUE ~ "Not significant"
    )
  )

map_b$lisa_class <- factor(
  map_b$lisa_class,
  levels = c("High-High", "Low-Low", "Not significant")
)

pB <- ggplot() +
  geom_sf(data = map_b, aes(fill = lisa_class), color = NA) +
  scale_fill_manual(
    values = c(
      "High-High" = "#ff1a1a",
      "Low-Low" = "#0a7d1c",
      "Not significant" = "#5c5c5c"
    ),
    drop = FALSE
  ) +
  labs(fill = "Cluster", title = "(B)") +
  coord_sf() +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
    axis.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0)
  )

# -------------------------
# 7) PANEL C - SATSCAN
# -------------------------

arq_col_txt <- achar_arquivo(dir_cluster, "\\.col\\.txt$")
tab_c <- ler_col_txt_permutacao(arq_col_txt)

shp_c <- ler_shp_col(dir_cluster)

clusters_c <- shp_c %>%
  left_join(tab_c, by = "cluster")

if (!mostrar_low) {
  clusters_c <- clusters_c %>% filter(cluster_type == "High")
}

if (st_crs(brasil) != st_crs(clusters_c)) {
  brasil_c <- st_transform(brasil, st_crs(clusters_c))
} else {
  brasil_c <- brasil
}

label_c <- clusters_c %>%
  st_centroid() %>%
  mutate(label = paste0("C", cluster))

pC <- ggplot() +
  geom_sf(data = brasil_c, fill = "grey82", color = "grey72", linewidth = 0.05) +
  geom_sf(
    data = clusters_c,
    aes(fill = cluster_type),
    color = NA,
    alpha = 0.25
  ) +
  geom_sf(
    data = clusters_c,
    aes(color = cluster_type),
    fill = NA,
    linewidth = 0.9
  ) +
  geom_sf_text(
    data = label_c,
    aes(label = label, color = cluster_type),
    size = 3,
    fontface = "bold",
    check_overlap = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "High" = "#ff4d4d",
      "Low" = "#4da6ff"
    ),
    guide = "none"
  ) +
  scale_color_manual(
    values = c(
      "High" = "#7a0019",
      "Low" = "#0b3c78"
    ),
    guide = "none"
  ) +
  labs(title = "(C)") +
  coord_sf() +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.3),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0)
  )

# -------------------------
# 8) COMBINE AND EXPORT
# -------------------------

figura2 <- (pA | pB) / pC +
  plot_layout(heights = c(1, 1.25)) &
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  "figura2_ABC.png",
  figura2,
  width = 11,
  height = 13,
  dpi = 500,
  bg = "white"
)

ggsave(
  "figura2_ABC.tiff",
  figura2,
  width = 11,
  height = 13,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  "figura2_ABC.pdf",
  figura2,
  width = 11,
  height = 13,
  bg = "white"
)

print("Arquivos gerados:")
print("figura2_ABC.png")
print("figura2_ABC.tiff")
print("figura2_ABC.pdf")
print(painel)
