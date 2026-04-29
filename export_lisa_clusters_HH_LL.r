# ============================================================
# EXPORTAR CSV DOS CLUSTERS ESPACIAIS LISA
# Saída: lisa_clusters_HH_LL.csv
# Campos:
#   mun_code
#   municipality
#   state
#   taxa_acum
#   lisa_p
#   cluster_type   (High-High / Low-Low)
# Sem geometria
# ============================================================

library(sf)
library(dplyr)
library(readr)
library(stringr)
library(spdep)

# -------------------------
# 1) ENTRADAS
# -------------------------
arq_painel <- "/home/ramos/Documentos/R/TUBERCULOSE/SATSCAN/painel_municipio_ano_ibge7_2004_2023_ok.csv"
arq_brasil <- "/home/ramos/Documentos/R/TUBERCULOSE/IBGE/BR_Municipios_2025/BR_Municipios_2025.shp"

ano_ini <- 2004
ano_fim <- 2023

col_mun   <- "CD_MUN_7"
col_ano   <- "ano"
col_cases <- "casos"
col_pop   <- "populacao"

mult_taxa <- 100000

# -------------------------
# 2) LER DADOS
# -------------------------
painel <- suppressMessages(readr::read_csv(arq_painel))
brasil <- st_read(arq_brasil, quiet = TRUE)

# detectar colunas do shapefile
nm <- names(brasil)

cd_col <- nm[str_detect(tolower(nm), "^cd_mun$|cod.*mun|geocod")]
if (length(cd_col) == 0) stop("Não encontrei coluna do código municipal no shapefile.")
cd_col <- cd_col[1]

mun_name_col <- nm[str_detect(tolower(nm), "^nm_mun$|nome.*mun|municip")]
if (length(mun_name_col) == 0) stop("Não encontrei coluna do nome do município no shapefile.")
mun_name_col <- mun_name_col[1]

uf_col <- nm[str_detect(tolower(nm), "^sigla_uf$|^nm_uf$|uf$|estado")]
if (length(uf_col) == 0) stop("Não encontrei coluna da UF/estado no shapefile.")
uf_col <- uf_col[1]

brasil <- brasil %>%
  mutate(
    mun_code = stringr::str_pad(as.character(.data[[cd_col]]), width = 7, side = "left", pad = "0"),
    municipality = as.character(.data[[mun_name_col]]),
    state = as.character(.data[[uf_col]])
  )

painel2 <- painel %>%
  mutate(
    mun_code   = stringr::str_pad(as.character(.data[[col_mun]]), width = 7, side = "left", pad = "0"),
    year       = as.integer(.data[[col_ano]]),
    cases      = as.numeric(.data[[col_cases]]),
    population = as.numeric(.data[[col_pop]])
  ) %>%
  filter(year >= ano_ini, year <= ano_fim)

# -------------------------
# 3) AGREGAÇÃO MUNICIPAL
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
# 4) LISA / MORAN LOCAL
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
    cluster_type = case_when(
      lisa_p <= 0.05 & x_center > 0 & lag_center > 0 ~ "High-High",
      lisa_p <= 0.05 & x_center < 0 & lag_center < 0 ~ "Low-Low",
      TRUE ~ NA_character_
    )
  )

# -------------------------
# 5) EXPORTAR CSV SEM GEOMETRIA
# -------------------------
saida_csv <- map_b %>%
  st_drop_geometry() %>%
  filter(!is.na(cluster_type)) %>%
  transmute(
    mun_code,
    municipality,
    state,
    cases_sum,
    pop_sum,
    taxa_acum,
    lisa_p,
    cluster_type
  ) %>%
  arrange(cluster_type, state, municipality)

readr::write_csv(saida_csv, "lisa_clusters_HH_LL.csv")

cat("Arquivo gerado:\n")
cat("lisa_clusters_HH_LL.csv\n")
cat("N municípios exportados:", nrow(saida_csv), "\n")
print(head(saida_csv, 10))
