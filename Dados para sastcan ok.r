rm(list = ls(all.names = TRUE))
gc()
closeAllConnections()

library(fst)
library(dplyr)
library(tidyr)
library(readr)
library(sf)

#========================================================
# 0. PATHS
#========================================================
dir_sinan <- "/home/ramos/Documentos/R/TUBERCULOSE/DADOS TRATADOS SINAN"
dir_ibge  <- "/home/ramos/Documentos/R/TUBERCULOSE/IBGE"
dir_out   <- "/home/ramos/Documentos/R/TUBERCULOSE/SATSCAN"
shape_dir <- "/home/ramos/Documentos/R/TUBERCULOSE/IBGE/BR_Municipios_2025"
path_shp  <- file.path(shape_dir, "BR_Municipios_2025.shp")

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

#========================================================
# 1. CASES: read and aggregate by municipality-year
#========================================================
setwd(dir_sinan)

dados <- read_fst("tuberculose 0 a 19 de 2003 a 2024.fst")

cat("\nYears in case file:\n")
print(sort(unique(dados$ano)))

dados2 <- dados %>%
  transmute(
    ano = as.integer(ano),
    # If you want municipality of residence instead, replace ID_MUNICIP by ID_MN_RESI
    CD_MUN_6 = sub("^([0-9]{6}).*$", "\\1", as.character(ID_MUNICIP))
  ) %>%
  filter(
    !is.na(CD_MUN_6),
    CD_MUN_6 != "",
    ano >= 2004,
    ano <= 2023
  )

casos_mun_ano <- dados2 %>%
  count(CD_MUN_6, ano, name = "casos") %>%
  arrange(CD_MUN_6, ano)

cat("\nMunicipalities with cases in 2004-2023:\n")
print(n_distinct(casos_mun_ano$CD_MUN_6))

#========================================================
# 2. SHAPEFILE: official municipal codes and coordinates
#========================================================
muni_sf <- st_read(path_shp, quiet = TRUE)

cat("\nColumns in shapefile:\n")
print(names(muni_sf))

# code column
possiveis_codigos <- c("CD_MUN", "cd_mun", "CODMUN", "cod_mun", "GEOCODIGO", "geocodigo")
col_codigo <- possiveis_codigos[possiveis_codigos %in% names(muni_sf)][1]

if (is.na(col_codigo)) {
  stop("Municipal code column not found in shapefile.")
}

# name column
possiveis_nomes <- c("NM_MUN", "nm_mun", "NOME", "nome", "NM_MUNICIP", "nm_municip")
col_nome <- possiveis_nomes[possiveis_nomes %in% names(muni_sf)][1]

if (is.na(col_nome)) {
  warning("Municipal name column not found. Proceeding without names.")
  muni_sf$NM_MUN_AUX <- NA_character_
  col_nome <- "NM_MUN_AUX"
}

# CRS
if (is.na(st_crs(muni_sf))) {
  stop("Shapefile has no CRS.")
}

if (st_crs(muni_sf)$epsg != 4326) {
  muni_sf <- st_transform(muni_sf, 4326)
}

# safer than centroid for irregular polygons
pts <- st_point_on_surface(muni_sf)
xy  <- st_coordinates(pts)

shape_ref <- muni_sf %>%
  st_drop_geometry() %>%
  transmute(
    CD_MUN_7 = as.character(.data[[col_codigo]]),
    NM_MUN   = as.character(.data[[col_nome]])
  ) %>%
  bind_cols(as.data.frame(xy)) %>%
  rename(
    longitude = X,
    latitude  = Y
  ) %>%
  mutate(
    CD_MUN_7 = sub("^.*?([0-9]{7}).*$", "\\1", CD_MUN_7),
    CD_MUN_6 = substr(CD_MUN_7, 1, 6)
  ) %>%
  distinct(CD_MUN_7, .keep_all = TRUE) %>%
  arrange(CD_MUN_7)

cat("\nMunicipalities in shapefile:\n")
print(nrow(shape_ref))

cat("\nUnique IBGE 7-digit municipal codes:\n")
print(n_distinct(shape_ref$CD_MUN_7))

cat("\nUnique derived DATASUS 6-digit codes:\n")
print(n_distinct(shape_ref$CD_MUN_6))

#========================================================
# 3. POPULATION: wide -> long
#========================================================
pop <- read.csv(
  file.path(dir_ibge, "pop0a19.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\nFirst population columns:\n")
print(names(pop)[1:min(10, ncol(pop))])

pop_long <- pop %>%
  mutate(
    CD_MUN_6 = sub("^([0-9]{6}).*$", "\\1", as.character(Mun))
  ) %>%
  select(CD_MUN_6, matches("^X?\\d{4}$")) %>%
  pivot_longer(
    cols = matches("^X?\\d{4}$"),
    names_to = "ano",
    values_to = "populacao"
  ) %>%
  mutate(
    ano = as.integer(sub("^X", "", ano)),
    CD_MUN_6 = as.character(CD_MUN_6),
    populacao = as.numeric(populacao)
  ) %>%
  filter(
    !is.na(CD_MUN_6),
    CD_MUN_6 != "",
    ano >= 2004,
    ano <= 2023
  ) %>%
  distinct(CD_MUN_6, ano, .keep_all = TRUE) %>%
  arrange(CD_MUN_6, ano)

#========================================================
# 4. FINAL PANEL
#========================================================
painel <- expand_grid(
  CD_MUN_7 = sort(unique(shape_ref$CD_MUN_7)),
  ano = 2004:2023
) %>%
  left_join(
    shape_ref %>% select(CD_MUN_7, CD_MUN_6, NM_MUN),
    by = "CD_MUN_7"
  ) %>%
  left_join(
    pop_long,
    by = c("CD_MUN_6", "ano")
  ) %>%
  left_join(
    casos_mun_ano,
    by = c("CD_MUN_6", "ano")
  ) %>%
  mutate(
    casos = if_else(is.na(casos), 0L, as.integer(casos))
  ) %>%
  arrange(CD_MUN_7, ano)

#========================================================
# 5. CHECKS
#========================================================
cat("\nPanel year range:\n")
print(range(painel$ano, na.rm = TRUE))

cat("\nYears in panel:\n")
print(table(painel$ano))

cat("\nUnique municipalities in panel (IBGE 7-digit):\n")
print(n_distinct(painel$CD_MUN_7))

cat("\nMissing population values:\n")
print(sum(is.na(painel$populacao)))

cat("\nMissing case values:\n")
print(sum(is.na(painel$casos)))

sem_pop <- painel %>%
  filter(is.na(populacao) | populacao <= 0) %>%
  distinct(CD_MUN_7, CD_MUN_6, NM_MUN)

cat("\nMunicipalities with invalid population in any year:\n")
print(sem_pop, n = Inf)

casos_sem_shape <- casos_mun_ano %>%
  anti_join(shape_ref %>% distinct(CD_MUN_6), by = "CD_MUN_6")

cat("\nCase rows without matching shapefile municipality:\n")
print(nrow(casos_sem_shape))

#========================================================
# 6. IDENTIFY AND REMOVE INVALID LOCALITIES FOR SATSCAN
#========================================================
bad_pop <- painel %>%
  transmute(
    localidade = CD_MUN_7,
    ano = ano,
    populacao = populacao
  ) %>%
  filter(is.na(populacao) | populacao <= 0)

bad_codes <- bad_pop %>%
  distinct(localidade)

cat("\nLocalities with invalid population for SaTScan:\n")
print(bad_codes, n = Inf)

# inspect total cases in problematic localities
casos_bad <- painel %>%
  filter(CD_MUN_7 %in% bad_codes$localidade) %>%
  group_by(CD_MUN_7, NM_MUN) %>%
  summarise(
    casos_totais = sum(casos, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nTotal cases in problematic localities:\n")
print(casos_bad, n = Inf)

# keep only valid municipality-years for population
painel_ok <- painel %>%
  filter(
    !CD_MUN_7 %in% bad_codes$localidade,
    !is.na(populacao),
    populacao > 0
  )

coord_csv_ok <- shape_ref %>%
  filter(!CD_MUN_7 %in% bad_codes$localidade) %>%
  select(CD_MUN_7, CD_MUN_6, NM_MUN, latitude, longitude)

cat("\nUnique municipalities after cleaning:\n")
print(n_distinct(painel_ok$CD_MUN_7))

cat("\nRemaining missing population values after cleaning:\n")
print(sum(is.na(painel_ok$populacao)))

#========================================================
# 7. EXPORT CLEAN CSV FILES
#========================================================
painel_csv_ok <- painel_ok %>%
  select(CD_MUN_7, CD_MUN_6, NM_MUN, ano, casos, populacao)

write_csv(
  painel_csv_ok,
  file.path(dir_out, "painel_municipio_ano_ibge7_2004_2023_ok.csv")
)

write_csv(
  coord_csv_ok,
  file.path(dir_out, "coordenadas_municipios_ibge7_ok.csv")
)

#========================================================
# 8. EXPORT CLEAN SATSCAN FILES
#========================================================
geo_satscan_ok <- coord_csv_ok %>%
  select(CD_MUN_7, latitude, longitude)

cas_satscan_ok <- painel_ok %>%
  transmute(
    localidade = CD_MUN_7,
    ano = ano,
    casos = casos
  )

pop_satscan_ok <- painel_ok %>%
  transmute(
    localidade = CD_MUN_7,
    ano = ano,
    populacao = populacao
  )

# final consistency checks
bad_pop_final <- pop_satscan_ok %>%
  filter(is.na(populacao) | populacao <= 0)

cat("\nInvalid population rows after cleaning:\n")
print(nrow(bad_pop_final))

cat("\nUnique localities in GEO:\n")
print(n_distinct(geo_satscan_ok$CD_MUN_7))

cat("\nUnique localities in CAS:\n")
print(n_distinct(cas_satscan_ok$localidade))

cat("\nUnique localities in POP:\n")
print(n_distinct(pop_satscan_ok$localidade))

write.table(
  geo_satscan_ok,
  file = file.path(dir_out, "coordinates_ibge7_2004_2023_ok.geo"),
  sep = " ",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  cas_satscan_ok,
  file = file.path(dir_out, "cases_ibge7_2004_2023_ok.cas"),
  sep = " ",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  pop_satscan_ok,
  file = file.path(dir_out, "population_ibge7_2004_2023_ok.pop"),
  sep = " ",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  na = ""
)

#========================================================
# 9. FINAL SUMMARY
#========================================================
cat("\n====================================\n")
cat("Files generated in:", dir_out, "\n")
cat("Clean panel CSV: ", file.path(dir_out, "painel_municipio_ano_ibge7_2004_2023_ok.csv"), "\n", sep = "")
cat("Clean coordinates CSV: ", file.path(dir_out, "coordenadas_municipios_ibge7_ok.csv"), "\n", sep = "")
cat("Clean CAS: ", file.path(dir_out, "cases_ibge7_2004_2023_ok.cas"), "\n", sep = "")
cat("Clean POP: ", file.path(dir_out, "population_ibge7_2004_2023_ok.pop"), "\n", sep = "")
cat("Clean GEO: ", file.path(dir_out, "coordinates_ibge7_2004_2023_ok.geo"), "\n", sep = "")
cat("Removed localities:", nrow(bad_codes), "\n")
cat("Remaining municipalities:", n_distinct(painel_ok$CD_MUN_7), "\n")
cat("Total cases in clean panel:", sum(painel_ok$casos, na.rm = TRUE), "\n")
cat("====================================\n")