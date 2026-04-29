rm(list = ls(all.names = TRUE))
gc()
closeAllConnections()

library(fst)
library(dplyr)
library(tidyr)
library(readr)
library(sf)

#========================================================
# 0. CAMINHOS
#========================================================
dir_sinan <- "/home/ramos/Documentos/R/TUBERCULOSE/DADOS TRATADOS SINAN"
dir_ibge  <- "/home/ramos/Documentos/R/TUBERCULOSE/IBGE"
dir_out   <- "/home/ramos/Documentos/R/TUBERCULOSE/SATSCAN"
shape <- "/home/ramos/Documentos/R/TUBERCULOSE/IBGE/BR_Municipios_2025/"
# ajuste para o nome real do shapefile
path_shp  <- file.path(shape, "BR_Municipios_2025.shp")

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

#========================================================
# 1. CASOS: ler e resumir por município-ano
#========================================================
setwd(dir_sinan)
dados <- read_fst("tuberculose 0 a 19 de 2003 a 2024.fst")

cat("\nAnos na base de casos:\n")
print(sort(unique(dados$ano)))

dados2 <- dados %>%
  transmute(
    ano = as.integer(ano),
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

cat("\nMunicípios com casos no período:\n")
print(n_distinct(casos_mun_ano$CD_MUN_6))

#========================================================
# 2. SHAPEFILE: ler códigos oficiais e coordenadas
#========================================================
muni_sf <- st_read(path_shp, quiet = TRUE)

cat("\nColunas do shapefile:\n")
print(names(muni_sf))

# identificar coluna do código do município
possiveis_codigos <- c("CD_MUN", "cd_mun", "CODMUN", "cod_mun", "GEOCODIGO", "geocodigo")
col_codigo <- possiveis_codigos[possiveis_codigos %in% names(muni_sf)][1]

if (is.na(col_codigo)) {
  stop("Não encontrei a coluna do código municipal no shapefile.")
}

# identificar coluna do nome do município
possiveis_nomes <- c("NM_MUN", "nm_mun", "NOME", "nome", "NM_MUNICIP", "nm_municip")
col_nome <- possiveis_nomes[possiveis_nomes %in% names(muni_sf)][1]

if (is.na(col_nome)) {
  warning("Não encontrei coluna de nome do município no shapefile. Vou seguir sem nome.")
  muni_sf$NM_MUN_AUX <- NA_character_
  col_nome <- "NM_MUN_AUX"
}

# garantir CRS geográfico
if (is.na(st_crs(muni_sf))) {
  stop("O shapefile está sem CRS definido.")
}

if (st_crs(muni_sf)$epsg != 4326) {
  muni_sf <- st_transform(muni_sf, 4326)
}

# calcular centroides
cent <- st_centroid(muni_sf)
xy <- st_coordinates(cent)

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

cat("\nMunicípios no shapefile:\n")
print(nrow(shape_ref))

cat("\nIBGE 7 dígitos únicos:\n")
print(n_distinct(shape_ref$CD_MUN_7))

cat("\nDATASUS 6 dígitos derivados do shape:\n")
print(n_distinct(shape_ref$CD_MUN_6))

#========================================================
# 3. POPULAÇÃO: wide -> long
#========================================================
pop <- read.csv(
  file.path(dir_ibge, "pop0a19.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cat("\nPrimeiras colunas da população:\n")
print(names(pop)[1:min(10, ncol(pop))])

pop_long <- pop %>%
  mutate(
    CD_MUN_6 = sub("^([0-9]{6}).*$", "\\1", Mun)
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
# 4. PAINEL FINAL
#    shape -> população -> casos
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
    casos = ifelse(is.na(casos), 0L, casos)
  ) %>%
  arrange(CD_MUN_7, ano)

#========================================================
# 5. CHECAGENS
#========================================================
cat("\nFaixa de anos do painel:\n")
print(range(painel$ano, na.rm = TRUE))

cat("\nTabela de anos do painel:\n")
print(table(painel$ano))

cat("\nMunicípios únicos no painel (IBGE 7):\n")
print(n_distinct(painel$CD_MUN_7))

cat("\nNAs em população:\n")
print(sum(is.na(painel$populacao)))

cat("\nNAs em casos:\n")
print(sum(is.na(painel$casos)))

sem_pop <- painel %>%
  filter(is.na(populacao)) %>%
  distinct(CD_MUN_7, CD_MUN_6, NM_MUN)

cat("\nMunicípios com pelo menos um ano sem população:\n")
print(nrow(sem_pop))

casos_sem_shape <- casos_mun_ano %>%
  anti_join(shape_ref %>% distinct(CD_MUN_6), by = "CD_MUN_6")

cat("\nLinhas de casos sem correspondência no shapefile:\n")
print(nrow(casos_sem_shape))

#========================================================
# 6. EXPORTAR CSV FINAL SEM COORDENADAS
#========================================================
painel_csv <- painel %>%
  select(
    CD_MUN_7,
    CD_MUN_6,
    NM_MUN,
    ano,
    casos,
    populacao
  )

write_csv(
  painel_csv,
  file.path(dir_out, "painel_municipio_ano_ibge7_2004_2023.csv")
)

#========================================================
# 7. EXPORTAR CSV DE COORDENADAS
#========================================================
coord_csv <- shape_ref %>%
  select(
    CD_MUN_7,
    CD_MUN_6,
    NM_MUN,
    latitude,
    longitude
  )

write_csv(
  coord_csv,
  file.path(dir_out, "coordenadas_municipios_ibge7.csv")
)

#========================================================
# 8. EXPORTAR ARQUIVOS SATSCAN
#========================================================
geo_satscan <- coord_csv %>%
  select(CD_MUN_7, latitude, longitude)

write.table(
  geo_satscan,
  file = file.path(dir_out, "coordinates_ibge7.geo"),
  sep = " ",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

cas_satscan <- painel %>%
  transmute(
    localidade = CD_MUN_7,
    ano = ano,
    casos = casos
  )

write.table(
  cas_satscan,
  file = file.path(dir_out, "cases_ibge7_2004_2023.cas"),
  sep = " ",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

pop_satscan <- painel %>%
  transmute(
    localidade = CD_MUN_7,
    ano = ano,
    populacao = populacao
  )

write.table(
  pop_satscan,
  file = file.path(dir_out, "population_ibge7_2004_2023.pop"),
  sep = " ",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

#========================================================
# 9. RESUMO FINAL
#========================================================
cat("\n====================================\n")
cat("Arquivos gerados em:", dir_out, "\n")
cat("Painel CSV: ", file.path(dir_out, "painel_municipio_ano_ibge7_2004_2023.csv"), "\n", sep = "")
cat("Coordenadas CSV: ", file.path(dir_out, "coordenadas_municipios_ibge7.csv"), "\n", sep = "")
cat("CAS: ", file.path(dir_out, "cases_ibge7_2004_2023.cas"), "\n", sep = "")
cat("POP: ", file.path(dir_out, "population_ibge7_2004_2023.pop"), "\n", sep = "")
cat("GEO: ", file.path(dir_out, "coordinates_ibge7.geo"), "\n", sep = "")
cat("Municípios no painel:", n_distinct(painel$CD_MUN_7), "\n")
cat("Casos totais:", sum(painel$casos, na.rm = TRUE), "\n")
cat("====================================\n")

bad_pop <- pop_satscan[is.na(pop_satscan$populacao) | pop_satscan$populacao <= 0, ]
bad_codes <- bad_pop %>%
  distinct(localidade)

print(bad_codes, n = Inf)

cas_satscan %>%
  filter(localidade %in% bad_codes$localidade) %>%
  count(localidade, ano, sort = TRUE)

geo_satscan %>%
  filter(CD_MUN_7 %in% bad_codes$localidade)
