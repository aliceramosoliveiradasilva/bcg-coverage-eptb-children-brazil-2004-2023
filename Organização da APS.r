# =========================================================
# APS - Equipes de Saúde por município e mês
# Leitura de 2008 a 2024, padronização e base final em long
# =========================================================

rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

setwd("/home/Ramos/Documentos/R/TUBERCULOSE/APS/")

# ---------------------------------------------------------
# 1) Definir anos e checar arquivos
# ---------------------------------------------------------
anos <- 2008:2024
arquivos <- paste0(anos, ".csv")

faltando <- arquivos[!file.exists(arquivos)]
if (length(faltando) > 0) {
  stop("Arquivos ausentes: ", paste(faltando, collapse = ", "))
}

# ---------------------------------------------------------
# 2) Função para ler e padronizar cada arquivo
# ---------------------------------------------------------
ler_aps_ano <- function(ano) {
  
  arq <- paste0(ano, ".csv")
  
  df <- read.csv2(
    arq,
    skip = 4,
    fileEncoding = "Latin1",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "-", "NA")
  )
  
  # padronizar nome da 1a coluna
  names(df)[1] <- "municipio_raw"
  
  # transformar para long
  df_long <- df %>%
    pivot_longer(
      cols = -municipio_raw,
      names_to = "competencia",
      values_to = "n_equipes"
    ) %>%
    mutate(
      n_equipes = as.numeric(n_equipes),
      ano = as.integer(str_extract(competencia, "^\\d{4}")),
      mes_abrev = str_extract(competencia, "(?<=/)\\w+"),
      mes = match(
        mes_abrev,
        c("Jan", "Fev", "Mar", "Abr", "Mai", "Jun",
          "Jul", "Ago", "Set", "Out", "Nov", "Dez")
      ),
      cod_mun = str_extract(municipio_raw, "^\\d{6}"),
      municipio = str_trim(str_remove(municipio_raw, "^\\d{6}\\s+")),
      uf = substr(cod_mun, 1, 2),
      data_competencia = as.Date(sprintf("%d-%02d-01", ano, mes)),
      ano_arquivo = ano
    ) %>%
    select(
      cod_mun, uf, municipio, ano, mes, mes_abrev,
      data_competencia, n_equipes, ano_arquivo, municipio_raw
    )
  
  # sanity check: o ano do nome da coluna deve bater com o arquivo
  anos_encontrados <- unique(df_long$ano)
  if (!all(anos_encontrados == ano)) {
    warning(
      "O arquivo ", arq,
      " contém colunas de ano diferentes do esperado: ",
      paste(anos_encontrados, collapse = ", ")
    )
  }
  
  return(df_long)
}

# ---------------------------------------------------------
# 3) Ler todos os anos e juntar
# ---------------------------------------------------------
aps_todos <- map_dfr(anos, ler_aps_ano)

# ---------------------------------------------------------
# 4) Ordenar e inspecionar
# ---------------------------------------------------------
aps_todos <- aps_todos %>%
  arrange(cod_mun, ano, mes)

glimpse(aps_todos)
head(aps_todos)
tail(aps_todos)

# ---------------------------------------------------------
# 5) Checagens de consistência
# ---------------------------------------------------------

# número de registros por ano
tab_ano <- aps_todos %>%
  count(ano)

print(tab_ano)

# número de registros por ano e mês
tab_ano_mes <- aps_todos %>%
  count(ano, mes) %>%
  arrange(ano, mes)

print(tab_ano_mes)

# checar duplicatas por município-ano-mês
dups <- aps_todos %>%
  count(cod_mun, ano, mes) %>%
  filter(n > 1)

print(dups)

# checar quantos meses por município-ano
tab_meses_por_mun_ano <- aps_todos %>%
  group_by(cod_mun, municipio, ano) %>%
  summarise(
    n_meses = sum(!is.na(mes)),
    .groups = "drop"
  )

print(summary(tab_meses_por_mun_ano$n_meses))

# ---------------------------------------------------------
# 6) Base anual: média de equipes no ano por município
# ---------------------------------------------------------
aps_anual_media <- aps_todos %>%
  group_by(cod_mun, uf, municipio, ano) %>%
  summarise(
    media_equipes = mean(n_equipes, na.rm = TRUE),
    mediana_equipes = median(n_equipes, na.rm = TRUE),
    min_equipes = min(n_equipes, na.rm = TRUE),
    max_equipes = max(n_equipes, na.rm = TRUE),
    n_meses_validos = sum(!is.na(n_equipes)),
    .groups = "drop"
  )

head(aps_anual_media)

# ---------------------------------------------------------
# 7) Base anual: valor de dezembro por município
#    útil se você quiser estoque ao fim do ano
# ---------------------------------------------------------
aps_anual_dezembro <- aps_todos %>%
  filter(mes == 12) %>%
  select(cod_mun, uf, municipio, ano, equipes_dezembro = n_equipes)

head(aps_anual_dezembro)

# ---------------------------------------------------------
# 8) Opcional: base pronta para merge com TB
#    mantendo apenas município-ano
# ---------------------------------------------------------
aps_merge_tb <- aps_anual_media %>%
  mutate(
    cod_mun = substr(cod_mun, 1, 6),
    ano = as.integer(ano)
  )

head(aps_merge_tb)

# ---------------------------------------------------------
# 9) Salvar
# ---------------------------------------------------------
write.csv(
  aps_todos,
  "aps_equipes_2008_2024_long.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  aps_anual_media,
  "aps_equipes_2008_2024_anual_media.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  aps_anual_dezembro,
  "aps_equipes_2008_2024_anual_dezembro.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


