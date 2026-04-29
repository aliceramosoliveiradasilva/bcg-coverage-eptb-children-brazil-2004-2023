library(data.table)

dir_cnes <- "/home/Ramos/Documentos/R/TUBERCULOSE/CAPACIDADE DIAGNOSTICA"

arquivos <- sort(list.files(
  path = dir_cnes,
  pattern = "\\.csv$",
  full.names = TRUE
))

safe_utf8 <- function(x) {
  y <- iconv(x, from = "Latin1", to = "UTF-8", sub = "")
  y[is.na(y)] <- ""
  y
}

ler_cnes_ano_long <- function(arq) {
  cat("Lendo:", basename(arq), "\n")
  
  con <- file(arq, open = "r", encoding = "Latin1")
  on.exit(close(con), add = TRUE)
  
  linhas <- readLines(con, warn = FALSE)
  linhas <- safe_utf8(linhas)
  linhas <- linhas[nzchar(trimws(linhas))]
  
  # localizar início da tabela
  n_sep <- vapply(
    gregexpr(";", linhas, fixed = TRUE, useBytes = TRUE),
    function(z) if (length(z) == 1L && z[1] == -1L) 0L else length(z),
    integer(1)
  )
  
  idx_header <- which(n_sep >= 5)[1]
  if (is.na(idx_header)) {
    stop(sprintf("Não consegui identificar o cabeçalho em: %s", arq))
  }
  
  linhas_tab <- linhas[idx_header:length(linhas)]
  
  # remover rodapé
  idx_fonte <- grep("^Fonte", linhas_tab, ignore.case = TRUE)
  if (length(idx_fonte) > 0) {
    linhas_tab <- linhas_tab[seq_len(idx_fonte[1] - 1)]
  }
  
  txt <- paste(linhas_tab, collapse = "\n")
  
  dt <- fread(
    text = txt,
    sep = ";",
    fill = TRUE,
    quote = "",
    header = TRUE,
    colClasses = "character",
    na.strings = c("", "NA", "...")
  )
  
  # manter apenas 13 colunas esperadas: municipio + 12 meses
  dt <- dt[, seq_len(min(13, ncol(dt))), with = FALSE]
  
  meses <- c("jan","fev","mar","abr","mai","jun","jul","ago","set","out","nov","dez")
  setnames(dt, c("municipio_raw", meses)[seq_len(ncol(dt))])
  
  # limpar municipio
  dt[, municipio_raw := safe_utf8(municipio_raw)]
  dt[, municipio_raw := gsub('^"|"$', "", municipio_raw)]
  dt[, municipio_raw := trimws(municipio_raw)]
  
  # extrair código e nome
  dt[, cod_mun6 := substr(municipio_raw, 1, 6)]
  dt[, nome_municipio := trimws(sub("^[0-9]{6}\\s+", "", municipio_raw))]
  
  # ano pelo nome do arquivo
  ano_arq <- as.integer(tools::file_path_sans_ext(basename(arq)))
  dt[, ano := ano_arq]
  
  # passar para longo
  dt_long <- melt(
    dt,
    id.vars = c("cod_mun6", "nome_municipio", "ano"),
    measure.vars = meses,
    variable.name = "mes",
    value.name = "n_estab"
  )
  
  # tratar ausências
  dt_long[, n_estab := trimws(n_estab)]
  dt_long[n_estab %chin% c("-", ""), n_estab := NA_character_]
  dt_long[, n_estab := as.integer(n_estab)]
  
  mapa_mes <- c(
    jan = 1L, fev = 2L, mar = 3L, abr = 4L, mai = 5L, jun = 6L,
    jul = 7L, ago = 8L, set = 9L, out = 10L, nov = 11L, dez = 12L
  )
  
  dt_long[, mes_num := unname(mapa_mes[mes])]
  dt_long[, competencia := sprintf("%d-%02d", ano, mes_num)]
  
  setcolorder(dt_long, c("cod_mun6", "nome_municipio", "ano", "mes", "mes_num", "competencia", "n_estab"))
  
  dt_long[]
}

cnes_diag_long <- rbindlist(
  lapply(arquivos, ler_cnes_ano_long),
  use.names = TRUE,
  fill = TRUE
)

library(data.table)

setDT(cnes_diag_long)

cnes_diag_mediana_ano_mun <- cnes_diag_long[
  ,
  .(
    mediana_anual = if (all(is.na(n_estab))) {
      NA_real_
    } else {
      as.numeric(median(n_estab, na.rm = TRUE))
    },
    n_meses_validos = as.integer(sum(!is.na(n_estab)))
  ),
  by = .(cod_mun6, nome_municipio, ano)
]

setorder(cnes_diag_mediana_ano_mun, ano, cod_mun6)

print(dim(cnes_diag_mediana_ano_mun))
print(head(cnes_diag_mediana_ano_mun, 20))

fwrite(
  cnes_diag_mediana_ano_mun,
  file = file.path(
    dir_cnes,
    "cnes_capacidade_diagnostica_mediana_anual_municipio_2006_2024.csv"
  ),
  sep = ",",
  na = ""
)

saveRDS(
  cnes_diag_mediana_ano_mun,
  file = file.path(
    dir_cnes,
    "cnes_capacidade_diagnostica_mediana_anual_municipio_2006_2024.rds"
  )
)