rm(list = ls())
gc()
library(read.dbc)
library(data.table)
library(stringr)
library(fst)
# Caminho da pasta com os arquivos .dbc
pasta <- "/home/Ramos/Documentos/R/TUBERCULOSE/DADOS BRUTOS SINAN/"

# Listar todos os .dbc
arquivos <- list.files(
  path = pasta,
  pattern = "\\.dbc$",
  full.names = TRUE
)

arquivos <- sort(arquivos)

# Função para ler cada arquivo e adicionar metadados
ler_um_dbc <- function(arq) {
  cat("Lendo:", basename(arq), "\n")
  
  dt <- as.data.table(read.dbc::read.dbc(arq))
  
  # Extrair ano do nome do arquivo, ex: TUBEBR04.dbc -> 2004
  yy <- as.integer(str_extract(basename(arq), "\\d{2}(?=\\.dbc$)"))
  ano <- ifelse(!is.na(yy), 2000L + yy, NA_integer_)
  
  dt[, arquivo_origem := basename(arq)]
  dt[, ano := ano]
  
  return(dt)
}

# Ler todos os arquivos em uma lista nomeada
bases_lista <- lapply(arquivos, ler_um_dbc)
names(bases_lista) <- tools::file_path_sans_ext(basename(arquivos))

# Unir todos em uma base só
base_tuberculose <- rbindlist(
  bases_lista,
  use.names = TRUE,
  fill = TRUE,
  idcol = "fonte"
)

# Verificação rápida
dim(base_tuberculose)
names(base_tuberculose)
table(base_tuberculose$ano, useNA = "ifany")

# Salvar em formato rápido
fwrite(base_tuberculose, file.path(pasta, "tuberculose_sinan_2003_2024.csv"))
saveRDS(base_tuberculose, file.path(pasta, "tuberculose_sinan_2003_2024.rds"))
write_fst(base_tuberculose, "tuberculose_sinan_2003_2024.fst")
