rm(list = ls())
gc()
library(read.dbc)
library(data.table)
library(stringr)
library(fst)
setwd("/home/Ramos/Documentos/R/TUBERCULOSE/DADOS TRATADOS SINAN/")
dados <- read_fst("tuberculose_sinan_2003_2024.fst")
colnames(dados)
unique(dados$NU_IDADE_N)

setDT(dados)

dados[, idade_unidade := fifelse(!is.na(NU_IDADE_N), NU_IDADE_N %/% 1000L, NA_integer_)]
dados[, idade_valor   := fifelse(!is.na(NU_IDADE_N), NU_IDADE_N %% 1000L, NA_integer_)]

# filtro: menores de 20 anos
dados_menor20 <- dados[
  !is.na(idade_unidade) & !is.na(idade_valor) &
    (
      idade_unidade %in% c(1L, 2L, 3L) |
        (idade_unidade == 4L & idade_valor < 20L)
    )
]

# checagem
dados_menor20[, .N, by = idade_unidade][order(idade_unidade)]

dados[, idade_label := fifelse(
  is.na(NU_IDADE_N), NA_character_,
  paste0(idade_valor, " ",
         fifelse(idade_unidade == 1, "hora(s)",
                 fifelse(idade_unidade == 2, "dia(s)",
                         fifelse(idade_unidade == 3, "mes(es)",
                                 fifelse(idade_unidade == 4, "ano(s)", "inválido")))))
)]


write_fst(dados_menor20, "tuberculose 0 a 19 de 2003 a 2024.fst")
