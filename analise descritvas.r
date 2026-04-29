library(data.table)
library(fst)

setwd("/home/Ramos/Documentos/R/TUBERCULOSE/DADOS TRATADOS SINAN/")

dados <- read_fst("tuberculose_0a19_2004_2023.fst", as.data.table = TRUE)

# Recriar componentes da idade, se necessário
dados[, idade_unidade := NU_IDADE_N %/% 1000L]
dados[, idade_valor   := NU_IDADE_N %% 1000L]

# Criar faixa etária
# Tudo que estiver em hora/dia/mês entra em 0-4 anos
dados[, faixa_etaria := fcase(
  idade_unidade %in% c(1L, 2L, 3L), "0-4",
  idade_unidade == 4L & idade_valor >= 0L  & idade_valor <= 4L,  "0-4",
  idade_unidade == 4L & idade_valor >= 5L  & idade_valor <= 9L,  "5-9",
  idade_unidade == 4L & idade_valor >= 10L & idade_valor <= 14L, "10-14",
  idade_unidade == 4L & idade_valor >= 15L & idade_valor <= 19L, "15-19",
  default = NA_character_
)]

dados[, faixa_etaria := factor(
  faixa_etaria,
  levels = c("0-4", "5-9", "10-14", "15-19")
)]

# Recodificar forma clínica
dados[, forma_tb := fcase(
  FORMA == 1, "Pulmonar",
  FORMA == 2, "Extrapulmonar",
  FORMA == 3, "Pulmonar + Extrapulmonar",
  default = NA_character_
)]

dados[, forma_tb := factor(
  forma_tb,
  levels = c("Pulmonar", "Extrapulmonar", "Pulmonar + Extrapulmonar")
)]

# Manter apenas registros válidos para a tabela
tab_base <- dados[
  !is.na(faixa_etaria) & !is.na(forma_tb)
]

# Tabela de frequências absolutas
tab_n <- tab_base[, .N, by = .(faixa_etaria, forma_tb)]

# Tabela com percentuais por faixa etária (linha)
tab_linha <- tab_base[, .N, by = .(faixa_etaria, forma_tb)]
tab_linha[, perc_linha := 100 * N / sum(N), by = faixa_etaria]

# Tabela com percentuais por forma clínica (coluna)
tab_coluna <- tab_base[, .N, by = .(faixa_etaria, forma_tb)]
tab_coluna[, perc_coluna := 100 * N / sum(N), by = forma_tb]

# Ver
tab_n[order(faixa_etaria, forma_tb)]
tab_linha[order(faixa_etaria, forma_tb)]
tab_coluna[order(faixa_etaria, forma_tb)]
colnames(dados)
