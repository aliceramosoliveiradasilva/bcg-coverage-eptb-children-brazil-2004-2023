# =========================================================
# TB 0-19 anos | Brasil | 2004-2023
# Script para gerar resultados descritivos e texto-base
# =========================================================

library(data.table)
library(fst)
library(openxlsx)
library(glue)

# ---------------------------------------------------------
# 1) Caminhos
# ---------------------------------------------------------
dir_base <- "/home/Ramos/Documentos/R/TUBERCULOSE/DADOS TRATADOS SINAN/"
arq_in   <- file.path(dir_base, "tuberculose_0a19_2004_2023.fst")

arq_out_xlsx <- file.path(dir_base, "resultados_tb_0a19_2004_2023.xlsx")
arq_out_txt  <- file.path(dir_base, "resultados_tb_0a19_2004_2023_texto_auto.txt")

# ---------------------------------------------------------
# 2) Funções auxiliares
# ---------------------------------------------------------
fmt_n <- function(x) {
  formatC(x, format = "d", big.mark = ",")
}

fmt_pct <- function(x, digits = 1) {
  formatC(round(x, digits), format = "f", digits = digits)
}

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) return("< 0.001")
  paste0("= ", formatC(p, format = "f", digits = 3))
}

mode_row <- function(dt, group_var, value_var) {
  tmp <- copy(dt)
  tmp <- tmp[!is.na(get(group_var)) & !is.na(get(value_var))]
  out <- tmp[, .N, by = c(group_var, value_var)]
  setorderv(out, cols = c(group_var, "N"), order = c(1, -1))
  out[, pct := 100 * N / sum(N), by = group_var]
  out[, rk := seq_len(.N), by = group_var]
  out[rk == 1][, rk := NULL]
}

chisq_or_fisher <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  
  tab <- table(x, y)
  
  if (length(x) == 0L || nrow(tab) < 2L || ncol(tab) < 2L) {
    return(data.table(
      test = NA_character_,
      statistic = NA_real_,
      p_value = NA_real_,
      n = sum(ok)
    ))
  }
  
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  
  if (any(chi$expected < 5)) {
    if (nrow(tab) == 2L && ncol(tab) == 2L) {
      ft <- fisher.test(tab)
      return(data.table(
        test = "Fisher's exact test",
        statistic = NA_real_,
        p_value = ft$p.value,
        n = sum(ok)
      ))
    } else {
      chi_sim <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 10000))
      return(data.table(
        test = "Chi-square test with simulated p-value",
        statistic = as.numeric(chi_sim$statistic),
        p_value = chi_sim$p.value,
        n = sum(ok)
      ))
    }
  }
  
  data.table(
    test = "Pearson's chi-square",
    statistic = as.numeric(chi$statistic),
    p_value = chi$p.value,
    n = sum(ok)
  )
}

make_cat_table <- function(dt, var, group_var) {
  x <- copy(dt)
  x[, categoria := as.character(get(var))]
  x[is.na(categoria) | categoria == "", categoria := "Missing"]
  
  out <- x[, .N, by = c(group_var, "categoria")]
  setnames(out, old = c(group_var, "categoria"), new = c("group", "category"))
  out[, pct := 100 * N / sum(N), by = group]
  out[, variable := var]
  setcolorder(out, c("variable", "category", "group", "N", "pct"))
  out[]
}

# ---------------------------------------------------------
# 3) Leitura
# ---------------------------------------------------------
dados <- read_fst(arq_in, as.data.table = TRUE)

# ---------------------------------------------------------
# 4) Recodificação
# ---------------------------------------------------------

# Idade
dados[, idade_unidade := NU_IDADE_N %/% 1000L]
dados[, idade_valor   := NU_IDADE_N %% 1000L]

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

# Forma clínica
dados[, forma_tb := fcase(
  FORMA == 1, "Pulmonary",
  FORMA == 2, "Extrapulmonary",
  FORMA == 3, "Pulmonary + extrapulmonary",
  default = NA_character_
)]

dados[, forma_tb := factor(
  forma_tb,
  levels = c("Pulmonary", "Extrapulmonary", "Pulmonary + extrapulmonary")
)]

# Grupo analítico principal
dados[, grupo_eptb := fcase(
  FORMA == 2, "Exclusive EPTB",
  FORMA == 3, "Pulmonary + EPTB",
  default = NA_character_
)]

dados[, grupo_eptb := factor(
  grupo_eptb,
  levels = c("Exclusive EPTB", "Pulmonary + EPTB")
)]

# Sexo
dados[, sexo := fcase(
  CS_SEXO %in% c("M", "1", 1), "Male",
  CS_SEXO %in% c("F", "2", 2), "Female",
  default = NA_character_
)]

# Raça/cor
dados[, raca := fcase(
  CS_RACA %in% c("1", 1), "White",
  CS_RACA %in% c("2", 2), "Black",
  CS_RACA %in% c("3", 3), "Asian",
  CS_RACA %in% c("4", 4), "Brown",
  CS_RACA %in% c("5", 5), "Indigenous",
  default = NA_character_
)]

# HIV: mantém só positivo/negativo como categorias analíticas
dados[, hiv_bin := fcase(
  HIV %in% c("1", 1), "Positive",
  HIV %in% c("2", 2), "Negative",
  default = NA_character_
)]

# AIDS
dados[, aids_bin := fcase(
  AGRAVAIDS %in% c("1", 1), "Yes",
  AGRAVAIDS %in% c("2", 2), "No",
  default = NA_character_
)]

# Álcool, drogas, tabaco
dados[, alcool := fcase(
  AGRAVALCOO %in% c("1", 1), "Yes",
  AGRAVALCOO %in% c("2", 2), "No",
  default = NA_character_
)]

dados[, drogas := fcase(
  AGRAVDROGA %in% c("1", 1), "Yes",
  AGRAVDROGA %in% c("2", 2), "No",
  default = NA_character_
)]

dados[, tabaco := fcase(
  AGRAVTABAC %in% c("1", 1), "Yes",
  AGRAVTABAC %in% c("2", 2), "No",
  default = NA_character_
)]

# Uso de substâncias: qualquer álcool, droga ou tabaco = Yes
dados[, uso_substancia := fcase(
  alcool == "Yes" | drogas == "Yes" | tabaco == "Yes", "Yes",
  alcool == "No"  & drogas == "No"  & tabaco == "No",  "No",
  default = NA_character_
)]

# Sítio extrapulmonar primário:
# usa EXTRAPU1_N; se vazio, usa EXTRAPU2_N
dados[, sitio_codigo := as.character(EXTRAPU1_N)]
dados[(is.na(sitio_codigo) | sitio_codigo == ""), sitio_codigo := as.character(EXTRAPU2_N)]
dados[sitio_codigo %in% c("", "NA", "<NA>"), sitio_codigo := NA_character_]

# ---------------------------------------------------------
# 5) Bases analíticas
# ---------------------------------------------------------
base_total <- dados[!is.na(forma_tb) & !is.na(faixa_etaria)]
base_eptb  <- dados[FORMA %in% c(2, 3)]
base_comp  <- dados[FORMA %in% c(2, 3) & !is.na(grupo_eptb)]

# ---------------------------------------------------------
# 6) Totais gerais
# ---------------------------------------------------------
n_total_eptb <- nrow(base_eptb)
n_exclusive  <- base_eptb[FORMA == 2, .N]
n_assoc      <- base_eptb[FORMA == 3, .N]

pct_exclusive <- 100 * n_exclusive / n_total_eptb
pct_assoc     <- 100 * n_assoc / n_total_eptb

# ---------------------------------------------------------
# 7) Distribuição por idade e forma clínica
# ---------------------------------------------------------
tab_n <- base_total[, .N, by = .(faixa_etaria, forma_tb)]
tab_linha <- copy(tab_n)
tab_linha[, pct_linha := 100 * N / sum(N), by = faixa_etaria]

tab_coluna <- copy(tab_n)
tab_coluna[, pct_coluna := 100 * N / sum(N), by = forma_tb]

# ---------------------------------------------------------
# 8) Distribuição de sexo e raça/cor entre casos EPTB
# ---------------------------------------------------------
sexo_dist <- base_eptb[!is.na(sexo), .N, by = sexo][order(-N)]
sexo_dist[, pct := 100 * N / sum(N)]

raca_dist <- base_eptb[!is.na(raca), .N, by = raca][order(-N)]
raca_dist[, pct := 100 * N / sum(N)]

top_sexo <- sexo_dist[1]
top_raca <- raca_dist[1]

# ---------------------------------------------------------
# 9) Tabela descritiva por grupo:
#    Exclusive EPTB vs Pulmonary + EPTB
# ---------------------------------------------------------
vars_table1 <- c(
  "sexo", "raca", "faixa_etaria", "hiv_bin",
  "aids_bin", "alcool", "drogas", "tabaco", "uso_substancia"
)

table1_long <- rbindlist(
  lapply(vars_table1, function(v) make_cat_table(base_comp, v, "grupo_eptb")),
  fill = TRUE
)

# p-valores
tests <- rbindlist(list(
  cbind(variable = "sexo",           chisq_or_fisher(base_comp$sexo,           base_comp$grupo_eptb)),
  cbind(variable = "raca",           chisq_or_fisher(base_comp$raca,           base_comp$grupo_eptb)),
  cbind(variable = "faixa_etaria",   chisq_or_fisher(base_comp$faixa_etaria,   base_comp$grupo_eptb)),
  cbind(variable = "hiv_bin",        chisq_or_fisher(base_comp$hiv_bin,        base_comp$grupo_eptb)),
  cbind(variable = "aids_bin",       chisq_or_fisher(base_comp$aids_bin,       base_comp$grupo_eptb)),
  cbind(variable = "alcool",         chisq_or_fisher(base_comp$alcool,         base_comp$grupo_eptb)),
  cbind(variable = "drogas",         chisq_or_fisher(base_comp$drogas,         base_comp$grupo_eptb)),
  cbind(variable = "tabaco",         chisq_or_fisher(base_comp$tabaco,         base_comp$grupo_eptb)),
  cbind(variable = "uso_substancia", chisq_or_fisher(base_comp$uso_substancia, base_comp$grupo_eptb))
), fill = TRUE)

tests[, p_value_fmt := vapply(p_value, fmt_p, FUN.VALUE = character(1))]

# ---------------------------------------------------------
# 10) Sítio anatômico extrapulmonar
# ---------------------------------------------------------
# Primeiro exporta a distribuição por código, sem assumir rótulos.
sitio_codigo_tab <- base_eptb[!is.na(sitio_codigo), .N, by = sitio_codigo][order(-N)]
sitio_codigo_tab[, pct := 100 * N / sum(N)]

# Dicionário para rotular os códigos do sítio anatômico
# PREENCHA DE ACORDO COM O DICIONÁRIO DO SEU BANCO.
# Enquanto não preencher, a aba de sítio ficará por código bruto.
site_lookup <- data.table(
  sitio_codigo = character(),
  sitio_label  = character()
)

if (nrow(site_lookup) > 0) {
  sitio_tab <- merge(
    sitio_codigo_tab,
    site_lookup,
    by = "sitio_codigo",
    all.x = TRUE
  )
  sitio_tab[is.na(sitio_label), sitio_label := paste0("Code ", sitio_codigo)]
  setcolorder(sitio_tab, c("sitio_codigo", "sitio_label", "N", "pct"))
} else {
  sitio_tab <- copy(sitio_codigo_tab)
}

# ---------------------------------------------------------
# 11) Bloco opcional para taxas de detecção
# ---------------------------------------------------------
# Espera um arquivo externo com denominadores populacionais.
# Exemplo UF-ano: colunas "ano", "SG_UF", "pop_0_19"
#
# Se você quiser calcular a mediana e a amplitude das taxas por UF-ano,
# basta descomentar e ajustar o caminho.
#
# arq_pop <- file.path(dir_base, "pop_0_19_uf_ano.csv")
# if (file.exists(arq_pop)) {
#   pop <- fread(arq_pop)
#   casos_uf_ano <- base_eptb[, .N, by = .(ano, SG_UF)]
#   taxas_uf_ano <- merge(casos_uf_ano, pop, by = c("ano", "SG_UF"), all.x = TRUE)
#   taxas_uf_ano[, taxa_100k := 100000 * N / pop_0_19]
#
#   resumo_taxas <- taxas_uf_ano[, .(
#     mediana = median(taxa_100k, na.rm = TRUE),
#     minimo  = min(taxa_100k, na.rm = TRUE),
#     maximo  = max(taxa_100k, na.rm = TRUE)
#   )]
# } else {
#   resumo_taxas <- data.table(
#     mediana = NA_real_,
#     minimo  = NA_real_,
#     maximo  = NA_real_
#   )
# }

# ---------------------------------------------------------
# 12) Extração de resultados-chave para texto
# ---------------------------------------------------------

# Top faixa etária dentro de cada forma clínica
idade_excl <- tab_coluna[forma_tb == "Extrapulmonary"][order(-pct_coluna)][1]
idade_assoc <- tab_coluna[forma_tb == "Pulmonary + extrapulmonary"][order(-pct_coluna)][1]

# p-valores de interesse
p_idade <- tests[variable == "faixa_etaria", p_value_fmt]
p_hiv   <- tests[variable == "hiv_bin", p_value_fmt]
p_subs  <- tests[variable == "uso_substancia", p_value_fmt]

# ---------------------------------------------------------
# 12) Extração de resultados-chave para texto
# ---------------------------------------------------------

# Seleciona diretamente a faixa etária mais frequente dentro de cada forma clínica
idade_excl <- tab_coluna[forma_tb == "Extrapulmonary"][order(-pct_coluna)][1]
idade_assoc <- tab_coluna[forma_tb == "Pulmonary + extrapulmonary"][order(-pct_coluna)][1]

# p-valores de interesse
p_sexo   <- tests[variable == "sexo", p_value_fmt]
p_raca   <- tests[variable == "raca", p_value_fmt]
p_idade  <- tests[variable == "faixa_etaria", p_value_fmt]
p_hiv    <- tests[variable == "hiv_bin", p_value_fmt]
p_aids   <- tests[variable == "aids_bin", p_value_fmt]
p_alcool <- tests[variable == "alcool", p_value_fmt]
p_drogas <- tests[variable == "drogas", p_value_fmt]
p_tabaco <- tests[variable == "tabaco", p_value_fmt]
p_subs   <- tests[variable == "uso_substancia", p_value_fmt]

texto_auto <- c(
  glue(
    "Between 2004 and 2023, {fmt_n(n_total_eptb)} cases of extrapulmonary tuberculosis (EPTB) were reported among children and adolescents aged 0 to 19 years in Brazil, including {fmt_n(n_exclusive)} cases with exclusively extrapulmonary involvement and {fmt_n(n_assoc)} cases with concomitant pulmonary disease."
  ),
  glue(
    "Among EPTB cases with information on sex, males accounted for {fmt_pct(top_sexo$pct)}% (n = {fmt_n(top_sexo$N)})."
  ),
  glue(
    "Among cases with information on race/skin colour, the brown category was the most frequent, accounting for {fmt_pct(raca_dist[raca == 'Brown', pct][1])}% (n = {fmt_n(raca_dist[raca == 'Brown', N][1])})."
  ),
  glue(
    "Regarding age distribution, adolescents aged {idade_excl$faixa_etaria} years accounted for the largest share of exclusive EPTB cases ({fmt_pct(idade_excl$pct_coluna)}%; n = {fmt_n(idade_excl$N)})."
  ),
  glue(
    "A similar pattern was observed among cases with combined pulmonary and extrapulmonary disease, in which adolescents aged {idade_assoc$faixa_etaria} years represented the largest proportion ({fmt_pct(idade_assoc$pct_coluna)}%; n = {fmt_n(idade_assoc$N)})."
  ),
  glue(
    "Statistically significant differences between exclusive EPTB and pulmonary + EPTB were observed for race/skin colour (p {p_raca}), age distribution (p {p_idade}), HIV status (p {p_hiv}), AIDS (p {p_aids}), alcohol use (p {p_alcool}), illicit drug use (p {p_drogas}), tobacco use (p {p_tabaco}), and any reported substance use (p {p_subs}), whereas no significant difference was observed for sex (p {p_sexo})."
  )
)

# ---------------------------------------------------------
# 13) Exportação
# ---------------------------------------------------------
wb <- createWorkbook()

addWorksheet(wb, "overall_summary")
writeDataTable(
  wb, "overall_summary",
  data.table(
    indicator = c(
      "Total EPTB cases",
      "Exclusive EPTB",
      "Pulmonary + EPTB",
      "Exclusive EPTB (%)",
      "Pulmonary + EPTB (%)"
    ),
    value = c(
      n_total_eptb,
      n_exclusive,
      n_assoc,
      round(pct_exclusive, 1),
      round(pct_assoc, 1)
    )
  )
)

addWorksheet(wb, "age_by_clinical_form_n")
writeDataTable(wb, "age_by_clinical_form_n", tab_n)

addWorksheet(wb, "age_by_clinical_form_rowpct")
writeDataTable(wb, "age_by_clinical_form_rowpct", tab_linha)

addWorksheet(wb, "age_by_clinical_form_colpct")
writeDataTable(wb, "age_by_clinical_form_colpct", tab_coluna)

addWorksheet(wb, "sex_distribution_eptb")
writeDataTable(wb, "sex_distribution_eptb", sexo_dist)

addWorksheet(wb, "race_distribution_eptb")
writeDataTable(wb, "race_distribution_eptb", raca_dist)

addWorksheet(wb, "table1_long")
writeDataTable(wb, "table1_long", table1_long)

addWorksheet(wb, "association_tests")
writeDataTable(wb, "association_tests", tests)

addWorksheet(wb, "site_distribution_code")
writeDataTable(wb, "site_distribution_code", sitio_codigo_tab)

addWorksheet(wb, "site_distribution_labeled")
writeDataTable(wb, "site_distribution_labeled", sitio_tab)

addWorksheet(wb, "text_auto")
writeData(wb, "text_auto", data.table(text = texto_auto))

saveWorkbook(wb, arq_out_xlsx, overwrite = TRUE)
writeLines(texto_auto, con = arq_out_txt)

cat("\nArquivos gerados com sucesso:\n")
cat(arq_out_xlsx, "\n")
cat(arq_out_txt,  "\n")