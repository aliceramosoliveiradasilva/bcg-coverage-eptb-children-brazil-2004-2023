# =========================================================
# BCG coverage and pediatric EPTB in Brazil
# Full revised script to address reviewers' comments
# =========================================================

rm(list = ls())
gc()

# ---------------------------------------------------------
# 0) Packages
# ---------------------------------------------------------
library(data.table)
library(fst)
library(MASS)

# ---------------------------------------------------------
# 1) Helper functions
# ---------------------------------------------------------

clean_names <- function(x) {
  x <- sub("^\ufeff", "", x)
  x <- sub("^ï..", "", x)
  trimws(x)
}

extract_code_from_label <- function(x) {
  x <- trimws(as.character(x))
  out <- substr(x, 1, 6)
  out[!grepl("^[0-9]{6}$", out)] <- NA_character_
  out
}

normalize_code6 <- function(x) {
  y <- suppressWarnings(as.integer(as.character(x)))
  ifelse(!is.na(y) & y >= 100000 & y <= 999999, sprintf("%06d", y), NA_character_)
}

to_num <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

to_num2 <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", ".", "..")] <- NA
  suppressWarnings(as.numeric(x))
}

fix_first_row_as_header <- function(dt) {
  dt <- copy(dt)
  setDT(dt)
  
  new_names <- as.character(unlist(dt[1]))
  setnames(dt, old = names(dt), new = new_names)
  
  dt <- dt[-1]
  names(dt) <- clean_names(names(dt))
  
  year_cols <- grep("^(19|20)\\d{2}$", names(dt), value = TRUE)
  if (length(year_cols) > 0) {
    setnames(dt, old = year_cols, new = paste0("X", year_cols))
  }
  
  dt[]
}

maybe_fix_header <- function(dt) {
  dt <- copy(dt)
  setDT(dt)
  names(dt) <- clean_names(names(dt))
  
  if (all(grepl("^V\\d+$", names(dt)))) {
    return(fix_first_row_as_header(dt))
  }
  
  first_row <- as.character(unlist(dt[1]))
  if (length(first_row) >= 2 &&
      first_row[1] == "Muni" &&
      any(grepl("^(19|20)\\d{2}$", first_row[-1]))) {
    return(fix_first_row_as_header(dt))
  }
  
  dt[]
}

wide_to_long_codeyear <- function(dt, value_name) {
  dt <- copy(dt)
  setDT(dt)
  names(dt) <- clean_names(names(dt))
  
  muni_col <- names(dt)[1]
  dt[, codigo := extract_code_from_label(dt[[muni_col]])]
  
  year_cols <- grep("^X\\d{4}$", names(dt), value = TRUE)
  if (length(year_cols) == 0) {
    stop("Nenhuma coluna de ano no formato X2004, X2005... foi encontrada.")
  }
  
  dt[, (year_cols) := lapply(.SD, to_num2), .SDcols = year_cols]
  
  out <- melt(
    dt,
    id.vars = "codigo",
    measure.vars = year_cols,
    variable.name = "ano",
    value.name = value_name
  )
  
  out[, ano := as.integer(sub("^X", "", ano))]
  out <- out[!is.na(codigo)]
  
  out[]
}

recode_yes_no <- function(x) {
  x <- to_num(x)
  fifelse(x == 1L, 1L, fifelse(x == 2L, 0L, NA_integer_))
}

or_table <- function(model) {
  est <- coef(model)
  se  <- sqrt(diag(vcov(model)))
  
  data.table(
    term = names(est),
    OR   = exp(est),
    LCL  = exp(est - 1.96 * se),
    UCL  = exp(est + 1.96 * se),
    p_value = 2 * pnorm(abs(est / se), lower.tail = FALSE)
  )
}

irr_table_nb <- function(model, var_name) {
  b  <- coef(model)[var_name]
  se <- sqrt(vcov(model)[var_name, var_name])
  
  data.table(
    variable = var_name,
    beta = b,
    IRR_1pp  = exp(b),
    LCL_1pp  = exp(b - 1.96 * se),
    UCL_1pp  = exp(b + 1.96 * se),
    IRR_10pp = exp(b * 10),
    LCL_10pp = exp((b - 1.96 * se) * 10),
    UCL_10pp = exp((b + 1.96 * se) * 10)
  )
}

# ---------------------------------------------------------
# 2) Read data
# ---------------------------------------------------------
dds       <- fread("/home/Ramos/Documentos/R/TUBERCULOSE/IBGE/mundo_onu_adh_municipio.csv", sep = ";")
tb        <- read_fst("/home/Ramos/Documentos/R/TUBERCULOSE/DADOS TRATADOS SINAN/tuberculose_sinan_2003_2024.fst", as.data.table = TRUE)
aps       <- fread("/home/Ramos/Documentos/R/TUBERCULOSE/DADOS TRATADOS SINAN/aps_equipes_2008_2024_anual_media.csv")
pop0_4    <- fread("/home/Ramos/Documentos/R/TUBERCULOSE/IBGE/pop0a4.csv")
pop5_9    <- fread("/home/Ramos/Documentos/R/TUBERCULOSE/IBGE/pop5a9.csv")
pop10_14  <- fread("/home/Ramos/Documentos/R/TUBERCULOSE/IBGE/pop10a14.csv")
pop15_19  <- fread("/home/Ramos/Documentos/R/TUBERCULOSE/IBGE/pop15a19.csv")
pop_total <- fread("/home/Ramos/Documentos/R/TUBERCULOSE/IBGE/pop_total.csv")
bcg       <- fread("/home/Ramos/Documentos/R/TUBERCULOSE/DADOS BRUTOS SINAN/BCG.csv")

# ---------------------------------------------------------
# 3) Fix headers where needed
# ---------------------------------------------------------
pop_total <- maybe_fix_header(pop_total)
pop0_4    <- maybe_fix_header(pop0_4)
pop5_9    <- maybe_fix_header(pop5_9)
pop10_14  <- maybe_fix_header(pop10_14)
pop15_19  <- maybe_fix_header(pop15_19)
bcg       <- maybe_fix_header(bcg)

# ---------------------------------------------------------
# 4) Population tables
# ---------------------------------------------------------
pop_total_long <- wide_to_long_codeyear(pop_total, value_name = "pop_total")
pop0_4_long    <- wide_to_long_codeyear(pop0_4,    value_name = "pop_0_4")
pop5_9_long    <- wide_to_long_codeyear(pop5_9,    value_name = "pop_5_9")
pop10_14_long  <- wide_to_long_codeyear(pop10_14,  value_name = "pop_10_14")
pop15_19_long  <- wide_to_long_codeyear(pop15_19,  value_name = "pop_15_19")

pop_0_19_long <- Reduce(
  function(x, y) merge(x, y, by = c("codigo", "ano"), all = TRUE),
  list(pop0_4_long, pop5_9_long, pop10_14_long, pop15_19_long)
)

pop_0_19_long[, pop_0_19 := rowSums(.SD, na.rm = TRUE),
              .SDcols = c("pop_0_4", "pop_5_9", "pop_10_14", "pop_15_19")]

pop_0_19_long <- pop_0_19_long[, .(
  codigo, ano, pop_0_4, pop_5_9, pop_10_14, pop_15_19, pop_0_19
)]

# ---------------------------------------------------------
# 5) APS rate per 100,000 inhabitants
# ---------------------------------------------------------
setDT(aps)
names(aps) <- clean_names(names(aps))
aps[, codigo := normalize_code6(cod_mun)]

aps_pop <- merge(
  aps[, .(codigo, ano, media_equipes, mediana_equipes, n_meses_validos)],
  pop_total_long,
  by = c("codigo", "ano"),
  all.x = TRUE
)

aps_pop[, aps_equipes_100k := (media_equipes / pop_total) * 100000]

mediana_aps_100k <- aps_pop[, median(aps_equipes_100k, na.rm = TRUE)]
mediana_aps_por_ano <- aps_pop[, .(
  mediana_aps_100k = median(aps_equipes_100k, na.rm = TRUE)
), by = ano][order(ano)]

print(mediana_aps_100k)
print(mediana_aps_por_ano)

# ---------------------------------------------------------
# 6) IDHM 2010
# ---------------------------------------------------------
setDT(dds)
names(dds) <- clean_names(names(dds))

dds_idhm <- dds[ano == 2010, .(
  codigo = substr(as.character(id_municipio), 1, 6),
  idhm
)]

print(summary(dds_idhm$idhm))

# ---------------------------------------------------------
# 7) BCG coverage
# ---------------------------------------------------------
setDT(bcg)
names(bcg) <- clean_names(names(bcg))
bcg[, codigo := extract_code_from_label(bcg[[names(bcg)[1]]])]

bcg_year_cols <- grep("^X\\d{4}$", names(bcg), value = TRUE)
bcg[, (bcg_year_cols) := lapply(.SD, to_num2), .SDcols = bcg_year_cols]

bcg_long <- melt(
  bcg,
  id.vars = "codigo",
  measure.vars = bcg_year_cols,
  variable.name = "ano",
  value.name = "bcg_raw"
)

bcg_long[, ano := as.integer(sub("^X", "", ano))]
bcg_long[, bcg_cap100 := pmin(bcg_raw, 100)]
bcg_long[, bcg_le120  := fifelse(bcg_raw <= 120, bcg_raw, NA_real_)]

# ---------------------------------------------------------
# 8) Prepare TB base
# ---------------------------------------------------------
setDT(tb)

# convert key fields from factor/character safely
tb[, DT_DIAG_chr     := as.character(DT_DIAG)]
tb[, DT_DIAG2        := as.IDate(DT_DIAG_chr)]
tb[, ano_base        := to_num(ano)]
tb[, id_mn_resi_num  := to_num(ID_MN_RESI)]
tb[, id_municip_num  := to_num(ID_MUNICIP)]
tb[, NU_IDADE_N_num  := to_num(NU_IDADE_N)]
tb[, FORMA_num       := to_num(FORMA)]
tb[, HIV_num         := to_num(HIV)]
tb[, CS_RACA_num     := to_num(CS_RACA)]
tb[, CS_GESTANT_num  := to_num(CS_GESTANT)]
tb[, AGRAVALCOO_num  := to_num(AGRAVALCOO)]
tb[, AGRAVDROGA_num  := to_num(AGRAVDROGA)]
tb[, AGRAVTABAC_num  := to_num(AGRAVTABAC)]
tb[, AGRAVDIABE_num  := to_num(AGRAVDIABE)]

# 8.1 Year of diagnosis
tb[, ano_diag := fifelse(
  !is.na(DT_DIAG2),
  as.integer(format(DT_DIAG2, "%Y")),
  ano_base
)]

# 8.2 Municipality code: residence first, notification fallback
tb[, codigo := fifelse(
  !is.na(id_mn_resi_num) & id_mn_resi_num >= 100000 & id_mn_resi_num <= 999999,
  sprintf("%06d", id_mn_resi_num),
  fifelse(
    !is.na(id_municip_num) & id_municip_num >= 100000 & id_municip_num <= 999999,
    sprintf("%06d", id_municip_num),
    NA_character_
  )
)]

tb[codigo %in% c("000000", "", "      "), codigo := NA_character_]

# 8.3 Age groups from NU_IDADE_N
tb[, idade_unidade := NU_IDADE_N_num %/% 1000L]
tb[, idade_valor   := NU_IDADE_N_num %% 1000L]

tb[, faixa_etaria := fcase(
  idade_unidade %in% c(1L, 2L, 3L), "0-4",
  idade_unidade == 4L & idade_valor >= 0L  & idade_valor <= 4L,  "0-4",
  idade_unidade == 4L & idade_valor >= 5L  & idade_valor <= 9L,  "5-9",
  idade_unidade == 4L & idade_valor >= 10L & idade_valor <= 14L, "10-14",
  idade_unidade == 4L & idade_valor >= 15L & idade_valor <= 19L, "15-19",
  default = NA_character_
)]

tb[, faixa_etaria := factor(
  faixa_etaria,
  levels = c("0-4", "5-9", "10-14", "15-19")
)]

# 8.4 Clinical form
tb[, forma_clinica := fcase(
  FORMA_num == 1L, "Pulmonary",
  FORMA_num == 2L, "Extrapulmonary",
  FORMA_num == 3L, "Pulmonary + extrapulmonary",
  default = NA_character_
)]

tb[, forma_clinica := factor(
  forma_clinica,
  levels = c("Pulmonary", "Extrapulmonary", "Pulmonary + extrapulmonary")
)]

tb[, eptb_any := fifelse(
  FORMA_num %in% c(2L, 3L), 1L,
  fifelse(FORMA_num == 1L, 0L, NA_integer_)
)]

# 8.5 Restrict to valid 0-19
tb_0_19 <- tb[
  !is.na(codigo) &
    !is.na(ano_diag) &
    !is.na(faixa_etaria)
]

# 8.6 Checks
cat("\nN tb_0_19:\n")
print(nrow(tb_0_19))

cat("\nDistribuição de FORMA_num em tb_0_19:\n")
print(tb_0_19[, .N, by = FORMA_num][order(FORMA_num)])

cat("\nDistribuição de eptb_any em tb_0_19:\n")
print(tb_0_19[, .N, by = eptb_any][order(eptb_any)])

cat("\nComprimento do código:\n")
print(tb_0_19[, .N, by = .(nchar_codigo = nchar(codigo))][order(nchar_codigo)])

if (tb_0_19[eptb_any == 1L, .N] == 0) {
  stop("tb_0_19 não contém casos com eptb_any == 1. Verifique FORMA/FORMA_num.")
}

# ---------------------------------------------------------
# 9) TB descriptive aggregation
# ---------------------------------------------------------
tb_mun_ano_fx_forma <- tb_0_19[!is.na(forma_clinica), .(
  n_tb = .N
), by = .(codigo, ano = ano_diag, faixa_etaria, forma_clinica)]

setorder(tb_mun_ano_fx_forma, codigo, ano, faixa_etaria, forma_clinica)

tb_mun_ano_forma_total <- tb_0_19[!is.na(forma_clinica), .(
  n_tb = .N
), by = .(codigo, ano = ano_diag, forma_clinica)]

tb_mun_ano_forma_total[, faixa_etaria := "TOTAL_0_19"]
setcolorder(tb_mun_ano_forma_total, c("codigo", "ano", "faixa_etaria", "forma_clinica", "n_tb"))

tb_mun_ano_fx_forma_all <- rbind(
  tb_mun_ano_fx_forma,
  tb_mun_ano_forma_total,
  fill = TRUE
)

# ---------------------------------------------------------
# 10) Ecological outcome: municipality-year count of EPTB
# ---------------------------------------------------------
eptb_mun_ano <- tb_0_19[eptb_any == 1L, .(
  n_eptb = .N
), by = .(codigo, ano = ano_diag)]

cat("\nResumo de eptb_mun_ano:\n")
print(eptb_mun_ano[1:10])
print(eptb_mun_ano[, .N])
print(eptb_mun_ano[, .(
  n_codigos = uniqueN(codigo),
  min_ano = min(ano),
  max_ano = max(ano)
)])

if (nrow(eptb_mun_ano) == 0) {
  stop("eptb_mun_ano ficou vazio.")
}

# ---------------------------------------------------------
# 11) HIV proxy at municipality-year level
# ---------------------------------------------------------
tb[, hiv_pos := fifelse(HIV_num == 1L, 1L,
                        fifelse(HIV_num == 2L, 0L, NA_integer_))]

hiv_proxy <- tb[
  !is.na(codigo) & !is.na(ano_diag),
  .(
    hiv_n_info = sum(!is.na(hiv_pos)),
    tb_hiv_prop = ifelse(sum(!is.na(hiv_pos)) >= 5, mean(hiv_pos, na.rm = TRUE), NA_real_)
  ),
  by = .(codigo, ano = ano_diag)
]

# ---------------------------------------------------------
# 12) Build municipality-year panel with true zeros
# ---------------------------------------------------------
anos_modelo <- Reduce(intersect, list(
  unique(pop_0_19_long$ano),
  unique(aps_pop$ano),
  unique(bcg_long$ano)
))
anos_modelo <- sort(anos_modelo)

painel <- CJ(
  codigo = unique(pop_0_19_long$codigo),
  ano = anos_modelo
)

painel <- merge(painel, pop_0_19_long, by = c("codigo", "ano"), all.x = TRUE)
painel <- merge(painel, eptb_mun_ano,   by = c("codigo", "ano"), all.x = TRUE)
painel[is.na(n_eptb) & !is.na(pop_0_19), n_eptb := 0L]

painel <- merge(
  painel,
  aps_pop[, .(codigo, ano, media_equipes, mediana_equipes, aps_equipes_100k, n_meses_validos)],
  by = c("codigo", "ano"),
  all.x = TRUE
)

painel <- merge(
  painel,
  bcg_long[, .(codigo, ano, bcg_raw, bcg_cap100, bcg_le120)],
  by = c("codigo", "ano"),
  all.x = TRUE
)

painel <- merge(painel, dds_idhm, by = "codigo", all.x = TRUE)
painel <- merge(painel, hiv_proxy, by = c("codigo", "ano"), all.x = TRUE)

# ---------------------------------------------------------
# 13) Final ecological datasets
# ---------------------------------------------------------
painel_main <- painel[
  !is.na(pop_0_19) & pop_0_19 > 0 &
    !is.na(aps_equipes_100k) &
    !is.na(idhm) &
    !is.na(bcg_le120)
]

painel_hiv <- painel_main[!is.na(tb_hiv_prop)]

print(range(painel_main$ano, na.rm = TRUE))
print(painel_main[, .(
  n = .N,
  zeros = sum(n_eptb == 0, na.rm = TRUE),
  nonzero = sum(n_eptb > 0, na.rm = TRUE),
  prop_zero = mean(n_eptb == 0, na.rm = TRUE),
  miss_hiv = sum(is.na(tb_hiv_prop))
)])

if (painel_main[, sum(n_eptb > 0, na.rm = TRUE)] == 0) {
  stop("painel_main ficou com 100% de zeros no desfecho.")
}

# ---------------------------------------------------------
# 14) Negative binomial ecological regressions
# ---------------------------------------------------------
mod_nb_orig <- glm.nb(
  n_eptb ~ bcg_le120 + factor(ano) + offset(log(pop_0_19)),
  data = painel_main
)

mod_nb_main <- glm.nb(
  n_eptb ~ bcg_le120 + factor(ano) + aps_equipes_100k + idhm +
    offset(log(pop_0_19)),
  data = painel_main
)

mod_nb_hiv <- glm.nb(
  n_eptb ~ bcg_le120 + factor(ano) + aps_equipes_100k + idhm + tb_hiv_prop +
    offset(log(pop_0_19)),
  data = painel_hiv
)

summary(mod_nb_orig)
summary(mod_nb_main)
summary(mod_nb_hiv)

# ---------------------------------------------------------
# 15) Practical significance of BCG effect
# ---------------------------------------------------------
irr_nb_orig <- irr_table_nb(mod_nb_orig, "bcg_le120")
irr_nb_main <- irr_table_nb(mod_nb_main, "bcg_le120")
irr_nb_hiv  <- irr_table_nb(mod_nb_hiv,  "bcg_le120")

print(irr_nb_orig)
print(irr_nb_main)
print(irr_nb_hiv)

bcg_summary_main <- painel_main[, .(
  n = .N,
  min = min(bcg_le120, na.rm = TRUE),
  p25 = quantile(bcg_le120, 0.25, na.rm = TRUE),
  median = median(bcg_le120, na.rm = TRUE),
  p75 = quantile(bcg_le120, 0.75, na.rm = TRUE),
  max = max(bcg_le120, na.rm = TRUE),
  iqr = IQR(bcg_le120, na.rm = TRUE)
)]

print(bcg_summary_main)

ref_year <- 2019

newdat_main <- data.table(
  bcg_le120 = c(70, 95),
  ano = ref_year,
  aps_equipes_100k = median(painel_main$aps_equipes_100k, na.rm = TRUE),
  idhm = median(painel_main$idhm, na.rm = TRUE),
  pop_0_19 = 100000
)

newdat_main[, pred_cases := predict(mod_nb_main, newdata = newdat_main, type = "response")]
newdat_main[, pred_rate_100k := pred_cases / pop_0_19 * 100000]

print(newdat_main)

# ---------------------------------------------------------
# 16) Individual-level model
# ---------------------------------------------------------
base_indiv <- tb_0_19[FORMA_num %in% c(1L, 2L, 3L)]

base_indiv[, sexo := fcase(
  CS_SEXO == "M", "Male",
  CS_SEXO == "F", "Female",
  default = NA_character_
)]

base_indiv[, raca5 := fcase(
  CS_RACA_num == 1L, "White",
  CS_RACA_num == 2L, "Black",
  CS_RACA_num == 3L, "Asian",
  CS_RACA_num == 4L, "Brown",
  CS_RACA_num == 5L, "Indigenous",
  default = NA_character_
)]

base_indiv[, raca5 := factor(
  raca5,
  levels = c("White", "Black", "Brown", "Asian", "Indigenous")
)]

base_indiv[, hiv_pos_ind  := fifelse(HIV_num == 1L, 1L, fifelse(HIV_num == 2L, 0L, NA_integer_))]
base_indiv[, alcool       := recode_yes_no(AGRAVALCOO_num)]
base_indiv[, drogas       := recode_yes_no(AGRAVDROGA_num)]
base_indiv[, tabaco       := recode_yes_no(AGRAVTABAC_num)]
base_indiv[, diabetes     := recode_yes_no(AGRAVDIABE_num)]

base_indiv[, gestante := fifelse(CS_GESTANT_num %in% c(1L, 2L, 3L, 4L), 1L,
                                 fifelse(CS_GESTANT_num == 5L, 0L, NA_integer_))]

vars_logit <- c("eptb_any", "faixa_etaria", "sexo", "raca5", "hiv_pos_ind",
                "alcool", "drogas", "tabaco", "diabetes")

base_indiv_main <- base_indiv[complete.cases(base_indiv[, ..vars_logit])]

mod_logit_main <- glm(
  eptb_any ~ faixa_etaria + sexo + raca5 + hiv_pos_ind +
    alcool + drogas + tabaco + diabetes,
  family = binomial(),
  data = base_indiv_main
)

summary(mod_logit_main)

base_indiv_nopreg <- base_indiv_main[!(sexo == "Female" & gestante == 1L)]

mod_logit_nopreg <- glm(
  eptb_any ~ faixa_etaria + sexo + raca5 + hiv_pos_ind +
    alcool + drogas + tabaco + diabetes,
  family = binomial(),
  data = base_indiv_nopreg
)

summary(mod_logit_nopreg)

or_logit_main   <- or_table(mod_logit_main)
or_logit_nopreg <- or_table(mod_logit_nopreg)

print(or_logit_main)
print(or_logit_nopreg)

# ---------------------------------------------------------
# 17) Optional subtype-specific ecological models
# Still depends on EXTRAPU codebook
# ---------------------------------------------------------
# severe_codes     <- c(...)
# pleural_ln_codes <- c(...)

# ---------------------------------------------------------
# 18) Save outputs
# ---------------------------------------------------------
dir_out <- "/home/Ramos/Documentos/R/TUBERCULOSE/RESULTADOS_REVISAO_BCG_EPTB/"
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

fwrite(tb_mun_ano_fx_forma,         file.path(dir_out, "tb_mun_ano_faixa_forma.csv"))
fwrite(tb_mun_ano_forma_total,      file.path(dir_out, "tb_mun_ano_total_forma.csv"))
fwrite(aps_pop,                     file.path(dir_out, "aps_pop.csv"))
fwrite(dds_idhm,                    file.path(dir_out, "idhm_2010.csv"))
fwrite(bcg_long,                    file.path(dir_out, "bcg_long.csv"))
fwrite(painel_main,                 file.path(dir_out, "painel_ecologico_main.csv"))
fwrite(irr_nb_orig,                 file.path(dir_out, "irr_nb_original.csv"))
fwrite(irr_nb_main,                 file.path(dir_out, "irr_nb_main.csv"))
fwrite(irr_nb_hiv,                  file.path(dir_out, "irr_nb_hiv.csv"))
fwrite(newdat_main,                 file.path(dir_out, "predicoes_bcg_70_95.csv"))
fwrite(or_logit_main,               file.path(dir_out, "or_logit_main.csv"))
fwrite(or_logit_nopreg,             file.path(dir_out, "or_logit_nopreg.csv"))

saveRDS(mod_nb_orig,      file.path(dir_out, "mod_nb_original.rds"))
saveRDS(mod_nb_main,      file.path(dir_out, "mod_nb_main.rds"))
saveRDS(mod_nb_hiv,       file.path(dir_out, "mod_nb_hiv.rds"))
saveRDS(mod_logit_main,   file.path(dir_out, "mod_logit_main.rds"))
saveRDS(mod_logit_nopreg, file.path(dir_out, "mod_logit_nopreg.rds"))

# ---------------------------------------------------------
# 19) Key outputs
# ---------------------------------------------------------
cat("\nAPS median per 100,000 inhabitants:\n")
print(mediana_aps_100k)

cat("\nBCG coverage summary in analytical dataset:\n")
print(bcg_summary_main)

cat("\nBCG effect in main negative binomial model:\n")
print(irr_nb_main)

cat("\nPredicted EPTB notification rates at 70% vs 95% BCG coverage:\n")
print(newdat_main)

cat("\nMain individual OR table:\n")
print(or_logit_main)

cat("\nOutputs saved to:\n", dir_out, "\n")
