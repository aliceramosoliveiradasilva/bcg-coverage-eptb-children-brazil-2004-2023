
# ============================================================
# SaTScan SPACE-TIME PERMUTATION: 10%, 20% e 50% lado a lado
# - detecta automaticamente .col.shp e .col.txt em cada pasta
# - usa os círculos dos clusters (col.shp)
# - classifica high/low a partir de Observed/Expected do col.txt
# - plota 3 painéis: 10%, 20%, 50%
# ============================================================

library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)

# -------------------------
# 1) CAMINHOS
# -------------------------
dir10 <- "/home/ramos/Documentos/R/TUBERCULOSE/SATSCAN/resultado da permutacao 10"
dir20 <- "/home/ramos/Documentos/R/TUBERCULOSE/SATSCAN/resultado permutacao 20"
dir50 <- "/home/ramos/Documentos/R/TUBERCULOSE/SATSCAN/resultado permutacao 50"

# shapefile base do Brasil
arq_brasil <- "/home/ramos/Documentos/R/TUBERCULOSE/IBGE/BR_Municipios_2025/BR_Municipios_2025.shp"

# -------------------------
# 2) FUNÇÕES AUXILIARES
# -------------------------
achar_arquivo <- function(dir, padrao) {
  arqs <- list.files(dir, pattern = padrao, full.names = TRUE, ignore.case = TRUE)
  if (length(arqs) == 0) stop(paste("Arquivo não encontrado em:", dir, "padrão:", padrao))
  arqs[1]
}

ler_col_txt_permutacao <- function(arquivo, analise) {
  suppressMessages({
    df <- readr::read_table2(
      file = arquivo,
      col_names = c(
        "cluster","loc_id","lat","lon","radius_km","span_km",
        "start_date","end_date","n_locations","test_stat","p_value",
        "observed","expected","ode","gini"
      ),
      locale = readr::locale(decimal_mark = ",", grouping_mark = ".")
    )
  })

  df %>%
    mutate(
      analysis = analise,
      cluster_type = if_else(ode > 1, "High", "Low"),
      period = paste(start_date, "to", end_date)
    )
}

ler_shp_col <- function(dir, analise) {
  shp <- achar_arquivo(dir, "\\.col\\.shp$")
  x <- st_read(shp, quiet = TRUE)

  # tenta padronizar o nome da coluna cluster
  nm <- names(x)
  id_col <- nm[str_detect(tolower(nm), "^cluster$|clust")]
  if (length(id_col) == 0) stop(paste("Não encontrei coluna de cluster em", shp))

  x %>%
    mutate(
      cluster = as.integer(.data[[id_col[1]]]),
      analysis = analise
    ) %>%
    select(cluster, analysis, geometry)
}

# -------------------------
# 3) LER TABELAS E SHAPEFILES
# -------------------------
col10 <- achar_arquivo(dir10, "\\.col\\.txt$")
col20 <- achar_arquivo(dir20, "\\.col\\.txt$")
col50 <- achar_arquivo(dir50, "\\.col\\.txt$")

tab10 <- ler_col_txt_permutacao(col10, "10%")
tab20 <- ler_col_txt_permutacao(col20, "20%")
tab50 <- ler_col_txt_permutacao(col50, "50%")

tab_all <- bind_rows(tab10, tab20, tab50)

shp10 <- ler_shp_col(dir10, "10%")
shp20 <- ler_shp_col(dir20, "20%")
shp50 <- ler_shp_col(dir50, "50%")

clusters_sf <- bind_rows(shp10, shp20, shp50) %>%
  left_join(tab_all, by = c("analysis", "cluster"))

# -------------------------
# 4) SELEÇÃO DOS CLUSTERS
#    opções:
#    A) todos os clusters significativos
#    B) top 3 HIGH + top 1 LOW por análise
# -------------------------

# --- opção B (mais limpa para artigo)
top_high <- tab_all %>%
  filter(cluster_type == "High", p_value <= 0.05) %>%
  group_by(analysis) %>%
  slice_max(order_by = test_stat, n = 3, with_ties = FALSE) %>%
  ungroup()

top_low <- tab_all %>%
  filter(cluster_type == "Low", p_value <= 0.05) %>%
  group_by(analysis) %>%
  slice_max(order_by = test_stat, n = 1, with_ties = FALSE) %>%
  ungroup()

clusters_keep <- bind_rows(top_high, top_low)

# se quiser TODOS, troque a linha abaixo por:
# clusters_plot <- clusters_sf %>% filter(p_value <= 0.05)
clusters_plot <- clusters_sf %>%
  semi_join(clusters_keep %>% select(analysis, cluster), by = c("analysis", "cluster"))

# -------------------------
# 5) MAPA BASE
# -------------------------
brasil <- st_read(arq_brasil, quiet = TRUE)

if (st_crs(brasil) != st_crs(clusters_plot)) {
  brasil <- st_transform(brasil, st_crs(clusters_plot))
}

label_df <- clusters_plot %>%
  st_centroid() %>%
  mutate(label = paste0("C", cluster))

# -------------------------
# 6) FIGURA PRINCIPAL
# -------------------------
p <- ggplot() +
  geom_sf(data = brasil, fill = "grey97", color = "grey82", linewidth = 0.04) +

  # preenchimento transparente
  geom_sf(
    data = clusters_plot,
    aes(fill = cluster_type),
    color = NA,
    alpha = 0.22
  ) +

  # contorno
  geom_sf(
    data = clusters_plot,
    aes(color = cluster_type),
    fill = NA,
    linewidth = 0.9
  ) +

  geom_sf_text(
    data = label_df,
    aes(label = label, color = cluster_type),
    size = 3,
    fontface = "bold",
    check_overlap = TRUE
  ) +

  facet_wrap(~analysis, ncol = 3) +

  scale_fill_manual(
    values = c(
      "High" = "#6D071A",  # vinho escuro
      "Low"  = "#0B2E59"   # azul escuro
    )
  ) +
  scale_color_manual(
    values = c(
      "High" = "#6D071A",
      "Low"  = "#0B2E59"
    )
  ) +

  labs(fill = "Cluster type", color = "Cluster type") +
  theme_void() +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("mapa_permutacao_10_20_50.png", p, width = 15, height = 6, dpi = 500, bg = "white")
ggsave("mapa_permutacao_10_20_50.pdf", p, width = 15, height = 6, bg = "white")

# -------------------------
# 7) VERSÃO SOMENTE CLUSTERS HIGH
# -------------------------
p_high <- ggplot() +
  geom_sf(data = brasil, fill = "grey60", color = "grey82", linewidth = 0.04) +

  geom_sf(
    data = clusters_plot %>% filter(cluster_type == "High"),
    fill = "#6D071A",
    color = NA,
    alpha = 0.22
  ) +

  geom_sf(
    data = clusters_plot %>% filter(cluster_type == "High"),
    color = "#6D071A",
    fill = NA,
    linewidth = 0.9
  ) +

  geom_sf_text(
    data = label_df %>% filter(cluster_type == "High"),
    aes(label = label),
    color = "#6D071A",
    size = 3,
    fontface = "bold",
    check_overlap = TRUE
  ) +

  facet_wrap(~analysis, ncol = 3) +
  theme_void() +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("mapa_permutacao_10_20_50_high_only.png", p_high, width = 15, height = 6, dpi = 500, bg = "white")
ggsave("mapa_permutacao_10_20_50_high_only.pdf", p_high, width = 15, height = 6, bg = "white")

# -------------------------
# 8) TABELA-RESUMO
# -------------------------
clusters_keep %>%
  select(
    analysis, cluster, cluster_type, period, n_locations,
    radius_km, observed, expected, ode, test_stat, p_value
  ) %>%
  arrange(analysis, cluster_type, desc(test_stat)) %>%
  write_csv("clusters_permutacao_10_20_50_resumo.csv")

print("Arquivos gerados:")
print("mapa_permutacao_10_20_50.png")
print("mapa_permutacao_10_20_50.pdf")
print("mapa_permutacao_10_20_50_high_only.png")
print("mapa_permutacao_10_20_50_high_only.pdf")
print("clusters_permutacao_10_20_50_resumo.csv")
