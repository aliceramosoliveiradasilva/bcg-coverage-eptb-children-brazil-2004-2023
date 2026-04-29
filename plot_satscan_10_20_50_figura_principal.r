
# ============================================================
# SaTScan 10%, 20% e 50% - figura limpa para artigo
# Mostra, em cada análise:
#   - os 3 clusters HIGH com maior LLR
#   - o 1 cluster LOW com maior LLR
# Círculos:
#   - vinho escuro = high
#   - azul escuro  = low
# ============================================================

library(sf)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(tidyr)
library(purrr)

# -------------------------
# 1) caminhos
# -------------------------
dir10 <- "/home/ramos/Documentos/R/TUBERCULOSE/SATSCAN/resultado 10"
dir20 <- "/home/ramos/Documentos/R/TUBERCULOSE/SATSCAN/resultado 20/"
dir50 <- "/home/ramos/Documentos/R/TUBERCULOSE/SATSCAN/resultado de 50/"

# AJUSTE para seu shapefile base
arq_brasil <- "/home/ramos/Documentos/R/TUBERCULOSE/IBGE/BR_Municipios_2025/BR_Municipios_2025.shp"

# -------------------------
# 2) função para ler col.txt
# -------------------------
ler_col_txt <- function(arquivo, analise){
  suppressMessages({
    df <- readr::read_table2(
      file = arquivo,
      col_names = c(
        "cluster","loc_id","lat","lon","radius_km","span_km",
        "start_date","end_date","n_locations","llr","p_value",
        "observed","expected","ode","rr","population","gini"
      ),
      locale = readr::locale(decimal_mark = ",", grouping_mark = ".")
    )
  })
  
  df %>%
    dplyr::mutate(
      analysis = analise,
      cluster_type = dplyr::if_else(rr > 1, "High", "Low"),
      period = paste(start_date, "to", end_date)
    )
}

# -------------------------
# 3) ler tabelas
# -------------------------
tab10 <- ler_col_txt(file.path(dir10, "Resultados 10.col.txt"), "10%")
tab20 <- ler_col_txt(file.path(dir20, "Resultados 20.col.txt"), "20%")
tab50 <- ler_col_txt(file.path(dir50, "Resultados 50.col.txt"), "50%")

tab_all <- bind_rows(tab10, tab20, tab50)

# -------------------------
# 4) selecionar clusters para a figura
#    3 HIGH + 1 LOW por análise
# -------------------------
top_high <- tab_all %>%
  filter(cluster_type == "High") %>%
  group_by(analysis) %>%
  slice_max(order_by = llr, n = 3, with_ties = FALSE) %>%
  ungroup()

top_low <- tab_all %>%
  filter(cluster_type == "Low") %>%
  group_by(analysis) %>%
  slice_max(order_by = llr, n = 1, with_ties = FALSE) %>%
  ungroup()

clusters_keep <- bind_rows(top_high, top_low) %>%
  arrange(analysis, cluster_type, desc(llr))

# exportar tabela resumida da figura
clusters_keep %>%
  select(
    analysis, cluster, cluster_type, period, n_locations,
    radius_km, observed, expected, rr, llr, p_value
  ) %>%
  write_csv("clusters_figura_principal_10_20_50.csv")

# -------------------------
# 5) ler shapefiles dos círculos
# -------------------------
shp10 <- st_read(file.path(dir10, "Resultados 10.col.shp"), quiet = TRUE) %>%
  mutate(cluster = as.integer(CLUSTER), analysis = "10%")

shp20 <- st_read(file.path(dir20, "Resultados 20.col.shp"), quiet = TRUE) %>%
  mutate(cluster = as.integer(CLUSTER), analysis = "20%")

shp50 <- st_read(file.path(dir50, "Resultados 50.col.shp"), quiet = TRUE) %>%
  mutate(cluster = as.integer(CLUSTER), analysis = "50%")

clusters_sf <- bind_rows(
  shp10 %>% select(cluster, analysis, geometry),
  shp20 %>% select(cluster, analysis, geometry),
  shp50 %>% select(cluster, analysis, geometry)
) %>%
  left_join(tab_all, by = c("analysis", "cluster")) %>%
  semi_join(clusters_keep %>% select(analysis, cluster), by = c("analysis", "cluster"))

# -------------------------
# 6) base do Brasil
# -------------------------
brasil <- st_read(arq_brasil, quiet = TRUE)

if (st_crs(brasil) != st_crs(clusters_sf)) {
  brasil <- st_transform(brasil, st_crs(clusters_sf))
}

# -------------------------
# 7) rótulo opcional dos clusters
# -------------------------
label_df <- clusters_sf %>%
  st_centroid() %>%
  mutate(
    label = paste0("C", cluster)
  )

# -------------------------
# 8) mapa principal
# -------------------------
p <- ggplot() +
  geom_sf(data = brasil, fill = "grey97", color = "grey82", linewidth = 0.04) +
  geom_sf(
    data = clusters_sf,
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
  scale_color_manual(
    values = c(
      "High" = "#6D071A",  # vinho escuro
      "Low"  = "#0B2E59"   # azul escuro
    )
  ) +
  labs(color = "Cluster type") +
  theme_void() +
  theme(
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("mapa_satscan_figura_principal_10_20_50.png", p, width = 15, height = 6, dpi = 500)
ggsave("mapa_satscan_figura_principal_10_20_50.pdf", p, width = 15, height = 6)

# -------------------------
# 9) versão sem LOW
#    útil se quiser destacar só alto risco
# -------------------------
p_high <- ggplot() +
  geom_sf(data = brasil, fill = "grey97", color = "grey82", linewidth = 0.04) +
  geom_sf(
    data = clusters_sf %>% filter(cluster_type == "High"),
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
    strip.text = element_text(face = "bold", size = 11)
  )

ggsave("mapa_satscan_figura_principal_high_only_10_20_50.png", p_high, width = 15, height = 6, dpi = 500)
ggsave("mapa_satscan_figura_principal_high_only_10_20_50.pdf", p_high, width = 15, height = 6)

# -------------------------
# 10) mensagem final no console
# -------------------------
print("Arquivos gerados:")
print("clusters_figura_principal_10_20_50.csv")
print("mapa_satscan_figura_principal_10_20_50.png")
print("mapa_satscan_figura_principal_10_20_50.pdf")
print("mapa_satscan_figura_principal_high_only_10_20_50.png")
print("mapa_satscan_figura_principal_high_only_10_20_50.pdf")
