library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

#=============================
# 1) caminho do arquivo CSV
#=============================
arq <- "/home/Ramos/Documentos/R/TUBERCULOSE/Joinpoint tb extr.csv"

#=============================
# 2) leitura
#=============================
dados <- read_csv(arq, show_col_types = FALSE)

#=============================
# 3) padronização e ordem das regiões
#=============================
ordem_regioes <- c("Brazil", "North", "North East", "Southeast", "South", "Midwest")

dados <- dados %>%
  mutate(
    Region = trimws(Region),
    Region = recode(
      Region,
      "Northeast" = "North East",
      "NorthEast" = "North East"
    ),
    Region = factor(Region, levels = ordem_regioes),
    `age range` = factor(`age range`, levels = c("0-19 years", "0-4 years"))
  )

# checagem rápida
print(unique(dados$Region))
print(table(dados$Region, useNA = "ifany"))
print(table(dados$`age range`, useNA = "ifany"))

#=============================
# 4) função para cada painel
#=============================
faz_grafico <- function(base, faixa, tag){
  ggplot(
    base %>% filter(`age range` == faixa),
    aes(x = Year)
  ) +
    geom_line(aes(y = Modeled), linewidth = 0.9) +
    geom_point(aes(y = Observed), size = 2) +
    facet_wrap(
      ~ Region,
      ncol = 3,
      nrow = 2,
      scales = "fixed",
      drop = FALSE
    ) +
    labs(
      x = "Year",
      y = "Notification Rate",
      tag = tag
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text = element_text(face = "bold"),
      plot.tag = element_text(face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

#=============================
# 5) painéis A e B
#=============================
pA <- faz_grafico(dados, "0-19 years", "A")
pB <- faz_grafico(dados, "0-4 years",  "B")

#=============================
# 6) figura final
#=============================
fig_final <- pA / pB

print(fig_final)

#=============================
# 7) exportação
#=============================
out_dir <- dirname(arq)

arquivo_png <- file.path(out_dir, "grafico_joinpoint_tb_extr_3x2.png")
arquivo_jpg <- file.path(out_dir, "grafico_joinpoint_tb_extr_3x2.jpg")

ggsave(
  filename = arquivo_png,
  plot = fig_final,
  width = 12,
  height = 14,
  dpi = 300
)

ggsave(
  filename = arquivo_jpg,
  plot = fig_final,
  width = 12,
  height = 14,
  dpi = 300
)

cat("Arquivos salvos em:\n")
cat(arquivo_png, "\n")
cat(arquivo_jpg, "\n")