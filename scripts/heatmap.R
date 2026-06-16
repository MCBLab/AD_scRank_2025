library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# ============================================================
# INPUT
# ============================================================
perb <- read.table("nf_scrank/results/rank_scores/perbscore_all_targets_mtg.txt",
                   header = TRUE, sep = "\t")

# ============================================================
# PASSO 1: Filtrar só antagonista
# ============================================================
perb <- perb %>% filter(binding == "antagonist")

# ============================================================
# PASSO 2: Extrair stage e cell type limpo
# ============================================================
perb <- perb %>%
  mutate(
    stage           = ifelse(grepl("^Early", cell_type), "Early", "Late"),
    cell_type_clean = gsub("^Early_Braak_|^Late_Braak_", "", cell_type)
  )

# ============================================================
# PASSO 3: Threshold global — top 2.5%
# ============================================================
threshold_global <- quantile(perb$perb_score, 0.975)
cat("Threshold global (top 2.5%):", threshold_global, "\n")

perb_top <- perb %>% filter(perb_score >= threshold_global)

cat("Total de linhas no top 2.5%:", nrow(perb_top), "\n")
cat("Targets únicos:", n_distinct(perb_top$target), "\n")
cat("Cell types únicos:", n_distinct(perb_top$cell_type_clean), "\n")

# ============================================================
# PASSO 4: Classificar targets
# ============================================================
result <- perb_top %>%
  group_by(target) %>%
  summarise(
    n_celltypes   = n_distinct(cell_type_clean),
    top_celltypes = paste(unique(cell_type_clean), collapse = ", "),
    n_stages      = n_distinct(stage),
    top_stages    = paste(unique(stage), collapse = "/"),
    max_score     = max(perb_score),
    mean_score    = mean(perb_score),
    .groups = "drop"
  ) %>%
  mutate(
    selectivity = case_when(
      n_celltypes <= 2 ~ "Seletivo",
      TRUE             ~ "Não seletivo"
    ),
    stage_pref = case_when(
      top_stages == "Early" ~ "Early only",
      top_stages == "Late"  ~ "Late only",
      TRUE                  ~ "Both stages"
    ),
    region = "MTG"
  )

# ============================================================
# PASSO 5: Resumo
# ============================================================
cat("\n=== RESUMO ===\n")
result %>%
  count(selectivity, stage_pref) %>%
  arrange(selectivity, stage_pref) %>%
  print()

cat("\n=== SELETIVOS ===\n")
result %>%
  filter(selectivity == "Seletivo") %>%
  arrange(stage_pref, desc(max_score)) %>%
  select(target, stage_pref, top_celltypes, max_score) %>%
  print(n = Inf)

cat("\n=== NÃO SELETIVOS ===\n")
result %>%
  filter(selectivity == "Não seletivo") %>%
  arrange(stage_pref, desc(max_score)) %>%
  select(target, stage_pref, top_celltypes, n_celltypes, max_score) %>%
  print(n = Inf)

# ============================================================
# PASSO 6: Preparar matriz para heatmap
# ============================================================
heatmap_df <- perb %>%
  mutate(col_name = paste0(cell_type_clean, "_", stage)) %>%
  select(target, col_name, perb_score) %>%
  filter(target %in% result$target) %>%
  pivot_wider(names_from = col_name, values_from = perb_score, values_fill = 0)

heatmap_mat <- as.matrix(heatmap_df[, -1])
rownames(heatmap_mat) <- heatmap_df$target

# Log10
heatmap_mat_log <- log10(heatmap_mat + 1e-20)

# Ordenar colunas: Early e Late de cada cell type lado a lado
col_order_sorted  <- sort(colnames(heatmap_mat_log))
heatmap_mat_log   <- heatmap_mat_log[, col_order_sorted]

# Nomes limpos das colunas
colnames_clean <- gsub("_Early$|_Late$", "", colnames(heatmap_mat_log))

# ============================================================
# PASSO 7: Ordenar linhas por stage_pref
# ============================================================
row_meta <- result %>%
  filter(target %in% rownames(heatmap_mat_log)) %>%
  arrange(factor(stage_pref, levels = c("Early only", "Both stages", "Late only")),
          desc(max_score))

heatmap_mat_log <- heatmap_mat_log[row_meta$target, ]

# ============================================================
# PASSO 8: Anotações
# ============================================================

# Colunas: Early vs Late
stage_annot <- ifelse(grepl("Early", colnames(heatmap_mat_log)), "Early", "Late")
col_annot <- HeatmapAnnotation(
  Stage = stage_annot,
  col   = list(Stage = c("Early" = "#2196F3", "Late" = "#F44336"))
)

# Linhas: Seletividade
row_annot <- rowAnnotation(
  Seletividade = row_meta$selectivity,
  col = list(
    Seletividade = c(
      "Seletivo"     = "#4CAF50",
      "Não seletivo" = "#FF9800"
    )
  )
)

# ============================================================
# PASSO 9: Heatmap dividido por stage_pref
# ============================================================
png("heatmap_MTG_antagonist_final.png", width = 1800, height = 1400, res = 130)
draw(Heatmap(
  heatmap_mat_log,
  name             = "log10(score)",
  top_annotation   = col_annot,
  left_annotation  = row_annot,
  
  # Dividir blocos por stage_pref
  row_split        = factor(row_meta$stage_pref,
                            levels = c("Early only", "Both stages", "Late only")),
  row_title_gp     = gpar(fontsize = 11, fontface = "bold"),
  
  cluster_columns      = FALSE,  # Early/Late de cada cell type lado a lado
  cluster_rows         = TRUE,   # clustering dentro de cada bloco
  cluster_row_slices   = FALSE,  # não mistura os blocos
  
  show_row_dend        = FALSE,
  column_labels        = colnames_clean,
  show_row_names       = TRUE,
  show_column_names    = TRUE,
  column_names_gp      = gpar(fontsize = 8),
  row_names_gp         = gpar(fontsize = 9),
  row_names_max_width  = unit(5, "cm"),
  column_title         = "MTG Antagonista — SELETIVO X NÃO SELETIVO"
))
dev.off()

# ============================================================
# PASSO 10: Exportar tabela
# ============================================================
write.csv(result, "target_classification_MTG_antagonist.csv", row.names = FALSE)
cat("\nArquivos gerados:\n")
cat("- heatmap_MTG_antagonist_final.png\n")
cat("- target_classification_MTG_antagonist.csv\n")