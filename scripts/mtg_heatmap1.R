library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

# Dados

perb <- read.table("AD_scRank_2025/data/perbscore_all_targets_mtg.txt",
                   header = TRUE, sep = "\t") %>%
  filter(binding == "antagonist") %>%
  mutate(
    stage           = ifelse(grepl("^Early", cell_type), "Early", "Late"),
    cell_type_clean = gsub("^Early_Braak_|^Late_Braak_", "", cell_type)
  )

#Threshold global — top 2.5%

threshold_global <- quantile(perb$perb_score, 0.975)
perb_top <- perb %>% filter(perb_score >= threshold_global)

#Classificar genes (seletividade)
result <- perb_top %>%
  group_by(target) %>%
  summarise(
    n_celltypes = n_distinct(cell_type_clean),
    max_score   = max(perb_score),
    .groups = "drop"
  ) %>%
  mutate(selectivity = ifelse(n_celltypes <= 2, "Seletivo", "Não seletivo"))

# Matriz de scores
heatmap_df <- perb %>%
  mutate(col_name = paste0(cell_type_clean, "_", stage)) %>%
  select(target, col_name, perb_score) %>%
  filter(target %in% result$target) %>%
  pivot_wider(names_from = col_name, values_from = perb_score, values_fill = 0)

heatmap_mat <- as.matrix(heatmap_df[, -1])
rownames(heatmap_mat) <- heatmap_df$target
heatmap_mat_log <- log10(heatmap_mat + 1e-20)

col_order_sorted <- sort(colnames(heatmap_mat_log))
heatmap_mat_log  <- heatmap_mat_log[, col_order_sorted]

#Labels
stage_of_col      <- ifelse(grepl("_Early$", colnames(heatmap_mat_log)), "E", "L")
colnames_clean    <- gsub("_Early$|_Late$", "", colnames(heatmap_mat_log))
col_labels_tagged <- paste0(colnames_clean, " [", stage_of_col, "]")

row_meta <- result %>% filter(target %in% rownames(heatmap_mat_log))
heatmap_mat_log <- heatmap_mat_log[row_meta$target, ]

# Asterisco em QUALQUER célula que passou o threshold 

heatmap_df_bool <- perb %>%
  mutate(
    col_name   = paste0(cell_type_clean, "_", stage),
    passou_top = perb_score >= threshold_global
  ) %>%
  select(target, col_name, passou_top) %>%
  filter(target %in% result$target) %>%
  pivot_wider(names_from = col_name, values_from = passou_top, values_fill = FALSE)

bool_mat <- as.matrix(heatmap_df_bool[, -1])
rownames(bool_mat) <- heatmap_df_bool$target
bool_mat <- bool_mat[rownames(heatmap_mat_log), colnames(heatmap_mat_log)]

#Anotação de colunas 
stage_annot <- ifelse(grepl("_Early$", colnames(heatmap_mat_log)), "Early", "Late")

col_annot <- HeatmapAnnotation(
  Stage = stage_annot,
  col   = list(Stage = c("Early" = "#2196F3", "Late" = "#F44336")),
  annotation_legend_param = list(
    Stage = list(title = "Stage")
  )
)

# ============================================================
# PASSO 7: Heatmap
# ============================================================
png("heatmap_seletivoxnseleyivo.png", width = 1800, height = 1400, res = 130)
draw(Heatmap(
  heatmap_mat_log,
  name                 = "log10(score)",
  top_annotation       = col_annot,              # ← barra de Stage de volta
  row_split            = row_meta$selectivity,
  row_title_gp         = gpar(fontsize = 10, fontface = "bold"),
  cluster_columns      = FALSE,
  cluster_rows         = TRUE,
  cluster_row_slices   = FALSE,
  show_row_dend        = FALSE,
  column_labels        = col_labels_tagged,      
  show_row_names       = TRUE,
  show_column_names    = TRUE,
  column_names_gp      = gpar(fontsize = 8),
  row_names_gp         = gpar(fontsize = 9),
  row_names_max_width  = unit(5, "cm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (bool_mat[i, j]) {
      grid.text("*", x, y, gp = gpar(fontsize = 14, fontface = "bold", col = "black"))
    }
  },
  column_title = "MTG Antagonistas - Seletivo x Não Seletivo"
), heatmap_legend_side = "right")
dev.off()