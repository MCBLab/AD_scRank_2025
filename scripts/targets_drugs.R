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

# Threshold global — top 2.5%
threshold_global <- quantile(perb$perb_score, 0.975)
perb_top <- perb %>% filter(perb_score >= threshold_global)

# Classificar genes — seletividade + stage_pref
result <- perb_top %>%
  group_by(target) %>%
  summarise(
    n_celltypes   = n_distinct(cell_type_clean),
    top_celltypes = paste(unique(cell_type_clean), collapse = ", "),
    top_stages    = paste(unique(stage), collapse = "/"),
    max_score     = max(perb_score),
    mean_score    = mean(perb_score),
    .groups = "drop"
  ) %>%
  mutate(
    selectivity = ifelse(n_celltypes <= 2, "Seletivo", "Não seletivo"),
    stage_pref  = case_when(
      top_stages == "Early" ~ "Early only",
      top_stages == "Late"  ~ "Late only",
      TRUE                  ~ "Both stages"
    ),
    region = "MTG"
  )

# Cruzar com tabela da Julia
drug_table <- read.csv("sead_project/data/ad_drug_targets.csv")

drug_info <- drug_table %>%
  filter(HGNC_Symbol != "Unknown") %>%
  group_by(HGNC_Symbol) %>%
  summarise(
    drugs                = paste(unique(Drug), collapse = "; "),
    functional_category  = paste(unique(Functional_Category), collapse = "/"),
    target_type          = paste(unique(Target_Type), collapse = "; "),
    max_phase            = min(Phase),
    .groups = "drop"
  ) %>%
  rename(target = HGNC_Symbol)

result <- result %>% left_join(drug_info, by = "target")

result_final <- result %>%
  filter(
    is.na(functional_category) |
      grepl("antagonist-like", functional_category)
  )

# Selecionar só os Seletivos
targets_finais <- result_final %>%
  filter(selectivity == "Seletivo") %>%
  pull(target)

# Matriz de scores
heatmap_df <- perb %>%
  mutate(col_name = paste0(cell_type_clean, "_", stage)) %>%
  select(target, col_name, perb_score) %>%
  filter(target %in% targets_finais) %>%
  pivot_wider(names_from = col_name, values_from = perb_score, values_fill = 0)

heatmap_mat <- as.matrix(heatmap_df[, -1])
rownames(heatmap_mat) <- heatmap_df$target
heatmap_mat_log <- log10(heatmap_mat + 1e-20)

col_order_sorted <- sort(colnames(heatmap_mat_log))
heatmap_mat_log  <- heatmap_mat_log[, col_order_sorted]

# Labels com [E] ou [L]
stage_of_col      <- ifelse(grepl("_Early$", colnames(heatmap_mat_log)), "E", "L")
colnames_clean    <- gsub("_Early$|_Late$", "", colnames(heatmap_mat_log))
col_labels_tagged <- paste0(colnames_clean, " [", stage_of_col, "]")

# Ordenar linhas por stage_pref
row_meta <- result_final %>%
  filter(target %in% rownames(heatmap_mat_log)) %>%
  arrange(factor(stage_pref, levels = c("Early only", "Both stages", "Late only")),
          desc(max_score))

heatmap_mat_log <- heatmap_mat_log[row_meta$target, ]

# Asterisco em QUALQUER célula que passou o threshold 
heatmap_df_bool <- perb %>%
  mutate(
    col_name   = paste0(cell_type_clean, "_", stage),
    passou_top = perb_score >= threshold_global
  ) %>%
  select(target, col_name, passou_top) %>%
  filter(target %in% targets_finais) %>%
  pivot_wider(names_from = col_name, values_from = passou_top, values_fill = FALSE)

bool_mat <- as.matrix(heatmap_df_bool[, -1])
rownames(bool_mat) <- heatmap_df_bool$target
bool_mat <- bool_mat[rownames(heatmap_mat_log), colnames(heatmap_mat_log)]

# Anotação de linha — só Fase
fase_vals <- as.character(row_meta$max_phase)
fase_vals[is.na(fase_vals)] <- "NA"

row_annot <- rowAnnotation(
  Fase = fase_vals,
  col  = list(
    Fase = c("1" = "#FFF9C4", "2" = "#FFB300", "3" = "#E65100", "NA" = "#BDBDBD")
  )
)

# Anotação de colunas — barra Early/Late colorida
stage_annot <- ifelse(grepl("_Early$", colnames(heatmap_mat_log)), "Early", "Late")

col_annot <- HeatmapAnnotation(
  Stage = stage_annot,
  col   = list(Stage = c("Early" = "#2196F3", "Late" = "#F44336"))
)

# Heatmap final
png("heatmap_seletivos_drug.png", width = 1800, height = 1400, res = 130)
draw(Heatmap(
  heatmap_mat_log,
  name             = "log10(score)",
  top_annotation   = col_annot,
  left_annotation  = row_annot,
  row_split        = factor(row_meta$stage_pref,
                            levels = c("Early only", "Both stages", "Late only")),
  row_title_gp     = gpar(fontsize = 11, fontface = "bold"),
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
  column_title = "MTG Antagonista — Seletivos|Drugs"
), heatmap_legend_side = "right")
dev.off()

# Exportar tabela final
write.csv(result_final, "target_classification_MTG_antagonist_final.csv", row.names = FALSE)