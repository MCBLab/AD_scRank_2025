library(dplyr)

MTG = read.csv("bpdrmorais/sead_project/SEAAD_MTG_RNAseq_all-nuclei_metadata.2024-02-13.csv")

#FILTRAR PACIENTES COM DEMENCIA 
MTG_DEM = MTG %>%
  filter(Cognitive.Status == "Dementia")

MTG_braak = MTG_DEM %>%
  mutate(Stage = ifelse(Braak %in% c("Braak IV","Braak V","Braak VI"),
                        "Late Braak",
                        "Early Braak")) %>%
  group_by(Sex, Stage, Subclass) %>%
  tally()

View(MTG_braak)

#QUANTIDADE DE CELULAS POR TIPO
MTG_braak %>%
  group_by(Subclass) %>%
  summarise(n_cells = n())

#SUBSET CELULAS 
set.seed(123)

subset_mtg <- MTG_DEM %>%
  group_by(Subclass) %>%
  filter(n() >= 200) %>% 
  slice_sample(prop = 1) %>%  
  slice_head(n = 1000) %>%     
  ungroup()

#SAMPLE IDS
selected_cells <- subset_mtg$sample_id

write.csv(selected_cells,
          "selected_cells_mtg.csv",
          row.names = FALSE)
