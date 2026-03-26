# Carregar bibliotecas
library(Seurat)
library(ggplot2)

# 1. Carregar o seu objeto 
seurat_ob <- readRDS("/home/bpdrmorais/sead_project/seurat_mtg.rds")

# 2. Criar a coluna 'Stage' (Late vs Early Braak)
seurat_ob$Stage <- ifelse(seurat_ob$Braak %in% c("Braak IV", "Braak V", "Braak VI"), 
                           "Late Braak", 
                           "Early Braak")
                          
                          
seurat_ob$celltype_simplified <- plyr::revalue(seurat_ob$Subclass, c(
  "Astrocyte"       = "Astro",
  "Chandelier"      = "Inh",
  "Endothelial"     = "Endo",
  "L4 IT"           = "Exc",
  "L5 ET"           = "Exc",
  "L5 IT"           = "Exc",
  "L5/6 NP"         = "Exc",
  "L6 CT"           = "Exc",
  "L6 IT Car3"      = "Exc",
  "L6 IT"           = "Exc",
  "L6b"             = "Exc",
  "Lamp5 Lhx6"      = "Inh",
  "Lamp5"           = "Inh",
  "Microglia-PVM"   = "Micro",
  "Oligodendrocyte" = "Oligo",
  "OPC"             = "OPC",
  "Pax6"            = "Inh",
  "Pvalb"           = "Inh",
  "Sncg"            = "Inh",
  "Sst"             = "Inh",
  "Vip"             = "Inh",
  "VLMC"            = "Endo",
  "Sst Chodl"       = "Inh"
))


# Coluna braal stage + cell type 
seurat_ob$cell_braak <- paste(seurat_ob$Stage, seurat_ob$celltype_simplified, sep = "_")

# 1. Criar a tabela cruzando Tipo Celular (linhas) e Estágio (colunas)
tabela <- table(seurat_ob$celltype_simplified, seurat_ob$Stage)

# 2. Mostrar a tabela na tela do R para você dar uma olhada rápida
print(tabela)




# 4. Plotar o UMAP colorido por essa nova coluna
seurat_ob <- NormalizeData(seurat_ob)

seurat_ob <- FindVariableFeatures(seurat_ob, selection.method = "vst", nfeatures = 2000)

seurat_ob <- ScaleData(seurat_ob)

seurat_ob <- RunPCA(seurat_ob)

# Usamos as 30 primeiras dimensões do PCA, que é o padrão ouro
seurat_ob <- RunUMAP(seurat_ob, dims = 1:30)

meu_umap <- DimPlot(seurat_ob, reduction = "umap", group.by = "cell_braak") + 
  ggtitle("MTG - Stage and Celltype")

plot(meu_umap)


# 6. Salvar o obeto FINAL e atualizado para o Nextflow
saveRDS(seurat_ob, "/home/bpdrmorais/sead_project/mtg_scrank.rds")
