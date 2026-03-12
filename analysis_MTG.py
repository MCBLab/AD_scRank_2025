import scanpy as sc
import pandas as pd
import numpy as np 

print("1. Bibliotecas carregadas")

# Carregar arquivo H5AD
print("2. Abrindo dataset MTG ...")
adata = sc.read_h5ad("/home/bpdrmorais/sead_project/SEAAD_MTG_RNAseq_final-nuclei.2024-02-13.h5ad", backed="r")
print("Dataset aberto com sucesso!")
print(adata)

# Carregar as IDs cells
print("3. Carregando lista de células MTG...")
cells = pd.read_csv("/home/bpdrmorais/sead_project/selected_cells_mtg.csv")
cell_ids = cells.iloc[:,0].tolist()
print(f"Número de células na lista: {len(cell_ids)}")

# Diagnóstico visual
print("--- DIAGNÓSTICO DE BARCODES ---")
print("5 primeiros IDs do H5AD:", adata.obs_names[:5].tolist())
print("5 primeiros IDs do CSV:", cell_ids[:5])
print("-------------------------------")

# Correspondecia das celulas
print("4. Encontrando células correspondentes e criando subset...")
celulas_validas = adata.obs_names.isin(cell_ids)
subset = adata[celulas_validas].to_memory()

print(f"Células encontradas e mantidas no dataset: {subset.n_obs}")

# Salvando o resultado 
print("5. Salvando o arquivo final...")
salvando = "/home/bpdrmorais/sead_project/subset_final_mtg.h5ad"
subset.write_h5ad(salvando)
print(f"Subset criado e salvo com sucesso em: {salvando}")