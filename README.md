# AD_scRank_2025 🧬

**MCBLab/AD_scRank_2025** applies gene perturbation modeling to single-cell 
RNA sequencing (scRNA-seq) data to identify prioritized therapeutic targets 
and characterize cell-type–specific vulnerability patterns throughout the 
progression of Alzheimer's disease (AD).

We use [scRank](https://www.sciencedirect.com/science/article/pii/S266637912400260X), 
a computational tool that simulates gene perturbations in co-expression networks 
and estimates how different cell types respond to changes in disease-associated genes. 
The entire workflow was implemented using [NF_scRank](https://github.com/MCBLab/nf_scrank), 
ensuring scalability, reproducibility, and modularity.

---

## Dataset

This project uses publicly available data from the 
[Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD)](https://sea-ad.org/):

| Dataset | Region | File |
|---------|--------|------|
| SEAAD MTG | Middle Temporal Gyrus | [Download](https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/RNAseq/SEAAD_MTG_RNAseq_all-nuclei.2024-02-13.h5ad) |
| SEAAD DLPFC | Dorsolateral Prefrontal Cortex | [Download](https://sea-ad-single-cell-profiling.s3.amazonaws.com/PFC/RNAseq/SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad) |

---

## Pipeline Summary

### 1. Pre-processing (R/Python/Scanpy)
Filtering patients with Dementia diagnosis, classifying samples into Early 
(Braak I–III) and Late (Braak IV–VI) stages, and selecting 200–1,000 cells 
per cell type. Exports a curated Seurat object (`.rds`) with a `cell_braak` 
metadata column for downstream analysis.

### 2. NF_scRank Pipeline
The curated object is passed to the [NF_scRank](https://github.com/MCBLab/nf_scrank) 
pipeline, which runs GENIE3 network inference and scRank perturbation scoring 
for each cell type and target gene.
```bash
nextflow run main.nf \
  --obj /path/to/SEAD_MTG_braakstage_symbol.rds \
  --column cell_braak \
  --species human \
  --n_cells 1000 \
  --n_cores 16 \
  --target input_bellenguez.txt \
  -profile singularity \
  --outdir results \
  -resume
```

### 3. Downstream Analysis (R)
Perturbation scores are loaded into R for visualization and interpretation:
- Heatmap of Late − Early Braak score differences per gene and cell type
- UMAP visualization of cell-type distribution
- Prioritization of 35 AD-associated therapeutic targets

---

## Inputs

| Parameter | Description |
|-----------|-------------|
| `--obj` | Processed Seurat object (`.rds`) with `cell_braak` metadata |
| `--target` | `input_bellenguez.txt` — 35 AD-associated target genes |
| `--column` | `cell_braak` — metadata column for group comparison |

---

## Outputs

- `rank_scores/perbscore_all_targets.txt` — perturbation scores per cell type and target
- Heatmap (Late − Early Braak)
- UMAP projections

---

## Credits

Developed by Beatriz Morais, Julia Apolonio, and Diego Coelho at 
[MCBLab](https://github.com/MCBLab).
Analysis performed on the NPAD/UFRN supercomputer.

## Citations
If you use this pipeline in your research, please cite:
https://www.sciencedirect.com/science/article/pii/S266637912400260X
