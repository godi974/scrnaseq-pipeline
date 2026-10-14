# scRNA-seq pipeline — PBMC 3k
**Data:** ~2,700 peripheral blood mononuclear cells (PBMCs) from one healthy donor, 10x Genomics Chromium 3′ (v1 chemistry), processed with Cell Ranger 1.1.0 (hg19). Loaded with `sc.datasets.pbmc3k()`.

**Goal:** run a complete single-cell RNA-seq analysis, from raw counts to annotated cell types, justifying each choice, before turning it into a reproducible pipeline.

**Steps:**
1. Load and explore the data
2. Quality control (MAD-based filtering of cells, filtering of genes)
3. Doublet detection (Scrublet)
4. Normalization and log transformation
5. Feature selection (highly variable genes)
6. Scaling and PCA
7. Neighbor graph, UMAP and Leiden clustering
8. Marker genes per cluster
9. Cell type annotation
