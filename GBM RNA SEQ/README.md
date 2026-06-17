Glioblastoma Transcriptomic Analysis (GSE4290)

Overview
Differential gene expression analysis of Glioblastoma (GBM) versus normal brain tissue using GEO dataset GSE4290.

Methods
- GEOquery
- limma
- EnhancedVolcano
- pheatmap
- clusterProfiler
- KEGG enrichment

Workflow
1. Download GSE4290 from GEO
2. Extract GBM and normal samples
3. Differential expression analysis
4. Volcano plot visualization
5. Heatmap of top DEGs
6. GO enrichment analysis
7. KEGG pathway enrichment

Key Results
- Significant DEGs: 7718
- Upregulated genes: 4153
- Downregulated genes: 3565

Files
- scripts/GBM_analysis_clean.R
- results/
- figures/

Author
Sanjana Kori