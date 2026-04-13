# scRNA-seq-Assignment4

# Introduction
Influenza A virus is a major and recurring cause of respiratory disease worldwide, contributing substantially to seasonal illness, hospitalization, and mortality each year. In addition to annual outbreaks, Influenza A viruses circulate in animal reservoirs such as birds and pigs, allowing genetic reassortment and cross-species transmission that can lead to pandemic strains with serious public health consequences [1]. Because of its continued global impact and ability to rapidly evolve, Influenza A remains an important model for studying host-pathogen interactions and immune defense mechanisms.

The nasal mucosa represents one of the first anatomical barriers encountered by inhaled respiratory pathogens and plays a central role in limiting infection before spread to the lower airways. Beyond its functions in olfaction, filtration, and air conditioning, the nasal cavity contains diverse specialized cell populations that contribute to immune surveillance and tissue protection [2]. These include epithelial cells that form the physical barrier, fibroblasts that support tissue structure and repair, endothelial cells involved in vascular regulation, neurons associated with olfactory sensing, and multiple immune cell populations such as macrophages, monocytes, neutrophils, natural killer (NK) cells, B cells, and T cells [2]. During viral infection, these cell types interact through cytokine signaling, antigen presentation, antibody production, and antiviral effector responses to coordinate pathogen clearance. Previous studies have shown that Influenza A infection can induce dynamic shifts in nasal cell composition, including early interferon-responsive neutrophils, recruitment of monocyte-derived macrophages, activation of lymphocytes, and establishment of tissue-resident immune memory following viral resolution [2].

Traditional bulk RNA sequencing methods measure average gene expression across mixed cell populations, masking important differences between individual cell types and transient cellular states. In contrast, single-cell RNA sequencing (scRNA-seq) profiles gene expression at single-cell resolution, allowing researchers to characterize cellular heterogeneity, identify rare populations, and detect condition-specific transcriptional responses that would otherwise be obscured [3,4]. This technology has become widely used in biomedical research because it enables detailed investigation of tissue microenvironments, developmental trajectories, and disease-associated immune responses [3,4].

Selecting an appropriate computational framework for scRNA-seq analysis requires balancing accessibility, reproducibility, analytical performance, and compatibility with modern datasets. Commonly used frameworks include Seurat and Scanpy, which together account for the majority of scRNA-seq analyses and implement highly similar standard workflows [5]. Standard scRNA-seq workflows include quality-control filtering, normalization, highly variable gene selection, dimensionality reduction, clustering, visualization, and marker gene detection [5]. Because scRNA-seq enables identification of both common and rare cell populations at single-cell resolution, it has become an important tool for cellular profiling and biomarker discovery [6].

Although Scanpy offers strong memory efficiency for Python-based pipelines, Seurat was selected for this study because it provides a widely adopted and integrated workflow for preprocessing, clustering, visualization, and differential expression within a single R environment. Seurat has been successfully applied to identify biologically meaningful cell populations using functions such as `RunPCA()`, `RunUMAP()`, `FindNeighbors()`, `FindClusters()`, and `FindAllMarkers()` [6]. Recent versions have also improved support for large datasets, multimodal integration, and pseudobulk analyses [5]. Its extensive documentation and compatibility with downstream tools such as SingleR and clusterProfiler further support reproducibility and ease of interpretation [6].

Quality-control thresholds were selected to remove low-quality or damaged cells while retaining biologically informative observations. Cells with fewer than 200 detected genes were excluded to reduce the influence of empty droplets or poor RNA capture, as low gene counts are indicative of poor-quality cells with insufficient transcript recovery [7]. Cells with mitochondrial transcript percentages greater than 15% were also removed, as elevated mitochondrial RNA content is widely associated with cellular stress, apoptosis, and membrane damage [8]. Filtering based on mitochondrial content is a standard step in scRNA-seq preprocessing pipelines, where cells with high mitochondrial proportions are routinely excluded to improve data quality [8]. These quality-control metrics, including gene counts and mitochondrial percentage, are commonly used to identify low-quality cells, although specific thresholds may vary depending on dataset characteristics [7]. 

Principal component analysis (PCA) is commonly used in scRNA-seq analysis as an initial dimensionality reduction step because gene expression datasets contain thousands of features per cell and substantial technical noise [9]. PCA summarizes the dominant sources of variation into a smaller number of components, reducing data complexity while retaining the most informative biological signal for downstream clustering and visualization [9].

Uniform Manifold Approximation and Projection (UMAP) was selected for visualization because nonlinear dimensionality reduction methods are better suited than linear approaches for representing complex single-cell datasets with thousands of measured features [10]. Although t-distributed stochastic neighbor embedding (t-SNE) has historically been widely used in single-cell analysis to identify distinct cell populations, it has several limitations, including loss of large-scale intercluster relationships, slower computation time, and reduced interpretability for very large datasets [10]. In comparison, UMAP has been shown to preserve local neighborhood structure while retaining more global relationships between clusters, resulting in embeddings that better reflect continuity between related cell states [1]. UMAP is also faster to compute, more reproducible, and scales efficiently to large scRNA-seq datasets [1]. These advantages make UMAP particularly suitable for visualizing cellular heterogeneity and cluster relationships in modern single-cell transcriptomic studies [10].

Cell type annotation is a critical step in scRNA-seq analysis because transcriptionally defined clusters must be assigned biological identities before meaningful interpretation can be made. Although manual annotation using canonical marker genes is widely used, it can be subjective, time-consuming, and less reliable when closely related populations share overlapping markers or when rare cell types are present [11]. Several automated supervised methods have been developed for this purpose, including Seurat label transfer, CellAssign, CHETAH, scmap, and SingleR [11,12]. SingleR was selected for this study because benchmarking analyses have shown it to be among the top-performing annotation methods overall, with particularly strong performance in distinguishing highly similar cell types, identifying rare populations, and maintaining accuracy as the number of cell classes increases [11]. Unlike anchor-based approaches such as Seurat, SingleR compares query expression profiles directly to reference transcriptomes using pseudo-bulk profiles, which can reduce noise and improve robustness for complex datasets [11]. Because the murine nasal mucosa contains diverse immune, stromal, epithelial, and neuronal populations with potential transcriptional overlap, SingleR provided a reproducible and objective framework for annotation, while marker gene visualization was retained as a complementary validation step. 

For downstream differential expression analysis, macrophage populations were selected because they are central mediators of innate antiviral immunity during influenza infection. Macrophages contribute to pathogen sensing, cytokine production, phagocytosis, antigen presentation, and orchestration of subsequent adaptive immune responses. Comparison between D02 and D14 samples was chosen to capture transcriptional differences between an earlier stage of infection and a later stage associated with recovery or immune memory. This temporal contrast allows investigation of how macrophage function changes over the course of infection.

To interpret differentially expressed genes in a broader biological context, over-representation analysis (ORA) was performed using clusterProfiler. Enrichment analysis identifies pathways and biological processes that are statistically overrepresented among significant genes, allowing clearer interpretation than examining gene lists alone. This approach is particularly useful for identifying coordinated antiviral, inflammatory, and immune regulatory programs associated with different stages of infection.

In this study, single-cell RNA sequencing data from murine nasal mucosa following Influenza A infection were analyzed to characterize cellular composition and transcriptional responses over time. By combining clustering, dimensionality reduction, differential expression analysis, and functional enrichment, this work aims to identify major cell populations within the nasal mucosa and determine how key immune cell types, particularly macrophages, respond during early versus later stages of infection.

# Methods
## Computational Environment
All analyses were performed in R (v4.5.1) using RStudio. Single-cell RNA sequencing analyses were conducted using Seurat (v5.4.0) [13]. Additional downstream analyses were performed using clusterProfiler(v4.16.0), SingleR(v2.10.0), celldex(v1.18.0), ggplot2(v4.0.2), dplyr(v1.2.0), tidyverse(v2.0.0), org.Mm.eg.db(v3.21.0), EnhancedVolcano(v1.26.0), SingleCellExperiment(v1.30.1), and enrichplot(v1.28.4). All analyses were performed in an R Markdown workflow to ensure reproducibility.

## Data Acquisition and Input

A preprocessed Seurat object containing murine nasal mucosa single-cell RNA sequencing data following Influenza A virus exposure was provided by Dr. Taylor. The dataset contained gene expression counts and associated metadata including sample identity (`orig.ident`), tissue region (`organ_custom`), and disease labels. Data were imported directly into R using the `readRDS()` function [14].

## Quality Control and Filtering

Initial cell quality assessment was performed using three standard metrics: number of detected genes per cell (`nFeature_RNA`), total transcript counts per cell (`nCount_RNA`), and percentage of mitochondrial transcripts (`percent.mt`) [15]. Mitochondrial percentages were calculated using genes beginning with the mouse mitochondrial prefix `mt-` [15]. Quality control metrics were visualized using violin plots and feature scatter plots. Cells with fewer than 200 detected genes or mitochondrial transcript percentages greater than 15% were removed to exclude low-quality, stressed, or dying cells.

## Downsampling, Normalization and Feature Selection

To improve computational efficiency while retaining balanced representation across samples, the filtered dataset was downsampled to 8,000 cells per sample prior to downstream analysis. Cells were log-normalized using the `LogNormalize` method with a scale factor of 10,000 [16]. Highly variable genes were identified using the variance stabilizing transformation (`vst`) method, retaining the top 2,000 most variable genes for downstream analysis. Expression values were then scaled prior to dimensionality reduction.

## Dimensionality Reduction

Principal component analysis (PCA) was performed using the variable genes identified during feature selection [7]. The number of informative principal components was assessed using an elbow plot. Based on the inflection point of the variance explained curve, the first 20 principal components were retained for downstream clustering and visualization.

## Graph-Based Clustering and UMAP Visualization

Cells were embedded into a shared nearest-neighbor graph using the first 20 principal components. Clusters were identified using Seurat’s graph-based clustering algorithm with a resolution parameter of 0.5. Uniform Manifold Approximation and Projection (UMAP) was then applied using the same 20 principal components to generate a two-dimensional representation of transcriptional similarity among cells [17]. UMAP projections were visualized by cluster identity and sample/timepoint (orig.ident) to assess clustering structure and potential batch effects.

## Automated Cell Type Annotation

The Seurat object was converted to a SingleCellExperiment object for compatibility with SingleR [18]. Automated cell type annotation was performed using SingleR with the MouseRNAseq reference dataset from celldex [19,20]. Predicted labels were added to the Seurat metadata and used as identities for downstream visualization. Annotation confidence was assessed using score heatmaps and delta distribution plots. Predicted cell populations included neurons, epithelial cells, fibroblasts, macrophages, monocytes, B cells, NK cells, granulocytes, and endothelial cells.

## Marker Gene Identification and Validation

To support annotation quality, marker genes were identified using Seurat’s `FindAllMarkers()` function on a downsampled subset of 500 cells [21]. Only positively enriched genes expressed in at least 25% of cells with log fold-change greater than 0.25 were retained. The top five markers per annotated cell type were extracted. Canonical marker genes were further visualized using `FeaturePlot()` and `VlnPlot()`, including Cd3d, Cd8a, Ms4a1, Lyz2, Csf1r, Epcam, Krt8, Ncr1, and Itgax [22,23].

## Cell Composition Analysis

Changes in cellular composition across samples were assessed by calculating both proportional abundance and absolute counts of annotated cell types per sample using metadata summaries generated with dplyr and visualized using stacked bar plots in ggplot2 [24].

## Differential Expression Analysis

Macrophage populations were selected for focused downstream analysis due to their established role in innate antiviral immunity. Cells annotated as macrophages were subsetted, and sample identities were reassigned according to collection timepoint using the `orig.ident` metadata field [15]. Differential expression analysis was performed between early infection (D02) and later infection (D14) using Seurat’s `FindMarkers()` function with `min.pct = 0.1` and `logfc.threshold = 0` [25]. Genes with adjusted p-values less than 0.05 were considered statistically significant. Differential expression results were visualized using volcano plots generated with EnhancedVolcano [26].

## Functional Enrichment Analysis

Significantly differentially expressed genes were separated into genes upregulated in D02 and D14 macrophages and analyzed for Gene Ontology (GO) Biological Process over-representation using clusterProfiler with the mouse annotation database org.Mm.eg.db [27,28]. Multiple-testing correction was performed using the Benjamini–Hochberg method, and terms with adjusted p-values below 0.05 were considered significantly enriched. Dot plots were used to visualize enriched pathways.

## Gene Set Enrichment Analysis

Gene Set Enrichment Analysis (GSEA) was additionally performed using ranked log2 fold-change values from all macrophage differential expression results with `gseGO()` [29]. Significant pathways were visualized using dot plots, ridge plots, and enrichment score plots generated with enrichplot.

## Heatmaps and Output Files

The top 20 differentially expressed macrophage genes were visualized using Seurat `DoHeatmap()`, grouped by sample identity [30]. Intermediate results, processed Seurat objects, and output tables were exported in  `.csv` and `.rds` formats for reproducibility and downstream reporting.

# Results

## Quality control and preprocessing

<img width="1388" height="851" alt="violin_plot_before_QC" src="https://github.com/user-attachments/assets/737df4ba-264c-417b-b5fd-501186370189" />

**Figure 1. Quality-control metrics prior to filtering of single-cell RNA-seq data.** Distributions of the number of detected genes per cell (`nFeature_RNA`), total transcript counts per cell (`nCount_RNA`), and percentage of mitochondrial transcripts (`percent.mt`) are shown across samples (D02, D05, D08, D14, and Naive) using violin plots. Each violin represents the density of cells across observed values for the indicated metric.

Quality-control metrics were examined prior to filtering to assess dataset consistency and identify low-quality cells (Figure 1). The distributions of detected genes and transcript counts were broadly similar across all samples, indicating comparable sequencing depth and capture efficiency between conditions. Most cells contained approximately 900–3500 detected genes and 3,500–6,000 transcript counts, although some higher-count outliers were present. Mitochondrial transcript percentages were generally low across samples, with the majority of cells below 5%. However, a subset of cells in each sample showed elevated mitochondrial content approaching 15%, consistent with stressed or damaged cells. These results supported the use of filtering thresholds to remove poor-quality observations prior to downstream analysis.

<img width="1388" height="854" alt="image" src="https://github.com/user-attachments/assets/afd0c693-5d8e-44a7-9e7a-2af0ca5707c0" />

**Figure 2. Relationship between transcript counts and detected genes prior to filtering.** Scatter plot showing the relationship between total transcript counts per cell (`nCount_RNA`) and the number of detected genes per cell (`nFeature_RNA`) across all samples. Each point represents a single cell and is colored according to sample identity (D02, D05, D08, D14, and Naive).

The relationship between sequencing depth and gene complexity was assessed prior to filtering (Figure 2). A clear positive association was observed, with cells containing higher transcript counts generally exhibiting larger numbers of detected genes, consistent with expected high-quality single-cell libraries. Cells from all five samples overlapped extensively within the main density region, indicating comparable data quality across experimental groups and no obvious sample-specific technical separation. A smaller population of outlier cells with high transcript counts but relatively low gene complexity was also observed, which may represent damaged cells, ambient RNA contamination, or multiplets. These results supported the application of downstream quality-control filtering.

<img width="1384" height="856" alt="violin_plot_after_QC" src="https://github.com/user-attachments/assets/47a8791c-6ac1-4bd9-9570-d0bddd2293ac" />

**Figure 3. Quality-control metrics after filtering of single-cell RNA-seq data.** Violin plots showing the distributions of detected genes per cell (`nFeature_RNA`), total transcript counts per cell (`nCount_RNA`), and percentage of mitochondrial transcripts (`percent.mt`) across samples (D02, D05, D08, D14, and Naive) following quality-control filtering. Each violin represents the density of cells for the indicated metric within each sample.

Quality-control metrics remained broadly consistent across all samples following filtering (Figure 3). Similar distributions of detected genes and transcript counts were observed between timepoints, indicating comparable sequencing depth and library complexity across conditions. Mitochondrial transcript percentages were low in all groups, with most cells below 5%, consistent with retention of high-quality cells. Compared with the pre-filtering dataset, only minor changes were observed in the overall distributions, suggesting that relatively few cells failed the applied thresholds and that the provided dataset had already undergone substantial preprocessing or consisted primarily of high-quality cells suitable for downstream analysis.

<img width="1382" height="854" alt="image" src="https://github.com/user-attachments/assets/2454ffb1-23b1-4b9d-bd0a-3fedde96fab4" />

**Figure 4. Elbow plot for principal component selection.** Elbow plot showing the standard deviation explained by each principal component (PC) following principal component analysis of the 2,000 highly variable genes. The first 30 PCs are displayed to assess the contribution of successive components to total transcriptional variance.

Principal component analysis was used to summarize the major sources of variation in the dataset prior to clustering and visualization (Figure X). The elbow plot showed a steep decline in variance explained across the first several PCs, followed by a more gradual decrease after approximately PC10 and a plateauing trend by around PC15–20. This pattern indicates that the earliest components captured the strongest biological signal, whereas later components contributed progressively less additional information. Based on this distribution, the first 20 principal components were retained for downstream nearest-neighbor graph construction, clustering, and UMAP visualization.

**Table 1. Summary of graph-based clustering parameters and outcomes.** Key parameters and outputs from unsupervised clustering analysis of the filtered single-cell RNA-seq dataset. Cells were clustered using a shared nearest-neighbor graph constructed from the first 20 principal components and partitioned using the Louvain community detection algorithm at a resolution of 0.5.

Unsupervised graph-based clustering was performed on 40,000 cells using the first 20 principal components (Table 1). Application of the Louvain algorithm at a resolution of 0.5 identified 30 transcriptionally distinct clusters. The high modularity score (0.9561) indicated strong separation between communities, supporting the presence of substantial cellular heterogeneity within the nasal mucosa dataset prior to downstream annotation and differential expression analyses.

<div align="center">
  
| Parameter | Value |
|-----------|-------|
| Cells analyzed | 40,000 |
| Principal components used | 20 |
| Clustering algorithm | Louvain |
| Resolution | 0.5 |
| Number of clusters identified | 30 |
| Modularity score | 0.9561 |
</div>

<img width="1390" height="856" alt="UMAP_cluster" src="https://github.com/user-attachments/assets/373a5988-ab35-4909-8c89-86151f12ac14" />

**Figure 5. UMAP visualization of unsupervised Seurat clusters.** Uniform Manifold Approximation and Projection (UMAP) of 40,000 filtered cells based on the first 20 principal components. Cells are colored according to graph-based Seurat cluster assignments generated using the Louvain algorithm at a resolution of 0.5. Cluster labels (0–29) indicate the 30 transcriptionally distinct communities identified.

UMAP visualization revealed clear transcriptional structure within the dataset, with cells separating into 30 discrete and well-defined clusters (Figure 5). Several large clusters contained substantial numbers of cells, while smaller peripheral clusters likely represented rarer cell populations. The clear separation between clusters is consistent with the high modularity score obtained during graph-based clustering and indicates substantial cellular heterogeneity within the nasal mucosa. These unsupervised clusters formed the basis for subsequent cell-type annotation using SingleR and downstream differential expression analyses.

<img width="1384" height="856" alt="UMAP_sample_identity" src="https://github.com/user-attachments/assets/8325ea06-112b-46a8-9bf4-0e7c8eb0beb2" />

**Figure 6. UMAP visualization colored by sample identity.** UMAP of the filtered single-cell dataset colored according to sample/timepoint identity (Naive, D02, D05, D08, and D14). Cells from all samples are projected onto the same low-dimensional embedding to assess mixing across biological conditions and potential batch effects.

Cells from all five samples were broadly distributed across the same UMAP structure, with substantial overlap observed within most clusters (Figure 6). This indicates that clustering was driven primarily by shared transcriptional cell states rather than by sample-specific technical effects. Although modest enrichment of certain samples was observed within selected regions, no major sample-exclusive clusters were detected. These results suggest minimal batch effects and support the suitability of the integrated dataset for comparative downstream analyses across infection timepoints.

**Table 2. Predicted cell-type abundance based on SingleR annotation**. Cell identities were assigned using the SingleR reference-based annotation framework with the MouseRNAseqData reference dataset. Values indicate the number of cells assigned to each predicted cell type across the filtered dataset.

Automated cell-type annotation identified a diverse mixture of epithelial, immune, stromal, and neuronal populations within the nasal mucosa dataset (Table 2). Neurons represented the largest predicted population (12,592 cells), followed by epithelial cells (5,521) and fibroblasts (5,290), indicating substantial contributions from structural and sensory cell types. Among immune populations, monocytes (3,300), B cells (2,806), macrophages (2,709), granulocytes (2,012), T cells (1,050), and NK cells (737) were also detected, consistent with active immune surveillance and inflammatory responses in nasal tissue. Smaller populations including endothelial cells, dendritic cells, erythrocytes, and rare glial-associated labels were also present. Overall, these results demonstrate marked cellular heterogeneity and support the presence of both resident tissue cells and infiltrating immune populations across infection timepoints

<div align="center">
  
| Cell Type | Number of Cells |
|-----------|----------------:|
| Neurons | 12,592 |
| Epithelial cells | 5,521 |
| Fibroblasts | 5,290 |
| Monocytes | 3,300 |
| B cells | 2,806 |
| Macrophages | 2,709 |
| Endothelial cells | 2,354 |
| Granulocytes | 2,012 |
| T cells | 1,050 |
| Hepatocytes | 963 |
| NK cells | 737 |
| Adipocytes | 236 |
| Microglia | 176 |
| Oligodendrocytes | 88 |
| Erythrocytes | 85 |
| Dendritic cells | 44 |
| Cardiomyocytes | 19 |
| Astrocytes | 18 |
</div>

<img width="1394" height="864" alt="UMAP_SingleR" src="https://github.com/user-attachments/assets/042c6533-240f-4404-9210-b7f58a3732e1" />

**Figure 7. UMAP visualization of SingleR-annotated cell populations.** UMAP projection of 40,000 filtered cells colored according to predicted cell identities assigned using the SingleR reference-based annotation framework with the MouseRNAseqData reference. Labels indicate the dominant predicted population occupying each transcriptionally distinct region.

Reference-based annotation assigned biologically meaningful identities to the unsupervised clusters identified by Seurat (Figure X). Major populations included neurons, epithelial cells, fibroblasts, monocytes, macrophages, B cells, endothelial cells, granulocytes, T cells, and NK cells. Large neuronal and epithelial compartments were evident, alongside multiple immune-associated populations distributed across distinct transcriptional regions. Closely related monocyte and macrophage populations occupied neighboring areas of the embedding, consistent with related myeloid lineages. Overall, the annotated UMAP confirmed substantial cellular heterogeneity within the nasal mucosa and provided a biologically interpretable framework for downstream comparative analyses.

**Table 3. Representative marker genes identified for annotated cell populations.** Top differentially expressed marker genes identified using Seurat `FindAllMarkers()` on a downsampled dataset. Marker genes are shown for major predicted cell populations and represent transcripts enriched relative to other cells in the dataset.

Marker gene analysis further supported the accuracy of automated cell-type annotation by identifying canonical lineage-associated transcripts within each predicted population (Table 3). B-cell populations were characterized by strong enrichment of **Cd79a**, **Cd79b**, **Igkc**, and **Ighm**, consistent with B-cell receptor signaling and immunoglobulin expression. Myeloid populations including monocytes and macrophages showed elevated expression of genes such as **Lyz2**, **Fcer1g**, and **C1qa**, while fibroblast populations were enriched for extracellular matrix-associated genes including **Col1a1** and **Col1a2**. Endothelial cells expressed vascular markers such as **Pecam1** and **Kdr**, whereas epithelial populations were marked by keratin genes including **Krt8** and **Krt18**. Together, these lineage-specific transcriptional signatures were consistent with known biology and validated the SingleR-derived cell identity assignments.

<div align="center">
  
| Cell Type | Representative Marker Genes |
|-----------|-----------------------------|
| B cells | Cd79a, Cd79b, Igkc, Ighm |
| T cells | Cd3d, Cd3e, Trbc2 |
| Macrophages | C1qa, C1qb, Lyz2, Fcer1g |
| Monocytes | S100a8, S100a9, Lyz2 |
| Fibroblasts | Col1a1, Col1a2, Sparc |
| Endothelial cells | Pecam1, Kdr, Emcn |
| Epithelial cells | Krt8, Krt18, Krt19 |
| Neurons | Snap25, Tubb3, Elavl4 |
</div>

<img width="1388" height="860" alt="image" src="https://github.com/user-attachments/assets/3bb3947b-5150-473c-af48-3e414fd65c0f" />

**Figure 8. Canonical marker gene expression across annotated cell populations.** Violin plots showing normalized expression distributions of representative lineage markers across SingleR-annotated cell types. **Lyz2** marks myeloid populations, **Cd3d** marks T cells, **Ms4a1** marks B cells, and **Epcam** marks epithelial cells.

Expression of canonical lineage markers further supported the accuracy of automated cell-type annotation (Figure 8). **Lyz2** expression was highest in macrophage, monocyte, and granulocyte populations. **Cd3d** was selectively enriched in T cells, while **Ms4a1** showed strong and specific expression in B-cell populations. **Epcam** expression was highest among epithelial populations and low in immune-associated clusters. These cell-type-restricted expression patterns provided independent validation of the SingleR-derived annotations.

<img width="1392" height="856" alt="cell_counts" src="https://github.com/user-attachments/assets/3ea45d8b-28a5-4886-a523-d36fb1c24f55" />

**Figure 9. Cell type composition across influenza infection timepoints.** Stacked bar plots show the proportional abundance of the ten most abundant SingleR-annotated cell types across samples (D02, D05, D08, D14, and Naive), with remaining low-abundance populations grouped as “Other.” Each bar represents one sample, and segment height corresponds to the proportion of cells assigned to each annotated cell type.

Cell type composition was broadly similar across all samples, with neuronal, fibroblast, epithelial, and B-cell populations representing the largest fractions of the dataset (Figure 9). Neurons comprised the dominant population in most samples, although their relative abundance was reduced in D14 compared with other timepoints. In contrast, D14 showed increased proportions of immune-associated populations, including macrophages, monocytes, granulocytes, and B cells, suggesting a stronger inflammatory or recovery-associated immune response at this later stage of infection. Earlier post-infection samples (D02–D08) displayed intermediate compositions, with relatively stable proportions of structural cell populations such as fibroblasts and epithelial cells. Overall, these results indicate dynamic shifts in cellular composition over time following influenza infection, with the most pronounced immune enrichment observed at D14

<p align="center">
  <img width="696" height="430" alt="image" src="https://github.com/user-attachments/assets/6e60bddc-32e2-4866-8b80-cfcf55bb33d3" />
</p>

**Figure 10. Differential gene expression in macrophages between D02 and D14 samples.** Volcano plot showing genes differentially expressed in macrophages comparing D02 and D14 timepoints. Each point represents one gene, positioned by log₂ fold change (x-axis) and adjusted significance (−log₁₀ adjusted P value; y-axis). Genes meeting both adjusted P < 0.05 and |log₂ fold change| > 0.5 are highlighted in red, while genes meeting only one threshold are shown in blue or orange. Selected significant genes are labelled.

Differential expression analysis revealed marked transcriptional differences in macrophages between early (D02) and later (D14) infection stages (Figure X). Multiple genes were significantly upregulated in D02 relative to D14, including several ribosomal protein genes (_Rps27_, _Rps28_,_ Rps21_, _Rps29_, _Rpl37_, and _Rpl39_), suggesting increased translational activity or heightened cellular activation during the early immune response. Additional genes such as _Hspa1a_ and _Gpr65_ were also elevated in D02, consistent with stress-response and inflammatory signaling pathways. In contrast, fewer genes showed strong upregulation in D14, with _Obp1a_ representing one of the most prominent later-stage markers. Overall, these results indicate that macrophages undergo substantial temporal transcriptional remodeling during influenza infection, shifting from an activated early-response state at D02 toward a different functional profile by D14.

<p>
  <img width="49%" src="https://github.com/user-attachments/assets/17abe23a-eee0-4522-a597-d2bb2fdaab13" alt="ORA_upregulated_D02" />
  <img width="49%" src="https://github.com/user-attachments/assets/d447eb14-fe83-469a-b617-d06a0588b17c" alt="ORA_upregulated_D14" />
</p>

**Figure 11. Gene Ontology biological process enrichment of differentially expressed macrophage genes.** Dot plots show over-representation analysis (ORA) of significantly upregulated genes in macrophages comparing D02 and D14 samples. The left panel shows biological processes enriched among genes upregulated in D02, while the right panel shows processes enriched among genes upregulated in D14. Dot size represents the number of genes associated with each term, x-axis indicates gene ratio, and color reflects adjusted P value.

Functional enrichment analysis revealed distinct biological programs associated with early and later infection stages in macrophages (Figure 11). Genes upregulated at D02 were enriched for processes related to cytoplasmic translation, gliogenesis, iron ion transport, and regulation of ferroptosis, suggesting an early response characterized by active protein synthesis, metabolic adaptation, and stress-related signaling. In contrast, genes upregulated at D14 were strongly enriched for ribonucleoprotein complex biogenesis, ribosome biogenesis, protein folding, heat response, and antigen processing and presentation via MHC class II pathways. These later-stage enrichments indicate a transition toward cellular recovery, proteostasis, and enhanced antigen presentation functions. Overall, the results support dynamic temporal remodeling of macrophage activity during influenza infection, with early inflammatory/metabolic responses followed by later immune regulatory and presentation-associated programs.


# Discussion

# References
[1] Fauci, A. S., & Collins, F. S. (2012). Benefits and Risks of Influenza Research: Lessons Learned. Science (American Association for the Advancement of Science), 336(6088), 1522–1523. https://doi.org/10.1126/science.1224305

[2] Kazer, S. W., Match, C. M., Langan, E. M., Messou, M.-A., LaSalle, T. J., O’Leary, E., Marbourg, J., Naughton, K., von Andrian, U. H., & Ordovas-Montanes, J. (2024). Primary nasal influenza infection rewires tissue-scale memory response dynamics. Immunity (Cambridge, Mass.), 57(8), 1955-1974.e8. https://doi.org/10.1016/j.immuni.2024.06.005

[3] Wang, Y., Li, Z., & Lu, J. (2024). Single-cell RNA sequencing reveals the epithelial cell, fibroblast, and key gene alterations in chronic rhinosinusitis with nasal polyps. Scientific Reports, 14(1), Article 2270. https://doi.org/10.1038/s41598-024-52341-8

[4] Pan, M., Jia, Z., Zhou, M., Chen, J., Qiu, H., Luo, X., Zhang, Y., Shi, Z., Wu, S., Wang, D., & Yang, Q. (2026). The Application of Single-cell RNA Sequencing Technology in the Research of Sino-Nasal Diseases. Clinical Reviews in Allergy & Immunology, 69(1), Article 18. https://doi.org/10.1007/s12016-026-09145-7

[5] Rich, J. M., Moses, L., Einarsson, P. H., Jackson, K., Luebbert, L., Booeshaghi, A. S., Antonsson, S., Sullivan, D. K., Bray, N., Melsted, P., & Pachter, L. (2026). The impact of package selection and versioning on single-cell RNA-seq analysis. Cell Systems, Article 101560. https://doi.org/10.1016/j.cels.2026.101560

[6] Wang, G., Zhang, E., Chen, A., & Meng, D. (2024). Single-cell RNA-seq analysis revealed the stemness of a specific cluster of B cells in acute lymphoblastic leukemia progression. PeerJ (San Francisco, CA), 12, Article e18296. https://doi.org/10.7717/peerj.18296

[7] Amezquita, R. A., Lun, A. T. L., Becht, E., Carey, V. J., Carpp, L. N., Geistlinger, L., Marini, F., Rue-Albrecht, K., Risso, D., Soneson, C., Waldron, L., Pagès, H., Smith, M. L., Huber, W., Morgan, M., Gottardo, R., & Hicks, S. C. (2020). Orchestrating single-cell analysis with Bioconductor. Nature Methods, 17(2), 137–145. https://doi.org/10.1038/s41592-019-0654-x

[8] Osorio, D., & Cai, J. J. (2021). Systematic determination of the mitochondrial proportion in human and mice tissues for single-cell RNA-sequencing data quality control. Bioinformatics, 37(7), 963–967. https://doi.org/10.1093/bioinformatics/btaa751

[9] Tsuyuzaki, K., Sato, H., Sato, K., & Nikaido, I. (2020). Benchmarking principal component analysis for large-scale single-cell RNA-sequencing. Genome Biology, 21(1), Article 9. https://doi.org/10.1186/s13059-019-1900-3

[10] Becht, E., McInnes, L., Healy, J., Dutertre, C.-A., Kwok, I. W. H., Ng, L. G., Ginhoux, F., & Newell, E. W. (2019). Dimensionality reduction for visualizing single-cell data using UMAP. Nature Biotechnology, 37(1), 38–44. https://doi.org/10.1038/nbt.4314

[11] Huang, Q., Liu, Y., Du, Y., & Garmire, L. X. (2021). Evaluation of Cell Type Annotation R Packages on Single-cell RNA-seq Data. Genomics, Proteomics & Bioinformatics, 19(2), 267–281. https://doi.org/10.1016/j.gpb.2020.07.004

[12] Cheng, C., Chen, W., Jin, H., & Chen, X. (2023). A Review of Single-Cell RNA-Seq Annotation, Integration, and Cell–Cell Communication. Cells (Basel, Switzerland), 12(15), 1970. https://doi.org/10.3390/cells12151970

[13] Satijalab. (n.d.). Satijalab/Seurat: R toolkit for single cell genomics. GitHub. https://github.com/satijalab/seurat 

[14] ReadRDS: Serialization interface for single objects. RDocumentation. (n.d.). https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/readRDS 

[15] Mary Piper, L. P. (2020, February 24). Single-cell RNA-seq: Quality control analysis. Introduction to Single-cell RNA-seq - ARCHIVED. https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html 

[16] Normalize raw data - lognormalize. - LogNormalize • Seurat. (n.d.). https://satijalab.org/seurat/reference/lognormalize 

[17] Run umap - runumap. - RunUMAP • Seurat. (n.d.). https://satijalab.org/seurat/reference/runumap 

[18] Drisso. (n.d.). Drisso/Singlecellexperiment: Clone of the Bioconductor Repository for the singlecellexperiment package, see https://bioconductor.org/packages/devel/bioc/html/singlecellexperiment.html for the official development version. GitHub. https://github.com/drisso/SingleCellExperiment 

[19] SingleR-inc. (n.d.). Singler-Inc/celldex: Collection of cell type reference datasets. GitHub. https://github.com/SingleR-inc/celldex 

[20] Dviraran. (n.d.). Dviraran/Singler: Singler: Single-cell RNA-seq cell types recognition (legacy version). GitHub. https://github.com/dviraran/singler 

[21] Gene expression markers for all identity classes - findallmarkers. - FindAllMarkers • Seurat. (n.d.). https://satijalab.org/seurat/reference/findallmarkers 

[22] Visualize “features” on a dimensional reduction plot - featureplot. - FeaturePlot • Seurat. (n.d.). https://satijalab.org/seurat/reference/featureplot 

[23] Single cell violin plot - vlnplot. - VlnPlot • Seurat. (n.d.). https://satijalab.org/seurat/reference/vlnplot 

[24] Holtz, Y. (n.d.). Data visualization with R and GGPLOT2: The R Graph Gallery. Data visualization with R and ggplot2 | the R Graph Gallery. https://r-graph-gallery.com/ggplot2-package.html 

[25] Gene expression markers of identity classes - findmarkers. - FindMarkers • Seurat. (n.d.). https://satijalab.org/seurat/reference/findmarkers 

[26] Kevinblighe. (n.d.). Kevinblighe/Enhancedvolcano: Publication-ready volcano plots with enhanced colouring and labeling. GitHub. https://github.com/kevinblighe/enhancedvolcano 

[27] YuLab-SMU. (n.d.). Yulab-SMU/Clusterprofiler: :bar_chart: A universal enrichment tool for interpreting OMICS data. GitHub. https://github.com/YuLab-SMU/clusterProfiler 

[28] Org.Mm.eg.db. Bioconductor. (n.d.). https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html 

[29] GSEGO: GSEGO. RDocumentation. (n.d.-a). https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/gseGO 

[30] Feature expression heatmap - doheatmap. - DoHeatmap • Seurat. (n.d.). https://satijalab.org/seurat/reference/doheatmap 
