
A collection of resources on analysis of scATAC+scRNA multi-omic data

## Table of contents

- [scATAC analysis](https://github.com/emdann/momicsTools#scatac-analysis)
- [Data structures for multi-omics](https://github.com/emdann/momicsTools#data-structures-for-multi-omics)
- [Joint dimensionality reduction](https://github.com/emdann/momicsTools#joint-dimensionality-reduction)
- [Peak-gene matching](https://github.com/emdann/momicsTools#peak-gene-matching)
- [Publicly available multi-omic datasets](https://github.com/emdann/momicsTools#publicly-available-multi-omic-datasets)

## scATAC analysis
<!-- 
#### Barcode Multiplets
This was raised as a problem in sn 10X protocols 
-->

### Picking features
Because Tn5 insertions can take place anywhere in the genome, after alignment for each cell barcode we will have a set of fragments covering different genomic positions. To allow comparison between cells we need to define a common set of genomic features. Once a feature set is defined, we count the number of fragments overlapping each feature for each single-cell, obtaining a (quasi-)binary feature x cell count matrix.

Different pipelines/papers use different strategies to define a feature set. 

* Binning the genome into equally sized windows (5-10kb) (e.g. [SnapATAC](https://github.com/r3fang/SnapATAC/blob/master/examples/10X_PBMC_15K/README.md#add_bmat)). Using a set of genomic bins is a good way to start exploring your data, however when working with a large genome (e.g. human and mouse) this generates a huge matrix with over 200k features, making many analysis steps much slower. In addition, epigenomics studies highlight that the tipical size of regulatory regions is closer to hundreds rather than thousands of bps. 

* Peak calling: take a pseudo-bulk coverage track of all the cells or a group of cells and identify peaks in this signal. The most used algorithm for peak calling is [MACS](https://github.com/macs3-project/MACS). 
    - [Peak calling from CellRanger](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/algorithms/overview#peaks) (big caveat: does peak calling for each sample independently)
    - [Cusanovich2018 approach](https://www.cell.com/cell/fulltext/S0092-8674(18)30855-9?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418308559%3Fshowall%3Dtrue): starts by binning the genome into fixed-size windows and building a bin x cell binary matrix. Bins that overlap ENCODE-defined blacklist regions are filtered out, and the top 20,000 most commonly used bins are retained. Then, the bins-by-cells binary matrix is normalized and rescaled using the term frequency-inverse document frequency (TF-IDF) transformation. Next, singular value decomposition (SVD) is performed to generate a PCs-by-cells LSI score matrix, which is used to cluster cells. Within each cluster, peak calling is performed on the aggregated scATAC-seq profiles

* Using a annotated set of enhancers/regulatory regions in your genome of interest (frequently seen for studies in Drosophila)

* Other scATAC-specific methods for feature extraction:

    - [Scregseg](https://github.com/BIMSBbioinfo/scregseg): uses an HMM with Dirichlet-Multinomial emission probabilities to segment the genome according to distinct relative cross-cell accessibility profiles ([paper](https://www.biorxiv.org/content/10.1101/2020.06.26.173377v1))
    - [BROCKMAN]: k-mer based accessibility ([paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2255-6))

### Quality control

* Calculating QC metrics ([Signac](https://satijalab.org/signac/articles/pbmc_vignette.html#computing-qc-metrics-1))([ArchR](https://www.archrproject.com/bookdown/plotting-sample-fragment-size-distribution-and-tss-enrichment-profiles-.html))
    - TSS enrichment: we expect an enrichment of fragments around transcription start sites
    - Fragment size distribution: we expect to see regular peaks, these correspond to fragments that span a different number of nucleosomes

* Filtering/refining peaks: this is dependent on the peak calling algorithms, but you might want to filter some peaks out (also to make the size of the dataset more manageable):
    - Remove peaks overlapping with [ENCODE blacklist](https://www.nature.com/articles/s41598-019-45839-z) (an R-ready annotation for model organism genomes can be found in [Signac](https://satijalab.org/signac/reference/index.html#section-data))
    - Remove peaks that are accessible in less than n% cells in any given cluster (strategy used in [Cusanovish2018](https://www.cell.com/cell/fulltext/S0092-8674(18)30855-9?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418308559%3Fshowall%3Dtrue))
    - Filter by peak width: peaks that are exactly as wide as the peak-calling width threshold are usually just a few reads from a few cells. On the other side, keeping very wide peaks might introduce artifacts (bigger peak, more reads)

* Filtering cells: 
    - Basic filtering of cells with very low coverage (these often mess up dimensionality reduction if unfiltered)
    - [ArchR doublet detection](https://www.archrproject.com/bookdown/how-does-doublet-identification-work-in-archr.html) 

### Dimensionality reduction

* TF-IDF + SVD (Latent Semantic Indexing) ([Signac](https://satijalab.org/signac/articles/pbmc_vignette.html#normalization-and-linear-dimensional-reduction-1)) ([ArchR](https://www.archrproject.com/bookdown/archrs-lsi-implementation.html))([muon](https://muon.readthedocs.io/en/latest/omics/atac.html#normalisation))
* Latent Dirichlet Allocation ([cisTopic](http://htmlpreview.github.io/?https://github.com/aertslab/cisTopic/blob/master/vignettes/WarpLDA_CompleteAnalysis.html))
* PeakVI ([paper](https://www.biorxiv.org/content/10.1101/2021.04.29.442020v1))([code](https://github.com/YosefLab/scvi-tools/blob/master/scvi/model/_peakvi.py))

Additional resources:
* [Blogpost on dimensionality reduction in scATAC-seq](http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/)

### Defining gene-level features
Summarising accessibility over gene promoters/bodies is often useful for data exploration and for qualitative or quantitative comparison with RNA. Different pipelines/papers use different strategies to reduce accessibility signal to a gene x cell matrix. 

- Counting fragments over gene bodies and promoters (implemented in [Signac](https://satijalab.org/signac/reference/GeneActivity.html))
- [ArchR gene scoring](https://www.archrproject.com/bookdown/calculating-gene-scores-in-archr.html): uses counts over gene bodies and promoter + weighting for distance of regulatory elements around the gene
- [Cicero gene activity score](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/##cicero-gene-activity-scores): based on measuring co-accessibility between peaks
- [Gene scores from cisTopic](https://www.embopress.org/doi/full/10.15252/msb.20209438): taking the average de-noised accessibility signal from peaks around each gene ([example implementation](https://github.com/emdann/scATAC_prep/blob/master/N2_add_cistopic.ipynb))

### Motif analysis 
- [chromVAR](https://github.com/GreenleafLab/chromVAR): determine variations in chromatin accessibility across peaks containing a set of TF binding motifs

## Data structures for multi-omics

- MUON data ([code](https://github.com/gtca/muon))([preprint](https://www.biorxiv.org/content/10.1101/2021.06.01.445670v1.full.pdf)) - python - extension of [AnnData](https://anndata.readthedocs.io/en/latest/)
- MultiAssayExperiment object ([code](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html)) - R/Bioconductor - extension of [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html)

## Joint dimensionality reduction

- j-SNE and j-UMAP ([paper](https://www.biorxiv.org/content/10.1101/2021.01.10.426098v1)) ([code](https://github.com/canzarlab/JVis-learn))
- Multi-omics factors analysis: ([paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1)) ([code](https://github.com/bioFAM/MOFA2)) ([vignette](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/10x_scRNA_scATAC.html))
- Seurat V4 Weighted Nearest Neighbor analysis ([paper](https://www.cell.com/cell/fulltext/S0092-8674%2821%2900583-3))([preprint](https://www.biorxiv.org/content/10.1101/2020.10.12.335331v1))([code](https://github.com/satijalab/seurat))([vignette](https://satijalab.org/seurat/v4.0/weighted_nearest_neighbor_analysis.html)) - extends to > 2 modalities
- Schema ([paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02313-2)) ([code](https://schema-multimodal.readthedocs.io/en/latest/overview.html))
- [Multigrate](https://icml-compbio.github.io/2021/papers/WCBICML2021_paper_44.pdf)
- MultiVI ([paper](https://www.biorxiv.org/content/10.1101/2021.08.20.457057v1.full))([vignette](https://docs.scvi-tools.org/en/stable/user_guide/notebooks/MultiVI_tutorial.html))

## Peak-gene matching

- DORC analysis ([example code](https://github.com/buenrostrolab/stimATAC_analyses_code/blob/master/R/runDORCs_stim.R))

## Publicly available multi-omic datasets

_(add newest on top)_

- [Fleck et al. 2021](https://www.biorxiv.org/content/10.1101/2021.08.24.457460v1) - Multiome data of human cerebral organoids
- [Mimitou et al. 2021](https://www.nature.com/articles/s41587-021-00927-2) - 10X genonics multiome + surface proteins, demonstrated on PBMCs
- [Dou et al. 2020](https://www.biorxiv.org/content/10.1101/2020.12.11.422014v1.full.pdf) 9383 mouse retina cells, used as gold standard for diagonal integration algorithm 
- [Trevino et al. 2020](https://www.biorxiv.org/content/10.1101/2020.12.29.424636v1.full.pdf) Chromatin and gene-regulatory dynamics of the developing human cerebral cortex at single-cell resolution (8981 cells)
- [Example data from 10X](https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets)
- [SHARE-seq data](https://www.cell.com/cell/fulltext/S0092-8674(20)31253-8?rss=yes) 
- [SNARE-seq data](https://www.nature.com/articles/s41587-019-0290-0)
- [sciCAR-data](https://science.sciencemag.org/content/361/6409/1380)












