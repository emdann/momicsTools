A collection of resources and knowledge for analysis of scATAC+scRNA
---
### Table of contents

- [scATAC analysis](https://github.com/emdann/momicsTools#scatac-analysis)
- [Joint dimensionality reduction](https://github.com/emdann/momicsTools#joint-dimensionality-reduction)
- [Publicly available multi-omic datasets](https://github.com/emdann/momicsTools#publicly-available-multi-omic-datasets)

---
### scATAC analysis
<!-- 
#### Barcode Multiplets
This was raised as a problem in sn 10X protocols 
-->

#### Picking features
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


#### Quality control

* Calculating QC metrics
    - 

* Filtering/refining peaks

* Filtering cells: 
    - [Doublet detection in ArchR](https://www.archrproject.com/bookdown/how-does-doublet-identification-work-in-archr.html) 

#### Dimensionality reduction
[blogpost on dimensionality reduction in scATAC-seq](http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/)

#### Defining gene-level features
Summarising accessibility over gene promoters/bodies is often useful for data exploration and for qualitative or quantitative comparison with RNA. Different pipelines/papers use different strategies to reduce accessibility signal to a gene x cell matrix. 

- Counting fragments over gene bodies and promoters (implemented in [Signac](https://satijalab.org/signac/reference/GeneActivity.html))
- [ArchR gene scoring](https://www.archrproject.com/bookdown/calculating-gene-scores-in-archr.html): uses counts over gene bodies and promoter + weighting for distance of regulatory elements around the gene
- [Cicero gene activity score](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/##cicero-gene-activity-scores): based on measuring co-accessibility between peaks
- [Gene scores from cisTopic](https://www.embopress.org/doi/full/10.15252/msb.20209438): taking the average de-noised accessibility signal from peaks around each gene ([example implementation](https://github.com/emdann/scATAC_prep/blob/master/N2_add_cistopic.ipynb))

#### Motif analysis 
- [chromVAR](https://github.com/GreenleafLab/chromVAR): determine variations in chromatin accessibility across peaks containing a set of TF binding motifs

---
### Joint dimensionality reduction

- Multi-omics factors analysis: ([paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1)) ([code](https://github.com/bioFAM/MOFA2)) ([vignette](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/10x_scRNA_scATAC.html))
- Seurat V4 Weighted Nearest neighbor analysis ([paper](https://www.biorxiv.org/content/10.1101/2020.10.12.335331v1))([code](https://github.com/satijalab/seurat))([vignette](https://satijalab.org/seurat/v4.0/weighted_nearest_neighbor_analysis.html))

---
### Publicly available multi-omic datasets

_(supposedly we will have more papers soon, add newest on top)_

- [Example data from 10X](https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets)
- [SHARE-seq data](https://www.cell.com/cell/fulltext/S0092-8674(20)31253-8?rss=yes) 
- [SNARE-seq data](https://www.nature.com/articles/s41587-019-0290-0)
- [sciCAR-data](https://science.sciencemag.org/content/361/6409/1380)












