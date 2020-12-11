A collection of resources and knowledge for analysis of scATAC+scRNA
---
### Table of contents

- [scATAC analysis](https://github.com/emdann/momicsTools#scatac-analysis)
- [Joint dimensionality reduction](https://github.com/emdann/momicsTools#joint-dimensionality-reduction)
- [Publicly available multi-omic datasets](https://github.com/emdann/momicsTools#publicly-available-multi-omic-datasets)

---
### scATAC analysis

#### Barcode Multiplets
This was raised as a problem in sn 10X protocols

#### Picking features
A. Using genomic bins

B. Peak calling

C. Other methods for feature extraction
    - [Scregseg](https://github.com/BIMSBbioinfo/scregseg): uses an HMM with Dirichlet-Multinomial emission probabilities to segment the genome according to distinct relative cross-cell accessibility profiles ([paper](https://www.biorxiv.org/content/10.1101/2020.06.26.173377v1))

#### Quality control
**Calculating QC metrics**

**Filtering/refining peaks**

###### Filtering cells: 
Doublet detection etc etc

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












