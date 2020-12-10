A collection of resources and knowledge for analysis of scATAC+scRNA

# scATAC preprocessing

## Barcode Multiplets
This was raised as a problem in sn 10X protocols

## Picking features

A. Using genomic bins

B. Peak calling

C. Other methods for feature extraction
    - [Scregseg](https://github.com/BIMSBbioinfo/scregseg): uses an HMM with Dirichlet-Multinomial emission probabilities to segment the genome according to distinct relative cross-cell accessibility profiles ([paper](https://www.biorxiv.org/content/10.1101/2020.06.26.173377v1))

## Quality control

### Filtering/refining peaks

# Filtering cells
Doublet detection etc etc

# scRNA preprocessing

Do we need this? There is quite a body of knowledge to this.

# Dimensionality reduction

## scRNA
Classic PCA + KNN graph

## scATAC

[blogpost on dimensionality reduction in scATAC-seq](http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/)

## Joint dimensionality reduction

- Multi-omics factors analysis [MOFA2](https://github.com/bioFAM/MOFA2)
- [Seurat V4 Weighted Nearest neighbor analysis](https://satijalab.org/seurat/v4.0/weighted_nearest_neighbor_analysis.html)

# Defining gene-level features

Summarising accessibility over gene promoters/bodies is often useful for data exploration and for qualitative or quantitative comparison with RNA. Different pipelines/papers use different strategies to reduce accessibility signal to a gene x cell matrix. 

- Counting fragments over gene bodies and promoters (implemented in [Signac](https://satijalab.org/signac/reference/GeneActivity.html))
- [ArchR gene scoring](https://www.archrproject.com/bookdown/calculating-gene-scores-in-archr.html): uses counts over gene bodies and promoter + weighting for distance of regulatory elements around the gene
- [Cicero gene activity score](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/#cicero-gene-activity-scores): based on measuring co-accessibility between peaks
- [Gene scores from cisTopic](https://www.embopress.org/doi/full/10.15252/msb.20209438): taking the average de-noised accessibility signal from peaks around each gene ([example implementation](https://github.com/emdann/scATAC_prep/blob/master/N2_add_cistopic.ipynb))


# Motif analysis 
- [chromVAR](https://github.com/GreenleafLab/chromVAR)











