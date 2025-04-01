# scToppR

An API wrapper for [ToppGene](https://toppgene.cchmc.org/)

scToppR is a package that allows seamless, workflow-based interaction with ToppGene, a portal for gene enrichment analysis. Researchers can use scToppR to directly query ToppGene's databases and conduct analysis with a few lines of code. 

Please note: The use of any data from ToppGene is governed by their [Terms of Use](https://toppgene.cchmc.org/navigation/termsofuse.jsp).

## Installation

This package can be installed from Bioonductor:
```
if (!requireNamespace('BiocManager', quietly = TRUE))
install.packages('BiocManager')

BiocManager::install('scToppR')
library(scToppR)
```

Additionally, the latest development version can be installed via GitHub:
```         
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github('BioinformaticsMUSC/scToppR')
library(scToppR)
```

## Usage

To query ToppGene and create a dataframe of results, use the function `toppFun`. The function takes as inputs a list or table of top markers (with clusters/cell types as columns) and a list of any ToppGene categories (e.g. "GeneOntologyMolecularFunction" and/or "ToppGene")

Input data: scToppR can handle 3 different types of inputs, please use the `type` parameter to select the correct format.
- `type=degs`:a dataframe containing differential expression results similar to FindAllMarkers (Seurat) or DESeq2 outputs. The dataframe needs columns to determine genes, groups of cells (e.g., clusters or celltypes, if applicable), average log fold changes, and p-values. If use bulk RNAseq data, please create a dummy cluster column (e.g. df$cluster = "bulk").
- `type=marker_df`: a dataframe with celltypes/cluster names as column names and marker genes for each as values. See `data("ifnb.markers.df")` as an example.
- `type=marker_list`: a vector of genes

The package includes example data in all three formats, using the IFNB dataset (Kang 2018) from the SeuratData package.

```         
data("ifnb.de")
toppData <- toppFun(ifnb.de,
                    gene_col = "gene",
                    cluster_col = "celltype",
                    p_val_col = "p_val_adj",
                    logFC_col = "avg_log2FC")
```

This results in a dataframe like the following (only showing the top 5 rows): 
| idx | Category | ID | Name | PValue | QValueFDRBH | QValueFDRBY | QValueBonferroni | TotalGenes | GenesInTerm | GenesInQuery | GenesInTermInQuery | Source | URL | Cluster | nlog10_fdr |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | GeneOntologyMolecularFunction | GO:0033885 | 10-hydroxy-9-(phosphonooxy)octadecanoate phosphatase activity | 0.000804141754370397 | 0.0128662680699264 | 0.0699043269565933 | 0.102930144559411 | 19912 | 1 | 16 | 1 |   |   | 0 | 1.8905474042743 |
| 2 | GeneOntologyMolecularFunction | GO:0038121 | C-C motif chemokine 21 receptor activity | 0.000804141754370397 | 0.0128662680699264 | 0.0699043269565933 | 0.102930144559411 | 19912 | 1 | 16 | 1 |   |   | 0 | 1.8905474042743 |
| 3 | GeneOntologyMolecularFunction | GO:0038117 | C-C motif chemokine 19 receptor activity | 0.000804141754370397 | 0.0128662680699264 | 0.0699043269565933 | 0.102930144559411 | 19912 | 1 | 16 | 1 |   |   | 0 | 1.8905474042743 |
| 4 | GeneOntologyMolecularFunction | GO:0047627 | adenylylsulfatase activity | 0.000804141754370397 | 0.0128662680699264 | 0.0699043269565933 | 0.102930144559411 | 19912 | 1 | 16 | 1 |   |   | 0 | 1.8905474042743 |
| 5 | GeneOntologyMolecularFunction | GO:0047352 | adenylylsulfate-ammonia adenylyltransferase activity | 0.000804141754370397 | 0.0128662680699264 | 0.0699043269565933 | 0.102930144559411 | 19912 | 1 | 16 | 1 |   |   | 0 | 1.8905474042743 |
| 6 | GeneOntologyMolecularFunction | GO:0035757 | chemokine (C-C motif) ligand 19 binding | 0.000804141754370397 | 0.0128662680699264 | 0.0699043269565933 | 0.102930144559411 | 19912 | 1 | 16 | 1 |   |   | 0 | 1.8905474042743 |

## Visualization

Once a toppData dataframe is created, it can be used to create a dotplot or balloon plot. Plots can be automatically saved, and you can enter any number of clusters and categories to produce several and save several plots at once.

Example code for dotplot

```         
toppPlot(toppData, 
         category = "GeneOntologyMolecularFunction", 
         clusters = "CD8 T")
```

![DotPlot of toppData results](/examples/toppplot_example.png)

If multiple clusters are including in the query, the function will return a list of ggplots, which can be shown using [Patchwork](https://patchwork.data-imaginist.com/) or another similar package.

Balloon plots can help visualize any overlap (or lack thereof) between the top terms for each group of cells.

Example code for balloon plot

```
toppBalloon(toppData,
            categories = "Pathway")
```

![Balloon plot of toppData results](/examples/balloon_example.png)
## Topp Categories

The available Topp Categories are:

```         
GeneOntologyMolecularFunction
GeneOntologyBiologicalProcess
GeneOntologyCellularComponent
HumanPheno
MousePheno
Domain
Pathway
Pubmed
Interaction
Cytoband
TFBS
GeneFamily
Coexpression
CoexpressionAtlas
ToppCell
Computational
MicroRNA
Drug
Disease
```

To capture these in R, run the command `get_ToppCats()`.

For more information, please see the vignettes.

## License
While scToppR is shared with a MIT License, please understand the Terms of Use set by ToppGene, which can be found [here.](https://toppgene.cchmc.org/navigation/termsofuse.jsp) Additionally, if data from ToppGene are used, please see the ToppGene citation notes [here.](https://toppgene.cchmc.org/help/publications.jsp)
