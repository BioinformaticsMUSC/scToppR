# scToppR

An API wrapper for [ToppGene](https://toppgene.cchmc.org/)

Currently, this package utilizes the ToppFun portion of the site. Lists of marker genes can be submitted to create a dataframe of results.

## Installation

This package can be installed from the Github repository:

```         
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github('BioinformaticsMUSC/scToppR')
library(scToppR)
```

## Usage

To query ToppGene and create a dataframe of results, use the function `toppFun`. The function takes as inputs a list or table of top markers (with clusters/cell types as columns) and a list of any ToppGene categories (e.g. "GeneOntologyMolecularFunction" and/or "ToppGene")

Input data: the input for `toppData` is a dataframe similar to FindAllMarkers (Seurat) or Wilcoxauc (Presto) outputs. The dataframe needs columns to determine genes, groups of cells (e.g., clusters or celltypes), average log fold changes, and p-values. 

The package includes example data in the FindAllMarkers format, using the PBMC 3K dataset from the Seurat guided clustering tutorial.

```         
data(pbmc.markers)
toppData <- toppFun(pbmc.markers, 
                    topp_categories = NULL,
                    cluster_col = "cluster",
                    gene_col = "gene"

            )
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
toppPlot(toppData, category = "ToppCell", clusters = "1", save = T, num_terms = 10)
```

![DotPlot of toppData results](/examples/toppplot_example.png)

If multiple clusters are including in the query, the function will return a list of ggplots, which can be shown using [Patchwork](https://patchwork.data-imaginist.com/) or another similar package.

Balloon plots can help visualize any overlap (or lack thereof) between the top terms for each group of cells.

Example code for balloon plot

```
toppBalloon(toppData,
            categories = c("GeneOntologyMolecularFunction"),
            balloons = 2,
            height = 6, 
            width=15)
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
