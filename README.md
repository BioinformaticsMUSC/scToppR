# scToppR

An API wrapper for [ToppGene](https://toppgene.cchmc.org/)

Currently, this package utilizes the ToppFun portion of the site. Lists of marker genes can be submitted to create a dataframe of results.

## Installation

This package can be installed from the Github repository:

```         
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github('BioinformaticsMUSC/scToppR')
```

## Usage

To query ToppGene and create a dataframe of results, use the function `toppFun`. The function takes as inputs a list or table of top markers (with clusters/cell types as columns) and a list of any ToppGene categories (e.g. "GeneOntologyMolecularFunction" and/or "ToppGene")

Input types: you can submit a single list or character vector of gene symbols, e.g. `c("IFNG", "FOXP3"...)` or a data.frame with column names as cluster/celltype labels. A table can also be provided. scToppR can also handle the "top_markers" data.frames from the package [Presto](https://github.com/immunogenomics/presto) as inputs.

If a table is provided, a data.frame of results for all clusters will be returned.

```         
toppData <- toppFun(top_markers, topp_categories = c("GeneOntologyMolecularFunction",
                                              "GeneOntologyBiologicalProcess"),
             max_results = 5)
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

Once a toppData dataframe is created, it can be used to create a dotplot.

Example code:

```         
toppPlot(toppData, category = "GeneOntologyMolecularFunction", clusters = "X0")
```

![DotPlot of toppData results](/examples/toppplot_example.png)

If multiple clusters are including in the query, the function will return a list of ggplots, which can be shown using [Patchwork](https://patchwork.data-imaginist.com/) or another similar package.

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
