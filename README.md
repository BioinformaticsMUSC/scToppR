# toppHat

An API wrapper for [ToppGene](https://toppgene.cchmc.org/)

Currently, this package utilizes the ToppFun portion of the site. Lists of marker genes can be submitted to create a dataframe of results.

## Usage

To query ToppGene and create a dataframe of results, use the function `toppFun`. The function takes as inputs a table of top markers (with clusters/cell types as columns) and a list of any ToppGene categories (e.g. "GeneOntologyMolecularFunction" and/or "ToppGene")

```
toppData <- toppFun(top_markers, topp_categories = c("GeneOntologyMolecularFunction",
                                              "GeneOntologyBiologicalProcess"),
             max_results = 5)
```

This results in a dataframe like the following (only showing the top 5 rows):
| Category | ID | Name | PValue | QValueFDRBH | QValueFDRBY | QValueBonferroni | TotalGenes | GenesInTerm | GenesInQuery | GenesInTermInQuery | Source | URL | Cluster |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GeneOntologyMolecularFunction | GO:0140662 | ATP-dependent protein folding chaperone | 4.94232275200795e-10 | 9.48925968385527e-08 | 5.53917783724434e-07 | 9.48925968385527e-08 | 19912 | 42 | 20 | 5 |   |   | X0 |
| GeneOntologyMolecularFunction | GO:0044183 | protein folding chaperone | 8.55821443516093e-09 | 8.21588585775449e-07 | 4.79586968560159e-06 | 1.6431771715509e-06 | 19912 | 73 | 20 | 5 |   |   | X0 |
| GeneOntologyMolecularFunction | GO:0051082 | unfolded protein binding | 1.34585354772594e-07 | 8.61346270544602e-06 | 5.0279477334899e-05 | 2.5840388116338e-05 | 19912 | 126 | 20 | 5 |   |   | X0 |
| GeneOntologyMolecularFunction | GO:0023026 | MHC class II protein complex binding | 2.49589173174481e-06 | 9.1957075756797e-05 | 0.000536782228519332 | 0.000479211212495003 | 19912 | 27 | 20 | 3 |   |   | X0 |
| GeneOntologyMolecularFunction | GO:0002135 | CTP binding | 2.87365861739991e-06 | 9.1957075756797e-05 | 0.000536782228519332 | 0.000551742454540782 | 19912 | 3 | 20 | 2 |   |   | X0 |


## Visualization

Once a toppData dataframe is created, it can be used to create a dotplot.

Example code:

```
toppPlot(toppData, category = "GeneOntologyMolecularFunction", clusters = "X0")
```
![DotPlot of toppData results](/examples/toppplot_example.pdf)