# comethyl
An R package for weighted region comethylation network analysis.

Comethyl builds upon the WGCNA package to identify and interpret modules of 
comethylated regions from whole-genome bisulfite sequencing data. Regions are 
defined from clusters of CpG sites or from genomic annotations, and then percent
methylation values are used to identify comethylation modules. 

Interesting modules are identified and explored by comparing with sample traits and 
examining functional enrichments. Results are then visualized with high-quality,
editable plots from ggplot2.

## Installation
You can install Comethyl from this repository and load it into your R session with the code below.

```
install.packages(c("BiocManager", "remotes"))
BiocManager::install("cemordaunt/comethyl")
library(comethyl)
```

## Documentation
Complete documentation for comethyl is available at [https://cemordaunt.github.io/comethyl/](https://cemordaunt.github.io/comethyl/).

- [Introduction to Comethyl](https://cemordaunt.github.io/comethyl/articles/comethyl.html)

- [Function reference](https://cemordaunt.github.io/comethyl/reference/index.html)

- [CpG cluster analysis vignette](https://cemordaunt.github.io/comethyl/articles/CpG_Cluster_Analysis.html)

- [Gene body analysis vignette](https://cemordaunt.github.io/comethyl/articles/Gene_Body_Analysis.html)

## Workflow
![Comethyl Workflow](man/figures/comethyl.png)

## Acknowledgements
Many thanks to Julia Mouat for creating the vignettes, and to both Ben Laufer and
Janine LaSalle for very helpful discussions.

