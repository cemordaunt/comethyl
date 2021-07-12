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
You can install comethyl from this repository with the code below.

```
install.packages(c("BiocManager", "remotes"))
BiocManager::install("cemordaunt/comethyl")
```

## Workflow
![Comethyl Workflow](man/figures/comethyl.png)
