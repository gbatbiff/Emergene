# Emergene - estimating emergence rate of bacterial traits of epidemiological interest

## Introduction
This package estimates the emergence rate of traits of epidemiological interest pairing ancestral state reconstruction and phyletic patterns of phylogenetic tree shape

## Pipeline description

- An ancestral state reconstruction is performed through stochastic character mapping using the ”make.simmap” function from the R package ”phytools” (Revell L, 2024).        Internal phylogenetic tree nodes that experienced a shift between the states i.e. when a trait is gained are spotted.
- Once that the history of a trait is reconstructed, the number of introductions of the trait across the population are detected, intended as Polyphyletic event
- Once a parent node is detected, its children node with opposite state become the MRCA and a sub-tree is generated. Then, the Entry rate is calculated based on
  the branch distance between the parent node and its child with opposite state i.e. branch duration associated with the trait introduction, where a transition between       demes occurs
- Finally, the Emergence Rate is computed, indicating the propagation time of the trait in the population that descend from the initial node that experienced a shift         between the states. (See Methods for details) 

## Installation
### Requirements
- Linux-based OS
- conda 

### Installation command

```bash
git clone https://github.com/gbatbiff/Emergene.git
cd Emergene/
conda env create -f environment.yml
conda activate Emergene
```
## Quick guide
The standard inputs are a time scaled phylogenetic tree (Bacdating, BEAST...), and the output from AMRFinder

The command with default settings is:
```bash
Rscript Emergene.R -t [treefile] -amr [AMRFinderPlus_output_table]
```

### Test run
```bash
Rscript Emergene.R -t data/treefile.nwk -amr data/AMRFinder.txt
```


## Notes
The tree labels have to match the strain name of AMRFinder output table 

## Filtering EMERGENe output
The EMERGENe script generates single outputs for each gene screened from AMRFinder results with the main following attributes:

' - amr
' - node lineages
' - tip lineages



### Citation

