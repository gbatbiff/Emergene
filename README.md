# EMERGENe - estimating emergence rate of bacterial traits of epidemiological interest

## Introduction
This package estimates the emergence rate of traits of epidemiological interest pairing ancestral state reconstruction and phyletic patterns of phylogenetic tree shape

## Pipeline description

- An ancestral state reconstruction is performed through stochastic character mapping using the ”make.simmap” function from the R package ”phytools” (Revell L, 2024).        Internal phylogenetic tree nodes that experienced a shift between the states i.e. when a trait is gained are spotted.
- Once that the history of a trait is reconstructed, the number of introductions of the trait across the population are detected, intended as Polyphyletic event
- Once a parent node is detected, its children node with opposite state become the MRCA and a sub-tree is generated. Then, the Entry rate is calculated based on
  the branch distance between the parent node and its child with opposite state i.e. branch duration associated with the trait introduction, where a transition between       demes occurs
- Finally, the Emergence Rate is computed, indicating the propagation time of the trait in the population that descend from the initial node that experienced a shift         between the states. (See methods section for details) 

## Installation
### Requirements
- Linux-based OS
- conda 

### Installation command

```bash
git clone https://github.com/gbatbiff/Emergene.git
cd EMERGENe/
conda env create -f environment.yml
conda activate EMERGENe
```
## Quick guide
The standard inputs are a time scaled phylogenetic tree (Bacdating, BEAST...), and the output from AMRFinder

The command with default settings is:
```bash
Rscript Emergene.R -t [treefile] -amr [AMRFinderPlus_output_table]
```

### Test run
```bash
Rscript emergene_test_run.R -t test/phylo/tree.nex -amr test/metadata/amr.csv
```

## Notes
The tree labels have to match the strain name of AMRFinder output table 

## EMERGENe output summary
The EMERGENe script generates single outputs for each gene screened from AMRFinder results with the main following attributes:

- `amr` - resistance/virulence genes identified by AMRFinder
- `node_lineages` - number of internal nodes descendants that gained the trait
- `tip_lineages` - number of terminal tips descendants that gained the trait
- `poly_parent` - internal node that experienced a polyphyletic event i.e. a shift between states where a trait is gained
- `polyphyly` - number of polyphyletic events
- `coalescent_interval` - branching distance calculated from polyphyletic parent node and its descendants
- `poly_parent_nodeheight` - node height (age) of the polyphyletic node
- `entry_rate` - trait introduction (see method section for details)
- `R_vs_S_prop` - percentage of descendants of the polyphyletic node that carry the trait
- `poly_parent_state_prob` - ancestral state probability of character assignment (Default >0.8)
- `emergence_rate` - post-introduction trait propagation (see method section for details)

## Output processing
Rscript to merge and filter multiple outputs for easy visualization through ggplot2
```bash
Rscript processing_EMERGENe.R [output_folder]
```


### Citation

