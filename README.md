# Emergene - an ancestral state coalescent based framework to estimate emergence rate of bacterial traits of epidemiological interest 

## Introduction
This package estimates the emergence rate of traits of epidemiological interest pairing ancestral state reconstruction and coalescence theory.

## Pipeline description

- An ancestral state reconstruction is performed to spot the internal nodes that experienced a shift in the states i.e. when a trait is gained
- Once a parent node is detected, its children node with opposite state become the MRCA and a sub-tree is generated
- The probability density of coalescence is calculated by dividing the sub-tree into fixed time frames calculated on the phylogenetic distance from the MRCA to the youngest terminal tip of 
  the sub-clade
