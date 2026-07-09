# Emergene - estimating emergence rate of bacterial traits of epidemioligical interest

## Introduction
This package estimates the emergence rate of traits of epidemiological interest pairing ancestral state reconstruction and coalescence theory.

## Pipeline description

- An ancestral state reconstruction is performed to spot the internal nodes that experienced a shift in the states i.e. when a trait is gained
- Once that the history of a trait is reconstructed, the number of introductions of the trait across the population are detected
- Once a parent node is detected, its children node with opposite state become the MRCA and a sub-tree is generated. Then, the Entry rate is calculated based on
  the branch distance between the parent node and its child with opposite state i.e. branch duration associated with the trait introduction, where a transition between       demes occurs
- Finally, the Emergence Rate is computed, indicating the propagation time of the trait in the population that descend from the initial node that experienced a shift         between the states. (See Method on paper for details) 
