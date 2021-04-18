# CostOfFuel-Pipeline

## About
This repository contains the implementation of the COmbinatorial and SemanTic analysis of FUnctional ELements (CostofFuel) pipeline presented in the paper: ```Leone M., Galeota E., Ceri S., Masseroli M., and Pellizzola M. "Identification, semantic annotation and comparison of chromatin regulatory state combinations in multiple biological condiions", 2021 ```. This pipeline is a flexible computational method to identify combinations of static and dynamic functional elements, and how they change among biological conditions. 
Given as input a set of ChIP-seq samples and the functional elements to be considered, the ```COSTofFUEL``` pipeline allow to:

- consider histone modifications, transcription factors and any other type of dynamic and static genomic features (e.g., CpG islands, transposable elements, ecc.). 
- rely on public repositories to identify from their metadata the biological conditions in which the dynamic functional elements of interest were charted
- identify combinations throughout the genome of the corresponding homogeneously analysed omic data using hidden Markov models. 
- focus on specific genomic regions, applying clustering to explore how significant combinations of the functional elements compare among cells and conditions. 
- perform functional enrichment analyses based on the genes found in regions with similar combinations of functional elements.
