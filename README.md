# CombSAFE

## About
This repository contains the implementation of the *Combinatorial and Semantic Analysis of Functional Elements* (CombSAFE) presented in: ```Leone M., Galeota E., Ceri S., Masseroli M., and Pellizzola M. "Identification, semantic annotation and comparison of regulatory functional element combinations in multiple biological conditions", 2021```. It is a flexible computational method to identify combinations of static and dynamic functional elements genome-wide, and how they change across semantically annotated biological conditions. 

Given as input a set of ChIP-seq dataset samples and the list of functional elements to be considered, the ```CombSAFE``` pipeline allows:
- considering histone modifications, transcription factors and any other type of dynamic and static genomic features (e.g., CpG islands, partially methylated domains, transposable elements, etc.)
- relying on public repositories to retrieve the considered static functional elements and identify from the input input dataset metadata the biological conditions in which the dynamic functional elements of interest were charted
- identifying combinations of static and dynamic functional elements  throughout the genome in the corresponding omics data using hidden Markov models 
- focusing on specific genomic regions, applying clustering to explore how significant combinations of the functional elements compare among cell types and biological conditions 
- performing functional enrichment analyses based on the genes found in genomic regions with similar combinations of functional elements.

## Cookbook

### Load input files
```combsafe.import_path(filepath)```<br/>
Parameters: 
* filepath, path object or file-like object

### Generate semantic annotations
```combsafe.generate_semantic_df(separator, encode_convert)```<br/>
Parameters: <br/>
  separator: str, default '\t' <br/>
Delimiter to use
* encode_convert: bool, default False
  * If true, id encode IDs are searched to be converted to GSM <br/>
Returns: 
  - DataFrame or TextParser
    - A comma-separated values (csv) file is returned as two-dimensional data structure with labeled axes.

### Data Analysis
```combsafe.plot_factor_freq(merged_df, n)```<br/>

![alt text](https://drive.google.com/uc?export=download&id=1WyFjK1eYM9nSbMKLht0dXp6ouscZ381P)

```combsafe.generate_fixed_factor_pool(semantic_df, ["H3K27ac", "H3K4me3", "H3K27me3"], 5)```

![alt text](https://drive.google.com/uc?export=download&id=1TD-wc-4rJ0DDagZebLu0BgFeLVJ8SwW0)

```combsafe.get_semantic_annotation_list(semantic_df, ["H3K27ac", "H3K4me3", "H3K27me3", "H3K4me1", "H3K36me3"])```

![alt text](https://drive.google.com/uc?export=download&id=1llQnJyeJku6evCgDaOymWuiIgCE5dYXO)


### Select features and combine samples

```combsafe.run_gmql(["H3K27ac", "H3K4me3", "H3K27me3", "H3K4me1", "H3K36me3"])```

### [Optional] Add custom tracks

```combsafe.add_custom_tracks(track_lable, path_to_custom_track, index)```

### Identify Chromatin States <> NON USIAMO PIù QUESTO TERMINE MA FUNCTIONAL STATES -> MODIFICARE 

```combsafe.identify_chromatin_states(number_of_states, n_core)```

### Whole Genome Analysis

...

### Single Gene Analysis

...

### [Optional] Semantic Analysis

...
