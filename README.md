# CombSAFE

## About
This repository contains the implementation of the *Combinatorial and Semantic Analysis of Functional Elements* (CombSAFE) presented in: ```Leone M, Galeota E, Masseroli M, Pellizzola M. "Identification, semantic annotation and comparison of combinations of functional elements in multiple biological conditions", 2021```. It is a flexible computational method to identify combinations of static and dynamic functional elements genome-wide, or in a specific genomic region, and how they change across semantically annotated biological conditions. 

Given as input a set of ChIP-seq dataset samples and the list of functional elements to be considered, the ```CombSAFE``` pipeline allows:
- considering histone modifications, transcription factors and any other type of dynamic or static genomic features (e.g., CpG islands, partially methylated domains, transposable elements, etc.)
- relying on public repositories to retrieve the considered static functional elements and identify from the input dataset metadata the biological conditions in which the dynamic functional elements of interest were charted
- leveraging natural language processing techniques and biomedical ontologies to complement the identified conditions with semantic annotations about tissue and disease types
- identifying combinations of static and dynamic functional elements throughout the genome in the corresponding omics data, using hidden Markov models 
- focusing on specific genomic regions, applying clustering to explore how significant combinations of the functional elements compare among semantically annoted cell types and desease/healthy conditions 
- performing functional enrichment analyses based on the genes found in genomic regions with similar combinations of functional elements.

## Structure
```
CombSAFE/
|-- README.md
|-- LICENSE
|-- .gitignore
|-- notebook/
|   |-- Functional_states_analysis.ipynb
|-- gene_list/
|   |-- MYC_associated.txt
|   |-- test_list.txt
|   |-- tumor_suppressor.txt
|-- CombSAFE/
|   |-- CombSAFE.py
|-- CombSAFE.yml
```

- `README.md` this file
- `LICENSE` MIT license file
- `.gitignore` standard .gitignore file for Python projects
- `notebook/Functional_states_analysis.ipynb` Python notebook to run a functional state analysis
- `gene_list/` folder where gene name lists are stored for the CombSAFE single gene analysis
- `gene_list/test_list.txt` list of random genes
- `gene_list/tumor_suppressor.txt` list of tumor suppressor genes
- `gene_list/MYC_associated.txt` list of MYC interacting genes
- `CombSAFE/CombSAFE.py` core Python routines called from within the notebook to perform the functional state analysis
- `CombSAFE.yml` dependency yaml file to load in order to perform the CombSAFE analysis


## How to install
In order to run the ```CombSAFE``` pipeline, please load the Conda environment with the command: ```conda env create -f CombSAFE.yml``` <br/>

NB: The `pyGMQL` package additionally requires Java. Please follow the installation procedure [here](https://pygmql.readthedocs.io/en/latest/installation.html). <br/>

NB2: Before using `PyEnsembl`, download and install Ensembl data. Please follow the installation procedure [here](https://pypi.org/project/pyensembl/).

## Cookbook
In the following, we show how to call the functions implemented to easily perform the different steps of our ```CombSAFE``` computational method, providing example resuls for some of them. 

### Load input files
```combsafe.import_path(filepath)```<br/>
Load the input path for the analysis.<br/>

Parameters: 
- ***filepath***: path object or file-like object 

Example:
```python
>> import_path("./path_to_files/")
```

---

### Generate semantic annotations
```combsafe.generate_semantic_df(separator, encode_convert)```<br/>
Generate semantic annotations about tissue and disease types from the input dataset.<br/>

Parameters: 
- ***separator***: str, default '\t'
  - Delimiter to use for the input .txt file
- ***encode_convert***: bool, default False
  - If true, encode IDs are searched to be converted to GSM

Returns: 
  - ***DataFrame*** or ***TextParser***
    - A comma-separated values (csv) file is returned as a two-dimensional data structure with labeled axes.

Example:
```python
>> semantic_df = generate_semantic_df(sep="\t", encode_convert=True)
```

---

### Data analysis
```combsafe.plot_factor_freq(dataframe, n)```<br/>
Vertical barplot of the factor frequency in the input dataset.<br/>

Parameters: 
- ***dataframe***: str, default '\t'
  - Dataframe of semantic annotations
- ***n***: int
  - Number of factors to diplay in the barplot

Example:
```python
>> plot_factor_freq(semantic_df, 30)
```

![alt text](https://drive.google.com/uc?export=download&id=1WyFjK1eYM9nSbMKLht0dXp6ouscZ381P) <br/>

---

```combsafe.generate_fixed_factor_pool(dataframe, factor_list, number_of_semantic_annotation)``` <br/>
Table containg lists of factors according to the selected parameters.<br/>

Parameters: 
- ***dataframe***: dataframe
  - Dataframe of semantic annotations
- ***factor_list***: list
  - List of factors to include in the analysis
- ***number_of_semantic_annotation***: int
  - Number of semantic annotations to include in the analysis

Example:
```python
>> generate_fixed_factor_pool(semantic_df, ["CTCF", "MYC"], 5)
```

![alt text](https://drive.google.com/uc?export=download&id=1Qc4W9vm2ekev_P13-56akRNpK_oY92BQ)

---

```combsafe.get_semantic_annotation_list(dataframe, factor_list)``` <br/>
List of semantic annotations according to the selected factors.<br/>

Parameters: 
- ***dataframe***: dataframe
  - Dataframe of semantic annotations
- ***factor_list***: list
  - List of factors to include in the analysis

Example:
```python
>> get_semantic_annotation_list(semantic_df, ["CTCF", "MYC", "POLR2A", "H3K4me3", "H3K27me3"])  <br/>
```

![alt text](https://drive.google.com/uc?export=download&id=1llQnJyeJku6evCgDaOymWuiIgCE5dYXO)

---

### Data extraction and replica combination
```combsafe.run_gmql(factor_list)```<br/>
Combine sample replicas of the listed factors and extract their semantic annotations regarding the conditions in which they were mapped.<br/>

Parameters: 
- ***factor_list***: list
  - List of factors to include in the analysis

Example:
```python
>> run_gmql(["CTCF", "MYC", "POLR2A", "H3K4me3", "H3K27me3"])
```
---

### [Optional] Add custom tracks

```combsafe.download_file("custom_tracks_link")```<br/>
Download custom tracks of static genomic elements from URL (e.g., UCSC) in the ./input_files/ folder. <br/>

Parameters: 
- ***custom_tracks_link***: url
  - url of custom track to be downloaded

Example:
```python
>> download_file("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz")
```

---

```combsafe.add_custom_tracks(track_label_name, path_to_custom_tracks, index)```<br/>
Add custom tracks of static genomic elements to the analysis (e.g., CpG islands). <br/>

Parameters: 
- ***track_label_name***: string
  - name of custom tracks
- ***path_to_custom_tracks***: path
  - UCSC path for downloading custom tracks
- ***index***: int
  - column to use for row labels of the DataFrame

Example:
```python
>> add_custom_tracks("CpG_Islands", "./input_files/cpgIslandExt.txt", index=1)
```

---

### Identification of combinations of genomic functional elements 

```combsafe.identify_functional_states(number_of_states, n_core)```<br/>
Add custom tracks of static genomic elements to the analysis (e.g., CpG islands). <br/>

Parameters: 
- ***number_of_states***: int
  - number of combinations of genomic functional elements
- ***n_core***: int
  - number of cores to use for the analysis

Example:
```python
>> identify_functional_states(track_lable_n, path_to_custom_track, index)
```

---

```combsafe.show_emission_graph(custom_palette=colors)```<br/>
Add custom tracks of static genomic elements to the analysis (e.g., CpG islands). <br/>

Parameters: 
- ***custom_palette***: list of exadecimals
  - optionally, add a list of customized colors in hexadecimal form to be assigned to the functional states

Example:
```python
>> colors = ['#c9f9ff', '#e6beff', '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231','#911eb4', '#bcf60c', '#f032e6', '#fffac8', '#fabebe', '#9a6324', '#46f0f0', '#008080']
>> show_emission_graph(custom_palette=colors)
```
![alt text](https://drive.google.com/uc?export=download&id=1Kk_vOm5wz_ski-fLvTxB48dhhu9TXcNY)

---




### Genome-wide analysis
```combsafe.data_driven_heatmap(functional_states_dataframe)```<br/>
Show a genome-wide heatmap with the most significant clusters of genomic regions based on their patterns of functional states. <br/>

Parameters: 
- ***functional_states_dataframe***: dataframe
  - dataframe of functional states

Return: 
- ***cluster_indices***: array
  - cluster integer labels for each data sample

Example:
```python
>> genome_wide_heatmap = data_driven_heatmap(functional_states_df)
```
![alt text](https://drive.google.com/uc?export=download&id=1jbyS_WY54SfJtCWQhw9tpiYW8vC2QJ_Q)

---


```combsafe.gene_ontology_enrichment_analysis(cluster_indices, functional_state_dataframe, significance_cut_off)```<br/>
Show a genome-wide heatmap with the most significant clusters of genomic regions based on their patterns of functional states. <br/>

Parameters:
- ***cluster_indices***: array
  - cluster integer labels for each data sample
- ***functional_state_dataframe***: dataframe
  - dataframe of functional states
- ***significance_cut_off***: int
  - threshold for p-value enrichment analysis

Example:
```python
>> gene_ontology_enrichment_analysis(genome_wide_heatmap, functional_states_df, 0.05)
```
---

### Single-gene analysis

```combsafe.single_gene_analysis(functional_state_dataframe, path_to_gene_list)```<br/>
Given a list of gene symbols in a textual file, the heatmap of the functional states of the related genomic regions is shown. <br/>

Parameters:
- ***functional_state_dataframe***: dataframe
  - dataframe of functional states
- ***path_to_gene_list***: path
  - path to the gene list 

Example:
```python
>> single_gene_analysis(functional_states_df, "path_to_gene_list/gene_list.txt")
```
![alt text](https://drive.google.com/uc?export=download&id=1zkj4DhgfR36UiIAM99ohF1byaXojKRNU)

---
