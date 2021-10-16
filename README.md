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
In order to run the ```CombSAFE``` pipeline, please follow the steps below: <br/>

1. Install the Anaconda package and environment manager from [here](https://docs.anaconda.com/anaconda/install/)
2. Load the CombSAFE environment with the command: ```conda env create -f CombSAFE.yml``` <br/>
3. Activate the CombSAFE environment with the command: ```conda activate CombSAFE```
4. Run the ```notebook/Functional_states_analysis.ipynb```

NB: The `pyGMQL` package additionally requires Java. Please follow the installation procedure [here](https://pygmql.readthedocs.io/en/latest/installation.html). <br/>
NB2: Before using `PyEnsembl`, download and install Ensembl data. Please follow the installation procedure [here](https://pypi.org/project/pyensembl/). <br/>
NB3: For windows users, Visual Studio v.14 or higher is required. Please follow the installation procedure [here](https://visualstudio.microsoft.com/visual-cpp-build-tools/)

## Cookbook
In the following, we show how to call the functions implemented to easily perform the different steps of our ```CombSAFE``` computational method, providing example resuls for some of them. 

### Load input files
```combsafe.import_path(path)```<br/>
Load the input path for the analysis. Input files must be structured as follows: <br/>

```
Input_folder/
|-- Chip_Files/
|   |-- 1.narrowPeak.bed
|   |-- 2.broadPeak.bed
|   |-- 3.narrowPeak.bed
|   |-- 4.broadPeak.bed
|   |-- 5.narrowPeak.bed
|   |-- 6.broadPeak.bed
|   |-- 7.narrowPeak.bed
|   |-- ...
|-- Textual_file.txt
```
- `Chip_Files`  a folder containing ChIP-Seq files
- `Textual_file.txt` a text file containing the following information:
  - `GSMID`, accession number related to a specific ChIP-Seq sample from GEO database
  - `Factor`, Transcription Factor or Histone Mark used for the analysis
  - `File`, Filename of the corresponding ChIP-Seq File in the Chip_Files folder

| GSMID     | Factor   | File             |
| :---------| :------- | :--------------- |
| GSM648494 | H3K4me1  | 1.narrowPeak.bed |
| GSM648495 | H3K4me3  | 2.broadPeak.bed  |
| GSM575280 | H3K27me3 | 3.narrowPeak.bed |
| ...       | ...      | ...              |	

Parameters: 
- ***path***: path object or file-like object   

Example:
```python
>> input_path = import_path("./Input_folder/")
```

---

### Generate semantic annotations
```combsafe.generate_semantic_annotations(input_path, sep, encode_convert)```<br/>
Generate semantic annotations about tissue and disease types from the input dataset.<br/>

Parameters:
- ***input_path*** str
  - path of the input dataset folder
- ***sep***: str, default '\t'
  - delimiter to use for the input .txt file
- ***encode_convert***: bool, default False
  - if true, encode IDs are searched to be converted to GSM

Returns: 
  - ***semantic_dataframe***
    - dataset of biological conditions from the input metadata

Example:
```python
>> semantic_df = generate_semantic_annotations(input_path, sep="\t", encode_convert=True)
```

---

### Data analysis
```combsafe.plot_factor_freq(semantic_dataframe, n)```<br/>
Vertical barplot of the factor frequency in the input dataset.<br/>

Parameters: 
- ***semantic_dataframe***: str, default '\t'
  - dataset of biological conditions from the input metadata
- ***n***: int
  - number of factors to display in the barplot

Example:
```python
>> plot_factor_freq(semantic_df, 30)
```

![alt text](https://drive.google.com/uc?export=download&id=1WyFjK1eYM9nSbMKLht0dXp6ouscZ381P) <br/>

---

```combsafe.generate_fixed_factor_pool(semantic_dataframe, factor_list, number_of_semantic_annotations)``` <br/>
Table containg lists of factors according to the selected parameters.<br/>

Parameters: 
- ***semantic_dataframe***: dataframe
  - dataset of biological conditions from the input metadata
- ***factor_list***: list
  - list of factors to include in the analysis
- ***number_of_semantic_annotations***: int
  - number of semantic annotations to include in the analysis

Example:
```python
>> generate_fixed_factor_pool(semantic_df, ["CTCF", "MYC"], 5)
```

![alt text](https://drive.google.com/uc?export=download&id=1Qc4W9vm2ekev_P13-56akRNpK_oY92BQ)

---

```combsafe.get_semantic_annotation_list(semantic_dataframe, factor_list)``` <br/>
List of semantic annotations according to the selected factors.<br/>

Parameters: 
- ***semantic_dataframe***: dataframe
  - dataset of biological conditions from the input metadata
- ***factor_list***: list
  - list of factors to include in the analysis

Example:
```python
>> get_semantic_annotation_list(semantic_df, ["CTCF", "MYC", "POLR2A", "H3K4me3", "H3K27me3"])
```

![alt text](https://drive.google.com/uc?export=download&id=1llQnJyeJku6evCgDaOymWuiIgCE5dYXO)

---

### Data extraction and replica combination
```combsafe.extract_data(factor_list)```<br/>
Combine sample replicas of the listed factors and extract their semantic annotations regarding the conditions in which they were mapped.<br/>

Parameters: 
- ***factor_list***: list
  - list of factors to include in the analysis

Example:
```python
>> extract_data(["CTCF", "MYC", "POLR2A", "H3K4me3", "H3K27me3"])
```
---

### [Optional] Add custom tracks

```combsafe.download_custom_tracks(custom_tracks_url)```<br/>
Download custom tracks of static genomic elements from URL (e.g., UCSC) in the ./input_files/ folder. <br/>

Parameters: 
- ***custom_tracks_url***: url
  - url of custom track to be downloaded

Example:
```python
>> download_custom_tracks("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz")
```

---

```combsafe.add_custom_tracks(tracks_label, path_to_custom_tracks, index)```<br/>
Add custom tracks of static genomic elements to the analysis (e.g., CpG islands). <br/>

Parameters: 
- ***tracks_label***: string
  - name of the custom tracks to add for the analysis
- ***path_to_custom_tracks***: path
  - UCSC path for downloading custom tracks
- ***index***: int, default=0
  - column to use for row labels of the DataFrame

Example:
```python
>> add_custom_tracks("CpG_Islands", "./input_files/cpgIslandExt.txt", index=1)
```

---

### Identification of combinations of genomic functional elements 

```combsafe.identify_functional_states(ChromHMM_path, number_of_states, n_core)```<br/>
identification of combinations of static and dynamic functional elements throughout the genome. <br/>

Parameters:
- ***ChromHMM_path***: path
  - path to the chromHMM software folder. It can be downloaded [here](http://compbio.mit.edu/ChromHMM/).
- ***number_of_states***: int
  - number of combinations of genomic functional elements
- ***n_core***: int
  - number of cores to use for the analysis

Example:
```python
>> identify_functional_states(chromhmm_path ="./ChromHMM/", number_of_states = 15, n_core = 20)
```

---

```combsafe.show_emission_graph(custom_palette=colors)```<br/>
Show emission parameters heatmap of genome functional states combination. <br/>

Parameters: 
- ***custom_palette***: list of exadecimals, default=None
  - add a list of customized colors in hexadecimal form to be assigned to the functional states. 

Example:
```python
>> colors = ['#c9f9ff', '#e6beff', '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231','#911eb4', '#bcf60c', '#f032e6', '#fffac8', '#fabebe', '#9a6324', '#46f0f0', '#008080']
>> show_emission_graph(custom_palette=colors)
```
![alt text](https://drive.google.com/uc?export=download&id=1Kk_vOm5wz_ski-fLvTxB48dhhu9TXcNY)

---

### Load dataframe of functional states
```combsafe.load_states_dataframe()```<br/>

Return:
- ***functional_states_dataframe***: dataframe
  - dataframe of functional states

Example:
```python
>> functional_states_df = load_functional_states_dataframe()
```

---

### Single-gene analysis

```combsafe.single_gene_analysis(functional_states_dataframe, path_to_gene_list_file)```<br/>
Given a list of gene symbols in a textual file, the heatmap of the functional states of the related genomic regions is shown. <br/>

Parameters:
- ***functional_states_dataframe***: dataframe
  - dataframe of functional states
- ***path_to_gene_list_file***: path
  - path to the gene list file

Example:
```python
>> single_gene_analysis(functional_states_df, "path_to_gene_list/gene_list.txt")
```
![alt text](https://drive.google.com/uc?export=download&id=1zkj4DhgfR36UiIAM99ohF1byaXojKRNU)

---


### Genome-wide analysis
```combsafe.genome_reduction(functional_states_dataframe)```<br/>
Reduce the initial functional state dataframe to visualize the functional states of the various semantic annotations in the form of a heatmap. <br/>
NB:  the proportions among the functional states are maintained as in the previous dataframe of functional states. <br/>

Parameters:
- ***functional_states_dataframe***: dataframe
  - dataframe of functional states

Return:
- ***reducted_dataframe***: dataframe
  - reducted dataframe of functional states
 
Example:
```python
>> reducted_df = genome_reduction(functional_states_df)
```


```combsafe.data_driven_heatmap(reducted_dataframe)```<br/>
Show a genome-wide heatmap with the most significant clusters of genomic regions based on their patterns of functional states. <br/>

Parameters: 
- ***reducted_dataframe***: dataframe
  - reducted dataframe of functional states

Return: 
- ***cluster_indices***: array
  - array of cluster integer labels for each data sample

Example:
```python
>> genome_wide_heatmap = data_driven_heatmap(reducted_df)
```
![alt text](https://drive.google.com/uc?export=download&id=1jbyS_WY54SfJtCWQhw9tpiYW8vC2QJ_Q)

---


```combsafe.gene_ontology_enrichment_analysis(cluster_indices, reducted_dataframe, significance_cut_off)```<br/>
Show a genome-wide heatmap with the most significant clusters of genomic regions based on their patterns of functional states. <br/>

Parameters:
- ***cluster_indices***: array
  - array of cluster integer labels for each data sample
- ***reducted_dataframe***: dataframe
  - reducted dataframe of functional states
- ***significance_cut_off***: int
  - threshold for p-value enrichment analysis, default = 0.05

Example:
```python
>> gene_ontology_enrichment_analysis(genome_wide_heatmap, reducted_df, 0.05)
```
---

