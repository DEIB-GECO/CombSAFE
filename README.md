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
3. Activate the CombSAFE environment with the command: ```conda activate CombSAFE``` ON Linux and macOS. On Windows systems digit  ```activate CombSAFE```
4. Run the ```notebook/Functional_states_analysis.ipynb```

NB: The `pyGMQL` package additionally requires Java. Please follow the installation procedure [here](https://pygmql.readthedocs.io/en/latest/installation.html). <br/>
NB2: The `PyEnsembl` package additionally requires Ensembl data. Please follow the installation procedure [here](https://pypi.org/project/pyensembl/). <br/>
NB3: For Windows users, Visual Studio v.14 or higher is required. Please follow the installation procedure [here](https://visualstudio.microsoft.com/visual-cpp-build-tools/)

## Cookbook
In the following, we show how to call the functions implemented to easily perform the different steps of our ```CombSAFE``` computational method, providing example resuls for some of them. 

### Generate input dataset from raw  data
```combsafe.create_dataset(sample_list_path, organism, threads=4, from_GEO=False)```<br/>
Run a ChIP-seq peak calling pipeline from input raw data. <br/><br/>
For single-end reads Input files must be structured as follows: <br/>

```
Input_folder/
|-- Raw_Reads/
|   |-- 1.rawreads.fastq
|   |-- 2.rawreads.fastq
|   |-- 3.rawreads.fastq
|   |-- 4.rawreads.fastq
|   |-- 5.rawreads.fastq
|   |-- 6.rawreads.fastq
|   |-- 7.rawreads.fastq
|   |-- ...
|-- Textual_file.txt
```
- `Raw_Reads`  a folder containing raw reads in fastq format
- `Textual_file.txt` a text file containing the following information:
  - `file`, filename of the corresponding raw reads file in the Raw_Reads folder
  - `factor`, transcription factor or histone mark used for the analysis
  - `description`, all available iformations of the biological source from which to extract terms for semantic annotations.

E.g.:

| file             | factor   | description                                                  |
| :--------------- | :------- | :----------------------------------------------------------- |
| 1.rawreads.fastq | CTCF     | low passage primary melanoma cultures                        |
| 2.rawreads.fastq | H3K4me3  | Bone marrow mononuclear cells                                |
| 3.rawreads.fastq | MYC      | human primary monocytes isolated from PBMC of healthy donors |
| ...              | ...      | ...                                                          |	

For paired-end reads Input files must be structured as follows: <br/> 

```
Input_folder/
|-- Raw_Reads/
|   |-- 1.forward_reads.fastq
|   |-- 1.reverse_reads.fastq
|   |-- 2.forward_reads.fastq
|   |-- 2.reverse_reads.fastq
|   |-- 3.forward_reads.fastq
|   |-- 3.reverse_reads.fastq
|   |-- ...
|-- Textual_file.txt
```
- `Raw_Reads`  a folder containing raw reads in fastq format
- `Textual_file.txt` a text file containing the following information:
  - `file_1`, filename of the corresponding forward raw reads file in the Raw_Reads folder
  - `file_2`, filename of the corresponding reverse raw reads file in the Raw_Reads folder
  - `factor`, transcription factor or histone mark used for the analysis
  - `description`, all available informations of the biological source from which to extract terms for semantic annotations.  

E.g.:

| file_1                | file_2                |factor    | description                                                       |
| :-------------------- | :---------------------|:-------- | :---------------------------------------------------------------- |
| 1.forward_reads.fastq | 1.reverse_reads.fastq | H3K4me1  | Human embryonic stem cells received from the John Doe laboratory  |
| 2.forward_reads.fastq | 2.reverse_reads.fastq | H3K4me3  | Nuclei derived from crude preps of adipose tissue                 |
| 3.forward_reads.fastq | 3.reverse_reads.fastq | H3K27me3 | Monocyte-derived macrophage                                       |
| ...                   | ...                   |...       | ...                                                               |	


If you want to start a functional state analysis on GEO experiments, set the ```from_GEO``` label as True. In that scenario, input files must be structured as follows: <br/>

```
Input_folder/
|-- Textual_file.txt
```
- `Textual_file.txt` a text file containing the following information:
  - `sample_id`, Id of the samples on GEO
  - `factor`, transcription factor or histone mark used for the analysis

E.g., 

| sample_id | factor   | 
| :---------| :------- | 
| GSM648494 | H3K4me1  | 
| GSM648495 | H3K4me3  | 
| GSM575280 | H3K27me3 | 
| ...       | ...      |

Parameters: 
- ***path***: path object or file-like object   
  - input path folder 
- ***organism***: string
  - reference genome assembly (e.g., "hg19", "hg38", "mm10", "mm39", "rn7", "danrer11", "dm6", "ce11", etc...)
- ***threads***: int, default 4
  - number of threads for the ChIp-seq pipeline.
- ***from_GEO***: bool, default False
  - if True, CombSAFE downloads raw reads and metadata from the GEO web page of the input GEO GSM Ids

Example:
```python
>> dataset = create_dataset("./Input_folder/", assembly="hg38", threads=20, from_GEO=False)
```

---

### Load input dataset
```combsafe.load_dataset(path, assembly, from_GEO=False)```<br/>
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
  - `sample_id`, univolcal id for each sample
  - `factor`, Transcription Factor or Histone Mark used for the analysis
  - `file`, filename of the corresponding ChIp-seq file in the Chip_Files folder
  - `description`, all available informations of the biological source from which to extract terms for semantic annotations. 

E.g.:

| sample_id    | factor   | file             | description                                                  |
| :----------  | :------- | :----------------| :----------------------------------------------------------- |
| 1            | CTCF     | 1.narrowPeak.bed | low passage primary melanoma cultures                        |
| 2            | H3K4me3  | 2.narrowPeak.bed | Bone marrow mononuclear cells                                |
| 3            | MYC      | 3.narrowPeak.bed | human primary monocytes isolated from PBMC of healthy donors |
|              | ...              | ...      | ...                                                          |


If your dataset is generate from GEO samples and you want to get the description from the GSM GEO webpage, set the ```from_GEO``` label as True. In that scenario, Textual_file.txt must be structured as follows: <br/>

- `Textual_file.txt` a text file containing the following information:
  - `sample_id`, Id of the samples on GEO
  - `factor`, Transcription Factor or Histone Mark used for the analysis
  - `file`, filename of the corresponding ChIp-seq file in the Chip_Files folder


| sample_id | factor   | file             |
| :---------| :------- | :--------------- |
| GSM648494 | H3K4me1  | 1.narrowPeak.bed |
| GSM648495 | H3K4me3  | 2.broadPeak.bed  |
| GSM575280 | H3K27me3 | 3.narrowPeak.bed |
| ...       | ...      | ...              |	

Parameters: 
- ***path***: path object or file-like object   
  - input path folder 
- ***organism***: string
  - reference genome assembly (e.g., "hg19", "hg38", "mm10", "mm39", "rn7", "danrer11", "dm6", "ce11", etc...)
- ***from_GEO***: bool, default False
  - if True, CombSAFE downloads raw reads and metadata from the GEO web page of the input GEO GSM Ids

Example:
```python
>> input_path = import_path("./Input_folder/", assembly="hg38", from_GEO=True)
```

---

### Generate semantic annotations
```combsafe.generate_semantic_annotations(dataset, ontology_1, ontology_2, disease = False, encode_convert=False)```<br/>
Generate semantic annotations about tissue and disease types from the input dataset.<br/>

Parameters:
- ***dataset*** dataset object
  - imported dataset object
- ***ontology_1***: str
  - url address of the first selected ontology
- ***ontology_2***: str
  - url address of the second selected ontology
- ***disease*** bool, default False
  - set True if one of the selected ontologies is involved in disease concepts
- ***encode_convert*** bool, default False
  - if true, encode IDs are searched to be converted to GSM

Returns: 
  - ***semantic_dataframe***
    - dataset of biological conditions from the input metadata

Example:
```python
>> ontology_1 = "https://raw.githubusercontent.com/obophenotype/cell-ontology/master/cl-basic.obo"
>> ontology_2 = "https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.obo"
>> semantic_df = generate_semantic_annotations(dataset, ontology_1, ontology_2, disease = True, encode_convert=False)
```

---

### Data analysis
```combsafe.plot_factor_freq(semantic_dataframe, n)```<br/>
Vertical barplot of the factor frequency in the input dataset.<br/>

Parameters: 
- ***semantic_dataframe***: dataframe
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

```combsafe.add_custom_tracks(tracks_label, path_to_custom_tracks, index)```<br/>
Add custom tracks of static genomic elements to the analysis (e.g., CpG islands). <br/>

Parameters: 
- ***tracks_label***: string
  - name of the custom tracks to add for the analysis
- ***path_to_custom_tracks***: path
  - UCSC path for downloading custom tracks

Example:
```python
>> add_custom_tracks("CpG_Islands", "./input_files/cpgIslandExt.txt")
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

Returns: 
  - ***functional_states_dataframe***
    - dataset of functional states for each biological conditions


```python
>> functional_states_df = identify_functional_states(chromhmm_path ="./ChromHMM/", number_of_states = 15, n_core = 20)
```

Alternatively, it is possible to load in house segmentated files from an other HMM segmentation tool and jump to the next step

```combsafe.load_custom_segments(input_segment_dir, num_states)```<br/>
load functional states files from input path. <br/>

Parameters:
- ***input_segment_dir***: path
  - path to the segmentated file folder.
- ***number_of_states***: int
  - number of combinations of genomic functional elements

Returns: 
  - ***functional_states_dataframe***
    - dataset of functional states for each biological conditions 
 

Example:
```python
>> functional_states_df = load_custom_segments(input_segment_dir ="./Input_folder/in_house_segmentated/", num_states=15)
```

Input files must be structured as follows: <br/>

```
Input_folder/
|-- Segmentated_Files/
|   |-- 1.semantic_annotation_15_segments.bed
|   |-- 2.semantic_annotation_15_segments.bed
|   |-- 3.semantic_annotation_15_segments.bed
|   |-- 4.semantic_annotation_15_segments.bed
|   |-- 5.semantic_annotation_15_segments.bed
|   |-- 6.semantic_annotation_15_segments.bed
|   |-- 7.semantic_annotation_15_segments.bed
|   |-- ...
|-- emissions.txt
```

- `Segmentated_Files`  a folder containing raw reads in fastq format
- `emissions.txt` a text file structured as follows:

| State     | H3K4me3     | POLR2A      |MYC          | H3K27me3    |
| :---------| :---------- | :---------- | :---------- | :---------- |
| 1         | 0.093326    | 0.457892    | 0.143540    | 0.924645    |
| 2         | 0.793153    | 0.658634    | 0.972344    | 0.487613    |
| 3         | 0.940996    | 0.000234    | 0.243758    | 0.187461    |
| 4         | 0.143540    | 0.763471    | 0.872346    | 0.104765    |
| ...       | ...         | ...         |...          | ...         |	

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

###  Distance Metric Heatmap

```combsafe.show_distance_matrix()```<br/>
Show distance matrix heatmap about functional states generated from the emission parameters file of an HMM model. <br/>

Example:
```python
>> show_distance_matrix()
```

![alt text](https://drive.google.com/uc?export=download&id=1TBzP1xohnIjRwjco5OTeLCEsDI47Xnz-)

---

### Single-gene analysis

```combsafe.single_gene_analysis(functional_states_dataframe, path_to_gene_list_file, distance_metric )```<br/>
Given a list of gene symbols in a textual file, the heatmap of the functional states of the related genomic regions is shown. <br/>

Parameters:
- ***functional_states_dataframe***: dataframe
  - dataframe of functional states
- ***path_to_gene_list_file***: path
  - path to the gene list file
- ***distance_metric***: function/string
  - distance metric to use for the data. CombSAFE auto-generate from the input ```emission.txt``` file a special metric ```functional_state_distance``` to weight distances among functional states according to the their function. Alternatively, ```hamming``` distance can be selected


Example:
```python
>> single_gene_analysis(functional_states_df, "path_to_gene_list/gene_list.txt", distance_metric = funtional_states_distance)
>> single_gene_analysis(functional_states_df, "path_to_gene_list/gene_list.txt", distance_metric = "hamming")
```
![alt text](https://drive.google.com/uc?export=download&id=1zkj4DhgfR36UiIAM99ohF1byaXojKRNU)

---

###  PCA Analysis

```combsafe.show_distance_matrix()```<br/>
Shows PCA heatmap among semantic annotation for selected components. <br/>

Parameters:
- ***functional_states_dataframe***: dataframe
  - dataframe of functional states
- ***number_of_components***: int
  -  number of components for the principal component analysis
 
Return:
- ***loadings***: dataframe
  - PCA loadings

Example:
```python
>> pca_analysis(functional_states_df, 10)
```

![alt text](https://drive.google.com/uc?export=download&id=1Ph4TnBsxYnlgUUGZ2gYWLNqU7M96QB6L)

---

### Genome-wide analysis
```combsafe.genome_reduction(functional_states_dataframe, threshold)```<br/>
Reduce the initial functional state dataframe to visualize the functional states of the various semantic annotations in the form of a heatmap. <br/>
NB:  the proportions among the functional states are maintained as in the previous dataframe of functional states. <br/>

Parameters:
- ***functional_states_dataframe***: dataframe
  - dataframe of functional states
- ***threshold***: int
  -  value adopted to avoid the over-representation of states that are specific for a small subset of the conditions or for one condition only. Select 100 for full representation
 
Return:
- ***reducted_dataframe***: dataframe
  - reducted dataframe of functional states
 
Example:
```python
>> reducted_df = genome_reduction(functional_states_df, threshold=90)
```


```combsafe.data_driven_heatmap(reducted_df, distance_metric, min_clust_size, min_sampl)```<br/>
Show a genome-wide heatmap with the most significant clusters of genomic regions based on their patterns of functional states. <br/>

Parameters: 
- ***reducted_df***: dataframe
  - reducted dataframe of functional states
- ***distance_metric***: function/string
  - distance metric to use for the data. CombSAFE auto-generate from the input ```emission.txt``` file a special metric ```functional_state_distance``` to weight distances among functional states according to the their function. Alternatively, ```hamming``` distance can be selected
- ***min_clust_size***: int
  - minimum number of clusters accetpted for the analysis
- ***min_sampl***: int
  - minimum number of samples per cluster accepted for the analysis

Return: 
- ***clustered_dataframe***: dataframe
  - dataframe of functional states ordered according to the cluster parameters

Example:
```python
>> clustered_heatmap = data_driven_heatmap(reducted_df, functional_states_distance, min_clust_size=10, min_sampl=2)
```
![alt text](https://drive.google.com/uc?export=download&id=1jbyS_WY54SfJtCWQhw9tpiYW8vC2QJ_Q)

---


```combsafe.gene_ontology_enrichment_analysis(clustered_dataframe, distance_metric, goea_tool)```<br/>
Show a genome-wide heatmap with the most significant clusters of genomic regions based on their patterns of functional states. <br/>

Parameters:
- ***clustered_dataframe***: dataframe
  - adataframe of functional states ordered according to the cluster parameters
- ***distance_metric***: function/string
  - distance metric to use for the data. CombSAFE auto-generate from the input ```emission.txt``` file a special metric ```functional_state_distance``` to weight distances among functional states according to the their function. Alternatively, ```hamming``` distance can be selected
- ***goea_tool***: str, "great" or "goatool"
  - tool for gene ontology enrichment analysis 

Example:
```python
>> gene_ontology_enrichment_analysis(clustered_heatmap, goea_tool = "great", distance_metric=functional_states_distance)
```
---


