combsafe_output = "./CombSAFE_output/"
GMQL_output_dir = "./"
n_states = 0
n_states_model = GMQL_output_dir + str(n_states) + "_state_model/"
chromHMM_output = n_states_model + "ChromHMM_output/"
graphs_path = n_states_model + "graphs/"
gmql_output = "./"
destination_path = "./"
custom_tracks = False
custom_description = False
custom_segments = False
custom_name = ""
organism=""
custom_path = "./"
emissions_path = "/"
file_path = "./"
clusters_path = "/"
bed_folder_path = "./"
nrm_df = False


color_list = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']



import os
import re
import gzip
import mygene
import psutil
import shutil
import urllib
import sqlite3
import logging
import warnings
import requests
import random
import hdbscan
import itertools
import subprocess
import gmql as gl
import numpy as np
import pandas as pd
from lxml import etree
import seaborn as sns
from PIL import Image
from scipy import stats
from pathlib import Path
from rpy2 import robjects
from urllib import request
import matplotlib.cm as cm
from tabulate import tabulate
import matplotlib.pyplot as plt
from fastcluster import linkage
from scipy.spatial import distance
from rpy2.robjects import pandas2ri
from pyensembl import EnsemblRelease
from scipy.stats import fisher_exact
from goatools.obo_parser import GODag
from sklearn.decomposition import PCA
from scipy.cluster import hierarchy as hc
from scipy.spatial import distance as ssd
from rpy2.robjects.packages import importr
from scipy.cluster.hierarchy import dendrogram
from sklearn.metrics import pairwise_distances
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj

from Bio import Entrez
#install the genome you need via pyensembl install --release 97 --species human where 97 can be replaced with any version you want
data = EnsemblRelease(97)


rpy2_logger.setLevel(logging.ERROR)
warnings.filterwarnings('ignore')
pandas2ri.activate()

def load_dataset(sample_list_path, assembly, from_GEO=False):
    
    global file_path
    global bed_folder_path
    global organism
    global custom_description
    
    if not from_GEO:
        custom_description = True

    file_path = path + next(os.walk(path))[2][0]
    organism = assembly
    bed_folder_path = path + next(x for x in next(os.walk(path))[1] if "." not in x)
    
    return(file_path, bed_folder_path)

def download_reference_genome(org, r_path):
    url = "http://hgdownload.cse.ucsc.edu/goldenPath/"+ org + "/bigZips/" + org + ".fa.gz"
    gz_file = url.split("/")[-1]
    fasta_file = r_path + gz_file.split(".gz")[0]
    base = fasta_file.split("/")[2].split(".fa")[0]
    r1 = request.urlretrieve(url=url, filename=gz_file)
    fasta_file = r_path + gz_file.split(".gz")[0]
    with gzip.open(gz_file, 'rb') as f_in:
        with open(fasta_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(gz_file)
    
    
def build_reference_index(ref_org):
    genomes = ["hg19", "hg38", "mm10", "mm39", "rn6", "rn7", "danrer10", "danrer11", "dm3", "dm6", "ce10", "ce11"]
    if ref_org not in genomes:
        raise ValueError("Wrong genome name")
    else:
        ref_path = "./processed_data/reference_genomes/" + ref_org + "/"
        subprocess.call(" ".join(["mkdir", "-p", ref_path]), shell=True)
        if not any(os.scandir(ref_path)):
            download_reference_genome(org = ref_org, r_path = ref_path)
            build_index = " ".join(["bowtie2-build", ref_path + ref_org + ".fa", ref_org]) 
            subprocess.call(build_index, shell=True)
            for file in os.listdir("./"):
                if file.endswith(".bt2"):
                    shutil.move(os.path.abspath(file), ref_path)

def generate_dataset(sample_list_path, organism, threads=4, from_GEO=False):
    
    n_threads=str(threads)
    build_reference_index(organism)
    samples = pd.read_csv(sample_list_path, "\t")
    broad = ["H3F3A", "H3K27ME3", "H3K36ME3", "H3K4ME1", "H3K79ME2", "H3K79ME3", "H3K9ME1", "H3K9ME2", "H4K20ME1"]
    
    out_dir = f"./processed_data/called_peaks/{organism}/peaks/"
    
    ref_path = f"./processed_data/reference_genomes/{organism}/"
    fasta_folder = f"{out_dir}/fastq_files/" 
    fastqc_folder = f"{out_dir}/fastqc/"
    bowtie_results = f"{out_dir}/bowtie2/"
    intermediate_bams = f"{out_dir}/bowtie2/intermediate_bams/"
    macs2_out = f"{out_dir}/macs2_out/"
    
    subprocess.call(" ".join(["mkdir", "-p", out_dir]), shell=True)
    file = open(out_dir.split("/peaks")[0] + "/" + organism + ".txt", "w")
    file.write("SampleID" + "\t" + "Factor" + "\t" + "Filename" + "\n")
    
    for i in range(samples.shape[0]):
        subprocess.call(" ".join(["mkdir", "-p", fasta_folder]), shell=True)
        subprocess.call(" ".join(["mkdir", "-p", fastqc_folder]), shell=True)
        subprocess.call(" ".join(["mkdir", "-p", intermediate_bams]), shell=True)
        
        if not from_GEO:
            global custom_description
            custom_description=True
            sample_name = samples.iloc[i,0].split(".fastq")[0]
        else:
            sample_name = samples.iloc[i,0]
        
        align_sam = f"{out_dir}/bowtie2/{sample_name}_m_unsorted.sam"
        align_bam = f"{out_dir}/bowtie2/{sample_name}_m_unsorted.bam"
        align_n_sorted = f"{out_dir}/bowtie2/{sample_name}_m_name_sorted.bam"
        align_fixmate = f"{out_dir}/bowtie2/{sample_name}_m_fixmate_sorted.bam"
        align_filtered = f"{out_dir}/bowtie2/{sample_name}_m_filtered_sorted.bam"
        align_sorted = f"{out_dir}/bowtie2/{sample_name}_sorted.bam"
                       
        if not from_GEO:
            input_fasta_dir = "/".join(sample_list_path.split("/")[0:-1]) + "/" + [f for f in next(os.walk("/".join(sample_list_path.split("/")[0:-1])))[1] if not f.startswith('.')][0] + "/"
            factor_name = samples.iloc[i,2]
            
            if (isinstance(samples.iloc[i, 1], float)):
                layout = "SINGLE"
                shutil.copy(input_fasta_dir + sample_name + ".fastq", fasta_folder + sample_name + ".fastq")
            else:
                layout = "PAIRED" 
                fastq_1 = fasta_folder + samples.iloc[i,0] 
                fastq_2 = fasta_folder + samples.iloc[i,1]
                shutil.copy(input_fasta_dir + sample_name + ".fastq", fastq_1)                
                shutil.copy(input_fasta_dir + sample_name + ".fastq", fastq_2)                
                
        else:
            factor_name = samples.iloc[i,1]
            info_query = ["esearch", "-db", "sra", "-query", sample_name, "|", "efetch", "-format", "runinfo"]
            proc = subprocess.Popen(" ".join(info_query), shell=True, stdout=subprocess.PIPE)
            serviceList  = proc.communicate()[0].decode('utf-8')
            
            srr_id=[]
            values = serviceList.rstrip().split('\n')
            for runs in range(1, len(values)):
                for attributes in range(0, len(values[0].split(","))):
                    if values[0].split(",")[attributes] == 'Run':
                        srr_id.append(values[runs].split(",")[attributes])
                    if values[0].split(",")[attributes] == "LibraryLayout":
                        layout = values[1].split(",")[attributes]
                    if values[0].split(",")[attributes] == "ScientificName":
                        organism_type = values[1].split(",")[attributes]
        
        factor_type = "--broad" if factor_name.upper() in broad else ""
                    
        if layout.upper() == 'SINGLE':
            #print("LAYOUT: SINGLE")
            
            fastq_name = fasta_folder + sample_name + ".fastq"            
                  
            if from_GEO:
                fastq_name = fasta_folder + sample_name + ".fastq"
                for srr in srr_id:
                    download_query = " ".join(["parallel-fastq-dump", "--sra-id", srr, "--threads", n_threads, "--outdir", fasta_folder])
                    subprocess.call(download_query, shell=True)        

                merge_files = " ".join(["cat", fasta_folder + "*.fastq", ">",  fastq_name])
                subprocess.call(merge_files, shell=True)
                
            fastqc_query = " ".join(["fastqc", fastq_name, "--extract"])
            subprocess.call(fastqc_query, shell=True)

            mv_fastqc_query = " ".join(['mv', fasta_folder + '*_fastqc*', fastqc_folder])
            subprocess.call(mv_fastqc_query, shell=True)

            fastqc_df = pd.read_csv(fastqc_folder + sample_name + "_fastqc/fastqc_data.txt", nrows=10, sep="\t")
            encoding = str(fastqc_df[fastqc_df["##FastQC"]=="Encoding"].iloc[:, 1]).split("\n")[0].split("  ")[-1]
            quality = "" if "Sanger" in encoding else "--phred64"

            bowtie2 = " ".join(['bowtie2', '-p', n_threads, quality, '-q', '-D', '20', '-R', '3', '-N', '1', '-L', '20', '-x', ref_path + organism, '-U', fastq_name, '-S', align_sam, '2> ./' + out_dir + '/log_file.txt'])
            subprocess.call(bowtie2, shell=True)

            drop_srr = " ".join(["rm", fasta_folder + "*SRR*"]) 
            subprocess.call(drop_srr, shell=True)

            bam_conversion = " ".join(["samtools", "view", "-h", "-S", "-b", "-@", n_threads, "-o", align_bam, align_sam])
            subprocess.call(bam_conversion, shell=True)

            bam_name_sorting = " ".join(["samtools", "sort", "-n", "-@", n_threads, "-o", align_n_sorted, align_bam])
            subprocess.call(bam_name_sorting, shell=True)

            fixmate = " ".join(["samtools", "fixmate", "-@", n_threads, align_n_sorted, align_fixmate])
            subprocess.call(fixmate, shell=True)

            bam_sorting = " ".join(["samtools", "sort", "-@", n_threads, "-o", align_sorted, align_fixmate])
            subprocess.call(bam_sorting, shell=True)

            filtered = " ".join(["samtools", "markdup","-r", "-@", n_threads, align_sorted, align_filtered])
            subprocess.call(filtered, shell=True)

            index_creation  = " ".join(["samtools", "index", "-@", n_threads, align_sorted])
            subprocess.call(index_creation, shell=True)

            move_intermediate  = " ".join(["mv", bowtie_results + "*_m_*", intermediate_bams])
            subprocess.call(move_intermediate, shell=True)

            macs2 = " ".join(["macs2", "callpeak", "-q", "0.01", "--keep-dup 1", "--extsize=150", "--nomodel", "-f", "BAM", "-g", "2.9e+9", factor_type, "-t", align_sorted, "-n", sample_name, "--outdir", macs2_out])
            subprocess.call(macs2, shell=True)

        if layout.upper() == 'PAIRED':
            #print("LAYOUT: PAIRED")
                               
            if from_GEO:
                fastq_1 = fasta_folder + sample_name + "_1" + ".fastq"
                fastq_2 = fasta_folder + sample_name + "_2" + ".fastq"
                for srr in srr_id:
                    download_query = ["parallel-fastq-dump", "-F", "--split-files", "--sra-id", srr, "--threads", n_threads, "--outdir", fasta_folder]
                    subprocess.call(" ".join(download_query), shell=True)

                r1 = " ".join(["cat", fasta_folder + "*_1.fastq", ">",  fastq_1])
                subprocess.call(r1, shell=True) 

                r2 = " ".join(["cat", fasta_folder + "*_2.fastq", ">",  fastq_2])
                subprocess.call(r2, shell=True)
                
            if not from_GEO:
                fastqc_query = " ".join(["fastqc", fasta_folder + "*.fastq*", "--extract"])            
            else:
                fastqc_query = " ".join(["fastqc", fasta_folder + "*GSM*", "--extract"])

            
            subprocess.call(fastqc_query, shell=True)

            mv_fastqc_query = " ".join(['mv', fasta_folder + '*_fastqc*', fastqc_folder])
            subprocess.call(mv_fastqc_query, shell=True)

            if not from_GEO:
                fastqc_df = pd.read_csv(fastqc_folder + sample_name + "_fastqc/fastqc_data.txt", nrows=10, sep="\t")            
            else:
                fastqc_df = pd.read_csv(fastqc_folder + sample_name + "_1_fastqc/fastqc_data.txt", nrows=10, sep="\t")
            
            encoding = str(fastqc_df[fastqc_df["##FastQC"]=="Encoding"].iloc[:, 1]).split("\n")[0].split("  ")[-1]
            quality = "" if "Sanger" in encoding else "--phred64"

            bowtie2 = " ".join(['bowtie2', '-p', n_threads, quality, '-q', '-D', '20', '-R', '3', '-N', '1', '-L', '20', '-x', ref_path + organism, '-1', fastq_1, '-2', fastq_2, '-S', align_sam, '2> ./' + out_dir + '/log_file.txt'])
            subprocess.call(bowtie2, shell=True)

            drop_srr = " ".join(["rm", fasta_folder + "*SRR*"]) 
            subprocess.call(drop_srr, shell=True)

            bam_conversion = " ".join(["samtools", "view", "-h", "-S", "-b", "-@", n_threads, "-o", align_bam, align_sam])
            subprocess.call(bam_conversion, shell=True)

            fixmate = " ".join(["samtools", "fixmate", "-@", n_threads, "-m", align_bam, align_fixmate])
            subprocess.call(fixmate, shell=True)

            bam_sorting = " ".join(["samtools", "sort", "-@", n_threads, "-o", align_sorted, align_fixmate])
            subprocess.call(bam_sorting, shell=True)

            filtered = " ".join(["samtools", "markdup","-r", "-@", n_threads, "-s", align_sorted, align_filtered])
            subprocess.call(filtered, shell=True)

            index_creation  = " ".join(["samtools", "index", "-@", n_threads, align_sorted])
            subprocess.call(index_creation, shell=True)

            move_intermediate  = " ".join(["mv", bowtie_results + "*_m_*", intermediate_bams])
            subprocess.call(move_intermediate, shell=True)

            macs2 = " ".join(["macs2", "callpeak", "-q", "0.01","-c", "--keep-dup 1", "--extsize=150", "--nomodel", "-f", "BAMPE", "-g", "2.9e+9", factor_type, "-t", align_sorted, "-n", sample_name, "--outdir", macs2_out])
            subprocess.call(macs2, shell=True) 

        bed_file_type = "broad" if factor_type.split("-")[-1] == 'broad' else "narrow"
        bed_file = macs2_out + sample_name + "_peaks." + bed_file_type +  "Peak"
        
        file.write(sample_name + "\t" + factor_name + "\t" + bed_file.split("/")[-1] + "\n")
        
        shutil.move(bed_file, out_dir + bed_file.split("/")[-1])
        
        for files in os.listdir(out_dir):
            if not (files.endswith("narrowPeak") or files.endswith("broadPeak")):
                if os.path.isdir(out_dir + files):
                    shutil.rmtree(out_dir + files)
                else:
                    os.remove(out_dir + files)     
    file.close()
    
    load_dataset(f"./processed_data/called_peaks/{organism}/", organism, from_GEO)

def collapse_annotations(a):
    if a is not None:
        terms = a.split(", ") 
        for term in a.split(", "):
            for term_2 in terms:
                if (term_2 in term) and term_2 != term:
                    terms.remove(term_2)
        if len(terms) == 1:
            return(terms[0]).replace("adenocarcinoma", "cancer").replace("carcinoma", "cancer")
        else:
            return(", ".join(terms).replace("adenocarcinoma", "cancer").replace("carcinoma", "cancer"))
    
        

def download_ontology(url, file_name):
    r = requests.get(url)
    with open(combsafe_output + 'annotations/ontologies/' + file_name, 'wb') as f:
        f.write(r.content)
 
'''
ontology_list = [
"po", 
"obi", #285/500 -> https://raw.githubusercontent.com/obi-ontology/obi/v2021-04-06/views/obi.obo
"pato",
"xao",
"zfa", 
"aism",
"bto", #170/500 -> https://raw.githubusercontent.com/BRENDA-Enzymes/BTO/master/bto.obo #8
"uberon", #105/500 -> https://raw.githubusercontent.com/obophenotype/uberon/master/basic.obo
"cl", #278/500 -> https://raw.githubusercontent.com/obophenotype/cell-ontology/master/cl.obo
"cl-basic", #225/500 -> https://raw.githubusercontent.com/obophenotype/cell-ontology/master/cl-basic.obo
"fbcv"
"doid" #89/500 ->https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/main/src/ontology/doid.obo
]
'''

def generate_semantic_annotations(dataset, ontology_1, ontology_2, disease = False, encode_convert=False):
    
    robjects.r('''

    if (!requireNamespace("BiocManager", quietly = TRUE)){
        install.packages("BiocManager", repos='http://cran.rstudio.com/')
    }    

    BiocManager::install(c("GEOmetadb","Onassis","OnassisJavaLibs"), update=FALSE)
    options(java.parameters = "-Xmx8g" )

    ''')
    
    geometadb = importr('GEOmetadb')
    onassis = importr('Onassis')
    utils = importr('utils')
    base= importr('base')
    dbi = importr('DBI')
    
    sql_database = Path(combsafe_output + 'annotations/GEOmetadb/GEOmetadb.sqlite')
    
    if not sql_database.is_file():
        os.makedirs(combsafe_output + 'annotations/GEOmetadb/')
        geometadb.getSQLiteFile(combsafe_output + 'annotations/GEOmetadb')

    if not os.path.exists(combsafe_output + 'annotations/ontologies/dictionary'):
        os.makedirs(combsafe_output + 'annotations/ontologies/dictionary')

    if custom_description:

        gsm_metadata_f = utils.read_csv(file_path, header = True, sep="\t", quote="")

    else:

        gsm_list = get_geoids_list(separator="\t", enc_convert=encode_convert)
        geo_con = onassis.connectToGEODB(str(sql_database))
        gsm_metadata = dbi.dbGetQuery(geo_con, base.paste0("select * from gsm where gsm in ('", base.paste(gsm_list, collapse="','"), "')"))
        gsm_metadata_f = gsm_metadata.loc[:, ['gsm', 'characteristics_ch1', 'source_name_ch1']]
    
    first_annotation = onassis.annotate(gsm_metadata_f, 'OBO', ontology_1, Stemmer = 'BIOLEMMATIZER') 
    first_annotation_f = onassis.filterconcepts(first_annotation, base.c('cell', 'donor', 'protein', 'tissue','function', 'immunoglobulin complex', 'chromatin', 'Homo sapiens', 'cell', 'organism', 'circulating', 'heterochromatin', 'quality', 'molecule', 'normal', 'group', 'base', 'inhibitor', 'acid', 'time' ,'signaling', 'localization', 'system', 'gene', 'binding', 'affinity', 'chromosome', 'structure', 'Mus musculus', 'Bos taurus', 'oxygen atom', 'atomic nucleus', 'nucleus', 'dforkhead box protein P3', 'mature', 'cell cycle', 'catalytic activity', 'Drosophila melanogaster', 'Drosophila <fruit fly, genus>', 'V(D)J recombination-activating protein 1', 'Mus <mouse', 'genus>', 'Drosophila <fruit fly', 'nutrient', 'cell line cell', 'Homo', 'chromatin', 'organism', 'protein', 'signaling', 'age', 'female', 'sex', 'male', 'undifferentiated', 'size', 'disease', 'viability', 'document', 'cell line', 'cell line lineage', 'cell type'))
    
    first_annotation_f = onassis.sim(first_annotation_f, pairconf='edge_rada_lca')
    collapsed_annotations = onassis.collapse(first_annotation_f, 0.7)

    for file in os.listdir("./"):
        if file.endswith(".a1"):
            os.remove(file)
            
    if disease:
        doid="TRUE"
    else:
        doid="FALSE"
        
    second_annotation = onassis.annotate(gsm_metadata_f, 'OBO', ontology_2, disease=doid)
    second_annotation = onassis.filterconcepts(second_annotation, base.c('disease'))

    for file in os.listdir("./"):
        if file.endswith(".a1"):
            os.remove(file)
            
    semantic_ontologies = onassis.mergeonassis(collapsed_annotations, second_annotation)
    semantic_matrix = onassis.entities(semantic_ontologies)
    
    semantic_matrix['term_name_1'] = semantic_matrix["term_name_1"].apply(lambda x: collapse_annotations(x))
    semantic_matrix['term_name_2'] = semantic_matrix["term_name_2"].apply(lambda x: collapse_annotations(x))
    
    similarity_matrix = onassis.simil(semantic_ontologies)
    s_matrix = pd.DataFrame(similarity_matrix)
    s_matrix.to_csv(combsafe_output + 'annotations/similarity_matrix.txt', sep="\t", index=False)
    semantic_matrix.to_csv(combsafe_output + 'annotations/semantic_matrix.txt', sep="\t", index=False)

    
    for src_filename in os.listdir():
        if src_filename.endswith(".xml"):
            dst_filename = os.path.join(combsafe_output + 'annotations/ontologies/dictionary/', os.path.basename(src_filename))
            shutil.move(src_filename, dst_filename)
        if src_filename.endswith(".obo"):
            dst_filename = os.path.join(combsafe_output + 'annotations/ontologies/', os.path.basename(src_filename))
            shutil.move(src_filename, dst_filename)
    
    info_file = pd.read_csv(dataset[0], "\t")
    ontology_file = pd.read_csv(combsafe_output + "annotations/semantic_matrix.txt", "\t")
    ontology_file_f = ontology_file.loc[1:,['sample_id', 'short_label_1', 'term_name_2']]
    ontology_file_f.columns = ['sample_id', 'Tissue', 'Disease']
    merged_df = info_file.merge(ontology_file_f, on=['sample_id']).fillna("NA")
    merged_df["Semantic_Annotation"] = merged_df['Tissue'].str.replace(r"[\(\[].*?[\)\]]", "", regex=True).str.rstrip().str.lower().replace(" , ", "_", regex=True) + "_000_" + merged_df['Disease'].str.lower().replace(', ','_', regex=True)
    merged_df.to_csv(combsafe_output + "annotations/annotation_matrix.txt", sep="\t", index=False)
    return(merged_df)


def save_dict_to_file(dic):
    f = open(n_states_model + 'dict.txt','w')
    f.write(str(dic))
    f.close()

def load_dict_from_file():
    f = open(n_states_model + 'dict.txt','r')
    data_information=f.read()
    f.close()
    return eval(data_information)


def get_geoids_list(separator="\t", enc_convert=False):
    info_file = pd.read_csv(file_path, sep=separator)
    geo_ids = []
    encode_ids = []
    ids_list = list(set(info_file.GSMID))
    tot_enc = sum('ENC' in s for s in ids_list)
    for ids in ids_list:
        if enc_convert:
            if "ENC" in ids:
                encode_ids.append(ids)
            else:
                geo_ids.append(ids)
                
                geo_ids_converted = []
                for names in range(len(encode_ids)):
                    name = encode_ids[names].split("_")[0]
                    info_query = " ".join(["esearch", "-db", "gds", "-query", name, "|", "efetch", "-format", "runinfo"])
                    proc = subprocess.Popen(info_query, shell=True, stdout=subprocess.PIPE)
                    serviceList  = proc.communicate()[0].decode('utf-8')
                    new_id = serviceList.rstrip().split("Accession")[-1].split("\t")[0][2:]
                    if bool(new_id):
                        geo_ids_converted.append(new_id) 
                    print(str(names) + " on a total of " + str(tot_enc) + " ids succesfully converted", end = "\r")
                    time.sleep(0.5)
        else:
            if "GSM" in ids:
                geo_ids.append(ids) 
   
    if enc_convert:
        print(str(len(geo_ids_converted)) + " on a total of " + str(len(encode_ids)) + " encode ids successfully converted")
        total_geo_ids = list(set(geo_ids + geo_ids_converted))
    else:
    	total_geo_ids = list(set(geo_ids))
    
    return(total_geo_ids)



def plot_factor_freq(semantic_dataframe, n):
    factor_frequency = semantic_dataframe['Factor'].value_counts().rename_axis('Factors').reset_index(name='Frequency')
    top_factors = factor_frequency.head(n)
    new = [tuple(r) for r in top_factors[['Factors', 'Frequency']].values]

    new = sorted(new, key=lambda score: score[1], reverse=True) 
    X, Y = zip(*new)    
    plt.figure(figsize = (15, 3))
    plt.title('Factor frequency')
    plt.ylabel('Number of samples')
    bars = plt.bar(range(len(X)), Y, 0.6, tick_label = X, color="red") 
    plt.xticks(rotation=90)
    hm=[]
    factors=[]
    for i, bar in enumerate(bars):
        if re.findall("^H[1-5]", X[i]):
            bar.set_color("blue")
            hm.append(X[i])
        else:
            factors.append(X[i])
    plt.show()



def generate_fixed_factor_pool(semantic_dataframe, factor_list, n):
    factor_frequency = semantic_dataframe['Factor'].value_counts().rename_axis('Factors').reset_index(name='Frequency')
    f_list = factor_frequency.head(20)['Factors']
    clean_set = pd.Series(list((set(f_list) - set(factor_list))))
    cols = ['factor_list', 'n_semantic_annotations']
    lst = []
    for n, comb in enumerate(itertools.combinations(clean_set, n - len(factor_list))):
        new = factor_list + list(comb)
        selection = semantic_dataframe[(semantic_dataframe.Semantic_Annotation != "lining cell_native cell_secretory cell_000_unknown") & (semantic_dataframe.Semantic_Annotation != "lining cell_native cell_secretory cell_000_breast cancer_cancer")]
        selection_f = selection[selection["Factor"].isin(new)]
        
        f_count = selection_f.loc[1:, ["Semantic_Annotation", "Factor"]].groupby("Semantic_Annotation")['Factor'].apply(lambda x: len(list(np.unique(x)))).to_frame()
        n_samples = f_count[f_count["Factor"] == len(new)].shape[0]
        lst.append([", ".join(new), n_samples])

    df1 = pd.DataFrame(lst, columns=cols)
    final_df = df1.sort_values(by=['n_semantic_annotations'], ascending=False).reset_index(drop=True).head(10)
    final_df.index = final_df.index + 1
    print(tabulate(final_df, headers='keys', tablefmt='psql'))
    

def get_semantic_annotation_list(semantic_dataframe, factor_list):
    v=[]
    selection = semantic_dataframe[semantic_dataframe["Factor"].isin(factor_list)]
    f_count = selection.loc[1:, ["Semantic_Annotation", "Factor"]].groupby("Semantic_Annotation")['Factor'].apply(lambda x: len(list(np.unique(x)))).to_frame()
    n_samples = f_count[f_count["Factor"] == len(factor_list)]
    n_samples_list = list(n_samples.index)
    for i in range(len(n_samples_list)):
        v.append((str(i+1) + " - " + n_samples_list[i]))
    for i in v:
        print(i)

def generate_schema(narrow=False):

    data_format = 'narrow' if narrow == True else 'broad' 
    output = combsafe_output + "gdm_metadata/" + data_format + "/schema.xml"
    
    # create XML 
    root = etree.Element('gmqlSchemaCollection')
    root.set("name", data_format)
    root.set("xmlns", "http://genomic.elet.polimi.it/entities")

    #create gmqlSchema
    gmqlSchema = etree.Element('gmqlSchema')
    gmqlSchema.set("type", "Peak")
    gmqlSchema.set("coordinate_system", "default")
    root.append(gmqlSchema)

    # create fields with text
    char = etree.Element('field')
    char.set("type", "STRING")
    char.text = 'chr'
    gmqlSchema.append(char)

    left = etree.Element('field')
    left.set("type", "LONG")
    left.text = 'left'
    gmqlSchema.append(left)

    right = etree.Element('field')
    right.set("type", "LONG")
    right.text = 'right'
    gmqlSchema.append(right)

    strand = etree.Element('field')
    strand.set("type", "CHAR")
    strand.text = 'strand'
    gmqlSchema.append(strand)

    name = etree.Element('field')
    name.set("type", "STRING")
    name.text = 'name'
    gmqlSchema.append(name)

    score = etree.Element('field')
    score.set("type", "DOUBLE")
    score.text = 'score'
    gmqlSchema.append(score)

    signal = etree.Element('field')
    signal.set("type", "DOUBLE")
    signal.text = 'signal'
    gmqlSchema.append(signal)

    pvalue = etree.Element('field')
    pvalue.set("type", "DOUBLE")
    pvalue.text = 'pvalue'
    gmqlSchema.append(pvalue)

    qvalue = etree.Element('field')
    qvalue.set("type", "DOUBLE")
    qvalue.text = 'qvalue'
    gmqlSchema.append(qvalue)
    
    if narrow:
        peak = etree.Element('field')
        peak.set("type", "DOUBLE")
        peak.text = 'peak'
        gmqlSchema.append(peak)

    with open(output, 'wb') as doc:
        doc.write(etree.tostring(root, pretty_print=True, xml_declaration=True, encoding="UTF-8"))



def generate_gdm_format(dataset):
    
    file = open(combsafe_output + "annotations/annotation_matrix.txt").readlines()

    if not os.path.exists(combsafe_output + 'gdm_metadata/narrow'):
        os.makedirs(combsafe_output + 'gdm_metadata/narrow')

    if not os.path.exists(combsafe_output + 'gdm_metadata/broad'):
        os.makedirs(combsafe_output + 'gdm_metadata/broad')

    sample_id = []
    keys = file[0].rstrip().split("\t")
    index_name = keys.index("File")

    for i in range(1, len(file)):
        attributes = (file[i].rstrip().split("\t"))
        name=attributes[index_name]
        sample_id.append(name)
        metafile=open(combsafe_output + "gdm_metadata/" + name + ".meta", "w")
        for n in range(0, len(keys)):
            metafile.write(keys[n] + "\t" + attributes[n].replace("_000_NA", "") + "\n")
        metafile.close()

    for bed in os.listdir(dataset[1] + "/"):
        if bed in sample_id:
            shutil.copy(dataset[1] + "/" + bed, combsafe_output + "gdm_metadata/")

    for file in next(os.walk(combsafe_output + 'gdm_metadata/'))[2]:
        if 'narrow' in file:
            shutil.move(combsafe_output + "gdm_metadata/" + file, combsafe_output + "gdm_metadata/narrow")
        else:
            shutil.move(combsafe_output + "gdm_metadata/" + file, combsafe_output + "gdm_metadata/broad")
            
    generate_schema(narrow=False)
    generate_schema(narrow=True)


def extract_data(dataset, factor_list):
    if not os.path.exists(combsafe_output + 'gdm_metadata/'):
        generate_gdm_format(dataset)
    global GMQL_output_dir

    GMQL_output_dir = "./CombSAFE_output/"+ "_".join(factor_list) + "_pipeline/"

    broad = gl.load_from_path(local_path="./CombSAFE_output/gdm_metadata/broad/", parser=gl.parsers.BroadPeakParser())
    narrow = gl.load_from_path(local_path="./CombSAFE_output/gdm_metadata/narrow/", parser=gl.parsers.NarrowPeakParser())
    f_broad = broad[broad['Factor'].isin(factor_list)]
    f_narrow = narrow[narrow['Factor'].isin(factor_list)]
    full_dataset = f_broad.union(f_narrow, left_name="broad", right_name="narrow")
    sem_ann = full_dataset.cover(minAcc=1, maxAcc="Any", groupBy=['Semantic_Annotation', 'Factor'])
    groups_c = sem_ann.group(meta=['Semantic_Annotation'], meta_aggregates = {'n_samp' : gl.COUNTSAMP()})
    final = groups_c[(groups_c['n_samp'] == len(factor_list))]
    results = final.materialize(GMQL_output_dir)

    global destination_path        
    destination_path = GMQL_output_dir + "bed_directory/"

    if not os.path.exists(destination_path):
        os.makedirs(destination_path)
    global gmql_output
    gmql_output = GMQL_output_dir + "files/"
    for file in os.listdir(gmql_output):    
        if file.endswith(".gdm"):
            new_name = file.split(".")[0] + ".bed"
            bedfile = pd.read_csv(gmql_output + file, sep = "\t", names = ["chrom", "start", "stop", "strand", "AccIndex", "JaccardIntersect", "JaccardResult"])
            bedfile = bedfile[bedfile['chrom'].str.len() < 6]
            bedfile.strand = "+"
            bedfile = bedfile[["chrom", "start", "stop", "strand"]]
            bedfile.to_csv(os.path.join(destination_path, new_name), sep="\t", index=False, header=False)
    
    
def load_extracted_data(pool):
    
    global GMQL_output_dir
    GMQL_output_dir = combsafe_output + "_".join(pool) + "_pipeline" + "/" 
    
    global gmql_output
    gmql_output = GMQL_output_dir + "files/"

def add_custom_tracks(tracks_label, path_to_custom_tracks, index = 0):
    
    global custom_path
    custom_path = path_to_custom_tracks
    
    global custom_name
    custom_name = tracks_label
    
    global custom_tracks
    custom_tracks = True
    
    destination_path = GMQL_output_dir + "bed_directory/"
    
    if not os.path.exists(destination_path):
        os.makedirs(destination_path)
    
    custom_file = pd.read_csv(path_to_custom_tracks, sep="\t", index_col=index, header=None).reset_index()
    custom_file.columns = ["chr", "start", "stop"]
    custom_file = custom_file[custom_file['chr'].str.len() < 6]
    custom_file["strand"] = "+"

    custom_file.to_csv(os.path.join(destination_path, custom_path.split("/")[-1].split(".")[0] + ".bed"), sep="\t", index=False, header=False)  
    


def download_custom_tracks(tracks_label, custom_tracks_url):
    
    global custom_name
    custom_name = tracks_label
    
    if not os.path.exists(combsafe_output + "downloaded_files/"):
        os.makedirs(combsafe_output + "downloaded_files/")

    # Download archive

    name = custom_tracks_url.split("/")[-1].split(".")[0]
    custom_file_path = combsafe_output + "downloaded_files/" + custom_tracks_url.split("/")[-1].split(".gz")[0]
    try:
      # Read the file inside the .gz archive located at url
      with urllib.request.urlopen(custom_tracks_url) as response:
            with gzip.GzipFile(fileobj=response) as uncompressed:
                file_content = uncompressed.read()

      # write to file in binary mode 'wb'
      with open(custom_file_path, 'wb') as f:
            f.write(file_content)

    except Exception as e:
        print(e)
    add_custom_tracks(tracks_label, custom_file_path, index = 1)

def segment_files(input_segment_dir=None, n_states=None, custom_segments = False):
    
    segments =[]
    
    if custom_segments:
        input_segments = input_segment_dir
    else:
        input_segments = chromHMM_output
    
    for file in os.listdir(input_segments):
        if not custom_segments:        
            if "segment" in file:
                segments.append(os.path.join(input_segments, file))
        else:
            segments.append(os.path.join(input_segments, file))        

    chrs = ['chr'+ str(x) for x in list(range(1,23))+['X', 'Y']]

    global n_states_model
    n_states_model = GMQL_output_dir + str(n_states) + "_state_model/"

    newDict = {} 
    for beds in segments:
        file = open(beds).readlines()
        for i2 in file:
            window2 = i2.rstrip().split("\t")
            if window2[0] in newDict:
                templist2 = newDict.get(window2[0])
                templist2.append(int(window2[1]))
                templist2.append(int(window2[2]))
            else:
                newDict[window2[0]] = [int(window2[1]), int(window2[2])]

    targetDict={}
    for key, value in newDict.items():
        targetDict[key] = sorted(list(set(value)))

    if not os.path.exists(n_states_model + 'segmentated_files/'):
        os.makedirs(n_states_model + 'segmentated_files/')

    for files in segments:
        file_n = open(files).readlines()
        myDict = {} 
        for i in file_n:
            window = i.rstrip().split("\t")
            if window[0] in myDict:
                templist = myDict.get(window[0])
                templist.append((int(window[1]), int(window[2]), window[3]))   
            else:
                myDict[window[0]] = [(int(window[1]), int(window[2]), window[3])]

        segmentated_file=open(n_states_model + 'segmentated_files/' + files.split("/")[-1], "w")
        for z in range(0, len(chrs)):
            n=0
            v = targetDict.get(chrs[z])
            v2= (list(set(v)))
            v2.sort()
            for reg in myDict.get(chrs[z]):
                while(reg[1] > v2[n]):
                    segmentated_file.write(chrs[z] + "\t" + str(v2[n]) + "\t" + str(v2[n+1]) + "\t" + reg[2] + "\n")
                    n+=1
        segmentated_file.close()     
    
def load_funct_states_metric():
    
    df = pd.read_csv(emissions_path, sep="\t", index_col=0)
    data = df.values.tolist()
    state_matrix=[]
    for i in range(df.shape[0]):
        single_state=[]
        for j in range(df.shape[0]):
            dst = 1 - stats.pearsonr(data[i], data[j])[0]
            if data[i] != data[j]:
                dst +=1
            single_state.append(dst)
        state_matrix.append(single_state)

    dst_matrix = pd.DataFrame.from_records(state_matrix)
    global nrm_df
    nrm_df=(dst_matrix-dst_matrix.min())/(dst_matrix.max()-dst_matrix.min())
    nrm_df.columns = list(df.index)
    nrm_df.index = list(df.index)
    nrm_df = ssd.squareform(nrm_df.to_numpy()[np.triu_indices_from(nrm_df, k=1)])
    return nrm_df


def functional_states_distance(seq1, seq2):
    distance = 0
    for x, y in zip(seq1, seq2):
        distance+=(nrm_df[int(x)-1][int(y)-1])
    return(distance)
    
def identify_functional_states(chromHMM_path, number_of_states, n_core):
    
    global n_states_model
    n_states_model = GMQL_output_dir + str(number_of_states) + "_state_model/"
    
    global n_states
    n_states = number_of_states
    
    txt_file = GMQL_output_dir + "cellmarkfiletable.txt"
    binarized_files = GMQL_output_dir + "binarized_files"
    cellmarkfiletable = open(txt_file, "w")
    se=[]
    
    gmql_output = GMQL_output_dir + "files/"
    for file in os.listdir(gmql_output):
        if file.endswith(".meta"):
            metafile = open(gmql_output + file).readlines()
            v=[]
            name = file.split(".")[0]
            for line in metafile:
                key = line.rstrip().split("\t")
                if key[0] == "Factor" or key[0] == "Semantic_Annotation":
                    v.append(key[1])
            se.append(v[1])
            cellmarkfiletable.write(v[1] + "\t" + v[0] + "\t" +  name + ".bed" + "\n")            
            
    if custom_tracks:
        for sem_ann in list(set(se)):
            cellmarkfiletable.write(sem_ann + "\t" + custom_name + "\t" +  custom_path.split("/")[-1].split(".")[0] + ".bed" + "\n")        
    cellmarkfiletable.close()
    
    binarizing = " ".join(['java', '-mx32600M', '-jar', chromHMM_path + 'ChromHMM.jar', 'BinarizeBed', '-center', './ChromHMM/CHROMSIZES/hg38.txt', destination_path, txt_file, binarized_files])
    subprocess.call(binarizing, shell=True)
    
    learning = " ".join(['java', '-mx32600M', '-jar', chromHMM_path + 'ChromHMM.jar', 'LearnModel', '-p ' + str(n_core),  '-r 10', binarized_files, n_states_model, str(number_of_states), 'hg38'])
    subprocess.call(learning, shell=True)
    
    global chromHMM_output
    chromHMM_output = n_states_model + 'ChromHMM_output/'
    
    if not os.path.exists(chromHMM_output):
        os.makedirs(chromHMM_output)    
    
    chrom_out = os.listdir(n_states_model)
    for f in chrom_out:
        if ('_' + str(number_of_states) in f):
            shutil.move(n_states_model + f, chromHMM_output)
    
    global graphs_path
    graphs_path = n_states_model + "graphs/"
    
    global emissions_path
    emissions_path = chromHMM_output + "emissions_" + str(n_states) + ".txt"
    
    if not os.path.exists(graphs_path):
        os.makedirs(graphs_path)
    segment_files(n_states=number_of_states)
    load_funct_states_metric()
    full_df = load_functional_states_dataframe()
    return full_df


def load_custom_segments(input_segment_dir, num_states):
        
    global n_states
    n_states = num_states
        
    global n_states_model
    n_states_model = GMQL_output_dir + str(num_states) + "_state_model/"  
    
    global custom_segments   
    custom_segments = True
    
    global emissions_path
    emissions_path = input_segment_dir + next(os.walk(input_segment_dir))[2][0]
    
    global graphs_path
    graphs_path = n_states_model + "graphs/"
        
    if not os.path.exists(graphs_path):
        os.makedirs(graphs_path)
    
    bed_segments_path = input_segment_dir + next(x for x in next(os.walk(input_segment_dir))[1] if "." not in x)
    
    segment_files(bed_segments_path, num_states, custom_segments=True)
    load_funct_states_metric()
    full_df = load_functional_states_dataframe()
    return full_df
    
    
    
def show_emission_graph(custom_palette = False):
    
    global color_list
    if custom_palette:
        color_list = custom_palette
        
    n_factor = GMQL_output_dir.split("/")[2].count("_")                 
    emis_txt = pd.read_csv(emissions_path, sep="\t", index_col=0)

    g_emi = sns.heatmap(emis_txt, cmap="Blues", cbar=False)
    plt.title('Emission Parameters', fontsize=15, pad = 20)
    plt.xticks(rotation = 90, fontsize=15)
    plt.yticks(rotation = 0, fontsize=15)
    plt.ylabel("State (Emission order)", fontsize = 15)
    g_emi.tick_params(left=False, bottom=False)
    fig = plt.gcf()

    if custom_tracks:
        fig.set_size_inches((n_factor+1)/3, (n_states+1)/3)
    else:
        fig.set_size_inches(n_factor/3, n_states/3)

    for l in g_emi.yaxis.get_ticklabels():
        a = l.get_text()
        l.set_backgroundcolor(color_list[int(a)-1])

    plt.savefig(graphs_path + 'chromatin_states_emissions_' + str(n_states) + '.jpeg', bbox_inches='tight', format="jpeg")
    plt.savefig(graphs_path + 'chromatin_states_emissions_' + str(n_states) + '.svg', bbox_inches='tight', format="svg")
    plt.show()

    emissions = Image.open(graphs_path + 'chromatin_states_emissions_' + str(n_states) + '.jpeg')
    emissions_path_jpeg = graphs_path + 'chromatin_states_emissions_' + str(n_states) + '.jpeg'
    em = Image.open(emissions_path_jpeg)
    em = em.resize((200, 350), Image.ANTIALIAS)


def load_functional_states_dataframe():
    #create dataframe with segmentated files
    labels = []
    files =[]
    for file in os.listdir(n_states_model + "segmentated_files/"):
        if "segment" in file:
            files.append(file)
            n_states = n_states_model.split("/")[3].split("_")[0]
            labels.append(file.split("_" + str(n_states) + "_")[0])

    functional_states_dataframe = pd.read_csv(n_states_model + 'segmentated_files/' + files[0], sep="\t", header=None)
    for i in range(1,len(files)): 
        states = pd.read_csv(n_states_model + 'segmentated_files/' + files[i], sep="\t",  usecols=[3], header=None)
        functional_states_dataframe[str(i+3)] = states

    functional_states_dataframe.columns = ('chr start stop'.split(" ") + labels)
    functional_states_dataframe = functional_states_dataframe.set_index(["chr", "start", "stop"]).applymap(lambda x: int(x[1:]))
    return(functional_states_dataframe)

def get_concat_h(im1, im2):
    dst = Image.new('RGB', (im1.width + im2.width, im2.height), color=(255,255,255,0))
    dst.paste(im2, (0, 0))
    dst.paste(im1, (im2.width, 0))
    return dst    
        
        
def single_gene_analysis(functional_states_dataframe, path_to_gene_list_file, distance_metric):
    #Run "pyensembl install --release 97 --species homo_sapiens"
    
    gene_list = open(path_to_gene_list_file, "r").readlines()
    nrm_df = load_funct_states_metric()
    plt.rcParams.update({'figure.max_open_warning': 0})

    list_name = path_to_gene_list_file.split("/")[-1].split(".")[0]
    single_gene_path = graphs_path + "single_gene_analysis/" + list_name + "/"

    if not os.path.exists(single_gene_path):
        os.makedirs(single_gene_path)

    for i in gene_list:
        gene_name = (i.rstrip())

        try:
            gene_data = data.genes_by_name(gene_name)[0]
        except:
            pass

        region = functional_states_dataframe.loc[["chr" + str(gene_data.contig)]].query(f'start >= {str(gene_data.start)} & stop <= {str(gene_data.end)}')
        number_of_states = int(n_states_model.split("/")[3].split("_")[0])    

        if type(distance_metric) is str:
            
            link=None
            column_cluster = distance_metric            

        else:
            distance_matrix_cols = pairwise_distances(region.T, metric=distance_metric)
            link = linkage(ssd.squareform(distance_matrix_cols), method="ward")
            column_cluster=True


        gene_coordinates = "chr" + str(gene_data.contig) + ":" + str(f"{gene_data.start:,}") + "-" +  str(f"{gene_data.end:,}")
        gene_heatmap = sns.clustermap(region, col_linkage=link, col_cluster = column_cluster, row_cluster=False, cmap=color_list[0:(number_of_states)], vmin=1, vmax=n_states, xticklabels=True, yticklabels=3, figsize=(10, 35), linewidths=.1)


        plt.suptitle("Gene: \"" + gene_name + "\"; " + "Coordinates: \"" + gene_coordinates + "\"; " + "Strand: \"" + gene_data.strand + "\"; ", fontsize=18, x=0.55, y=1)
        gene_heatmap.cax.set_visible(False)
        gene_heatmap.ax_heatmap.set_xticklabels(gene_heatmap.ax_heatmap.get_xmajorticklabels(), fontsize = 12)
        
        plt.savefig(single_gene_path + gene_name + "_data_driven_heatmap_temp.jpeg", bbox_inches='tight', format="jpeg", dpi=100)

        hitmap_path = single_gene_path + gene_name + "_data_driven_heatmap_temp.jpeg"
        emission_path = graphs_path + 'chromatin_states_emissions_' + str(number_of_states) + '.jpeg'
        em_big = Image.open(emission_path)
        em = em_big.resize((int(em_big.size[0]/1.0), int(em_big.size[1]/1.0)), Image.ANTIALIAS)
        hit_big = Image.open(hitmap_path)
        hit = hit_big.resize((int(hit_big.size[0]), int(hit_big.size[1])), Image.ANTIALIAS)
        get_concat_h(em, hit).save(single_gene_path + gene_name + "_data_driven_heatmap.jpeg")
        os.remove(single_gene_path + gene_name + "_data_driven_heatmap_temp.jpeg")
        plt.cla()
        
def load_reducted_df():
    cols = ['chr','start','stop']
    reducted_df = pd.read_csv(n_states_model + "segmentated_files/reducted_genome.txt", sep="\t")
    reducted_df.columns = cols + reducted_df.columns.tolist()[len(cols):]
    reducted_df.set_index(["chr", "start", "stop"], inplace=True)
    return(reducted_df)

def show_emission_graph_with_coverage(new_df, custom_palette = False):
    
    global color_list
    if custom_palette:
        color_list = custom_palette

    n_factor = GMQL_output_dir.split("/")[2].count("_")
    n_states = int(n_states_model.split("/")[3].split("_")[0])
    emis_txt = pd.read_csv(emissions_path, sep="\t", index_col=0)

    occurrencies = new_df.apply(pd.Series.value_counts).fillna(0)
    occurrencies['chromatin state percentage'] = occurrencies.sum(axis=1)
    x = np.array(list(occurrencies.index))
    y = np.array(occurrencies['chromatin state percentage'])
    porcent = 100.*y/y.sum()
    emis_txt['Coverage %'] = porcent.round(decimals=2)

    #fig = plt.figure(figsize=(3,7))
    ax1 = plt.subplot2grid((1,10), (0,0), colspan=8)
    ax2 = plt.subplot2grid((1,10), (0,8), colspan=2)

    g = sns.heatmap(emis_txt[list(emis_txt.columns[0:-1])], ax=ax1, cmap="Blues", cbar=False, linewidths=.5)
    g.set_title('Emission Parameters', fontsize=20, pad = 20)
    g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 15)
    g.set_xticklabels(g.get_xticklabels(), rotation = 90, fontsize = 15)
    g.set_ylabel("State (Emission order)",fontsize=15, labelpad=20)
    #g.set_xlabel("Marks",fontsize=15, labelpad=20)
    g.tick_params(left=False, bottom=False)

    fig = plt.gcf()
    if custom_tracks:
        fig.set_size_inches((n_factor+(n_factor+1))/3, (n_states+(n_states))/4)
    else:
        fig.set_size_inches(n_factor/2, n_states/2.5)
    for l in g.yaxis.get_ticklabels():
        a = l.get_text()
        l.set_backgroundcolor(color_list[int(a)-1])

    g2 = sns.heatmap(emis_txt[list(emis_txt.columns)[-1]].to_frame(), ax=ax2,  annot=True, cmap="Blues", cbar=False, annot_kws={"size": 10}, linewidths=.5)
    g2.tick_params(left=False, bottom=False)
    g2.set_xticklabels(g2.get_xticklabels(), rotation = 90, fontsize = 15)
    g2.set_ylabel('')
    g2.set_yticklabels('')

    plt.savefig(graphs_path + 'chromatin_states_emissions_coverage_' + str(n_states) + '.jpeg', bbox_inches='tight', format="jpeg" , dpi=100)
    plt.savefig(graphs_path + 'chromatin_states_emissions_coverage_' + str(n_states) + '.svg', bbox_inches='tight', format="svg" , dpi=100)
    plt.show()
    
    
def genome_reduction(full_df, treshold):
    
    list_of_dataframe=[]
    intra_chrom_dict = {}
    chrs = ['chr'+ str(x) for x in list(range(1,23))+['X', 'Y']]
    inactive_chromatin = np.argmax(np.bincount(full_df.values.flat))
    nrm_df = load_funct_states_metric()
    distance_metric = "hamming"
    for chrom in chrs:
        temp_df = full_df.loc[chrom].drop_duplicates()
        print("\r", chrom, end=" ")
        temp_df2 = temp_df[(temp_df==inactive_chromatin).sum(axis=1)/len(temp_df.columns) <= float(treshold/100)]
        data_ = temp_df2.values.tolist()
        
        if type(distance_metric) is str:
            min_samples_=2
        else:
            min_samples_=20

        clusterer_f = hdbscan.HDBSCAN(metric=distance_metric, min_cluster_size=5, min_samples=min_samples_)
        clusterer_f.fit(data_)

        max_clusters = (clusterer_f.labels_).max()
        indices = (clusterer_f.labels_)
        indices = np.where(indices==-1, max_clusters+1, indices)

        temp_df2 = temp_df2.copy()
        temp_df2.reset_index(inplace=True)
        temp_df2['chr'] = chrom
        temp_df2 = temp_df2.set_index(["chr", "start", "stop"])
        for clusters in list(set(indices)):
            res_list = list(filter(lambda x: indices[x] == clusters, range(len(indices))))
            c_index = np.random.choice(res_list, 1)
            c_index = list(set(c_index))
            c_index.sort()     
            for positions in c_index:
                list_of_dataframe.append(temp_df2.iloc[int(positions)])
            cluster_of_regions = temp_df2.iloc[res_list].index.tolist()
            list_of_genes = []
            for region in cluster_of_regions:
                chr_n = region[0][3:]
                if(chr_n.isnumeric()):
                    chr_n = int(chr_n)
                gene_names = data.gene_names_at_locus(contig=chr_n , position=region[1])
                if len(gene_names) > 0:
                    list_of_genes.append(gene_names[0])
            key = "_".join(map(str, temp_df2.iloc[int(c_index[0])].name))
            intra_chrom_dict[key] = (list(set(list_of_genes)))
    reducted_df = pd.concat(list_of_dataframe, axis=1).T

    new = []
    for i in reducted_df.columns:
        new.append(i.split("_" + n_states_model.split("/")[3].split("_")[0] + "_")[0])
    reducted_df.columns = new
    reducted_df.to_csv(n_states_model + 'segmentated_files/' + "reducted_genome.txt", sep='\t', encoding='utf-8')
    save_dict_to_file(intra_chrom_dict)
    return(reducted_df)

def data_driven_heatmap(reducted_df, distance_metric, min_clust_size, min_sampl):
    copy_df = reducted_df.copy()
    nrm_df = load_funct_states_metric()
    distance_matrix_cols = pairwise_distances(copy_df.T, metric=distance_metric)
    link = linkage(ssd.squareform(distance_matrix_cols), method="ward")

    distance_matrix_rows = pairwise_distances(copy_df, metric=distance_metric)
    if type(distance_metric) is str:
        clusterer_f = hdbscan.HDBSCAN(metric=distance_metric)
    else:
        clusterer_f = hdbscan.HDBSCAN(metric="hamming", min_cluster_size=min_clust_size, min_samples=min_sampl)

    clusterer_f.fit(distance_matrix_rows)
    indices_f = clusterer_f.labels_
    set_of_indices = set(indices_f)
    set_of_indices -= {0,-1}

    copy_df['cluster'] = indices_f
    copy_df.drop(copy_df[copy_df.cluster < 1].index, inplace=True)
    copy_df = copy_df.sort_values('cluster')
    copy_df2 = copy_df.copy()
    cluster = copy_df.pop("cluster")

    if len(cluster.unique()) <= 20:
        side_color = sns.color_palette("tab20")
    else:
        side_color = sns.crayons.values()

    lut = dict(zip(cluster.unique(), side_color))
    row_colors = cluster.map(lut)

    n_states = int(n_states_model.split("/")[3].split("_")[0])   


    try:
        g = sns.clustermap(copy_df, row_colors = row_colors, col_linkage=link, row_cluster=False, vmin=1, vmax=n_states,cmap=color_list[0:n_states], xticklabels=True, yticklabels=False, figsize=(10, 25), linewidths=0, cbar_pos = (0.1, 0.65, 0, 0))
    except ValueError:  #raised if `y` is empty.
        pass
    #g.cax.set_visible(False)
    g.ax_row_dendrogram.set_visible(False)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = 12)
    g.ax_heatmap.set_ylabel("")

    from matplotlib.patches import Patch
    legend_handles =[]
    cmap = list(side_color)[0:len(set_of_indices)]
    for i in range(len(cmap)):
        legend_handles.append(Patch(color=cmap[i], label='Cluster_' + str(i+1)))

    plt.legend(handles=legend_handles, title = "cluster type", ncol=2, bbox_to_anchor=[-0.1, 0.01], loc='right', fontsize=12, handlelength=.9)

    fig = plt.gcf()

    cax = fig.axes[-1]
    cax.tick_params(labelsize=1, size=0)

    fig.savefig(graphs_path + "data_driven_heatmap_temp.jpeg", bbox_inches='tight', format="jpeg", dpi=100)

    hitmap_path = graphs_path + "data_driven_heatmap_temp.jpeg"
    emissions_c_path = graphs_path + 'chromatin_states_emissions_coverage_' + str(n_states) + '.jpeg'
    em_big = Image.open(emissions_c_path)
    em = em_big.resize((int(em_big.size[0]/1.5), int(em_big.size[1]/1.5)), Image.ANTIALIAS)
    hit_big = Image.open(hitmap_path)
    hit = hit_big.resize((int(hit_big.size[0]), int(hit_big.size[1])), Image.ANTIALIAS)

    get_concat_h(em, hit).save(graphs_path + "data_driven_heatmap.jpeg")
    os.remove(graphs_path + "data_driven_heatmap_temp.jpeg")
    return(copy_df2)

def gene_ontology_enrichment_analysis(copy_df2, distance_metric, goea_tool):
    
    if goea_tool == "great":

        robjects.r('install.packages("BiocManager", repos="http://cran.r-project.org")')
        robjects.r('BiocManager::install("rGREAT", update=FALSE)')
        great = importr('rGREAT')    

    mg = mygene.MyGeneInfo()

    global clusters_path
    clusters_path = graphs_path + "GOEA/" #<-
    nrm_df = load_funct_states_metric()
    if not os.path.exists(clusters_path):
        os.makedirs(clusters_path)

    ens_list = [] 


    if len(copy_df2["cluster"].unique()) <= 20:
        side_color = sns.color_palette("tab20")
    else:
        side_color = sns.crayons.values()

    for cluster_number in list(set(copy_df2['cluster'])):
        clustered_df = copy_df2[copy_df2["cluster"] == cluster_number]
        clustered_df.pop("cluster")
        cluster_name = "Cluster_" + f"{cluster_number:02}"
        cluster_color = list(side_color)[cluster_number-1]
        cluster_of_regions = clustered_df.index.tolist()

        heatmap = sns.clustermap(clustered_df, vmin=1, vmax=n_states, cmap=color_list[0:n_states], col_cluster=distance_metric, row_cluster=distance_metric, xticklabels=True, yticklabels=False, figsize=(10, 30), linewidths=.1, annot=True)
        heatmap.cax.set_visible(False)
        #heatmap.ax_row_dendrogram.set_visible(False)
        heatmap.ax_heatmap.set_xticklabels(heatmap.ax_heatmap.get_xmajorticklabels(), fontsize = 14)
        fig = plt.gcf()
        fig.text(0.52, 0.99, cluster_name, ha="left", va="bottom", size="x-large",color=cluster_color, weight="bold")

        heatmap_path = clusters_path + cluster_name + "_heatmap_temp.jpeg"    
        fig.savefig(heatmap_path, bbox_inches='tight', format="jpeg", dpi=100)         
        heat_big = Image.open(heatmap_path)
        heat = heat_big.resize((int(heat_big.size[0]), int(heat_big.size[1])), Image.ANTIALIAS)
        emission_path = graphs_path + 'chromatin_states_emissions_' + str(n_states) + '.jpeg'
        em_big = Image.open(emission_path)
        em = em_big.resize((int(em_big.size[0]/1), int(em_big.size[1]/1)), Image.ANTIALIAS)
        get_concat_h(em, heat).save(clusters_path + cluster_name + "_heatmap.jpeg")
        os.remove(clusters_path + cluster_name + "_heatmap_temp.jpeg")

        intra_chrom_dict = load_dict_from_file()

        list_of_genes = []
        for reg in cluster_of_regions:
            key = "_".join(map(str, reg))
            if str(key) in intra_chrom_dict:
                list_of_genes += intra_chrom_dict[str(key)]
        
        if list_of_genes:
            
            if goea_tool == "great":

                #robjects.r('install.packages("BiocManager", repos="http://cran.r-project.org")')
                #robjects.r('BiocManager::install("rGREAT", update=FALSE)')
                #great = importr('rGREAT')

                bed_file = open(clusters_path + cluster_name + ".bed", "w")
                for gene in list_of_genes:
                    gene_info = data.genes_by_name(gene)
                    bed_file.write("chr" + gene_info[0].contig + "\t" + str(gene_info[0].start) + "\t" + str(gene_info[0].end) + "\n")
                bed_file.close()

                bed_sample = pd.read_csv(clusters_path + cluster_name + ".bed", sep="\t", names=["chr", "start", "stop"])

                job = great.submitGreatJob(bed_sample, species = organism)
                tb = great.getEnrichmentTables(job)
                tb[1].head(20).to_csv(clusters_path + cluster_name + "_bp_enrichment_analysis.txt", sep='\t')
                tb[0].head(20).to_csv(clusters_path + cluster_name + "_mf_enrichment_analysis.txt", sep='\t')

            if goea_tool == "goatools":

                obo_fname = download_go_basic_obo()
                fin_gene2go = download_ncbi_associations()
                obodag = GODag("go-basic.obo")

                objanno = Gene2GoReader(fin_gene2go, taxids=[9606]) ################################################
                ns2assoc = objanno.get_ns2assc()
                for nspc, id2gos in ns2assoc.items():
                    print("{NS} {N:,} annotated human genes".format(NS=nspc, N=len(id2gos)))

                background=[]
                Entrez.email = 'A.N.Other@example.com'
                handle = Entrez.esearch(db='gene', term='"Homo sapiens"[Organism] AND ("genetype protein coding"[Properties] AND gene_nucleotide_pos[filter])')
                record = Entrez.read(handle)
                count = record['Count']
                handle = Entrez.esearch(db='gene', term='"Homo sapiens"[Organism] AND ("genetype protein coding"[Properties] AND gene_nucleotide_pos[filter])', retmax = count)
                record = Entrez.read(handle)
                for id_genes in record['IdList']:
                    background.append(int(id_genes))

                goeaobj = GOEnrichmentStudyNS(
                        background, # List of human protein-coding genes
                        ns2assoc, # geneid/GO associations
                        obodag, # Ontologies
                        propagate_counts = False,
                        alpha = 0.05, # default significance cut-off
                        methods = ['fdr_bh']) # defult multipletest correction method

                cluster_file = open(clusters_path + cluster_name + "_genes.txt", "w")
                for i in list_of_genes:
                    cluster_file.write(i +  "\n")
                cluster_file.close()            

                gene_query = mg.querymany(list_of_genes, scopes='symbol', species=9606, entrezonly=True)
                gene_id_ = []
                for i in gene_query:
                    if '_id' in i:
                        if 'ENS' not in i['_id']:
                            gene_id_.append(int(i['_id']))
                        else:
                            ens_list.append(i['_id'])
                    else:
                        pass

                geneids_study = gene_id_
                goea_results_all = goeaobj.run_study(geneids_study, prt=None)
                goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]

                #goeaobj.wr_txt(clusters_path + cluster_name + "_gene_ontology_enrichment_analysis.txt", goea_results_sig)
                goeaobj.wr_xlsx(clusters_path + cluster_name + "_gene_ontology_enrichment_analysis.xlsx", goea_results_sig)
                #plot_results(clusters_path + cluster_name + "_gene_ontology_enrichment_analysis.png", goea_results_sig)
                
def pca_analysis(df, n):
    pca = PCA(n_components=n)
    pca.fit(df)
    map_ = pd.DataFrame(pca.components_, columns=df.columns)
    g = sns.clustermap(map_, cmap="Blues", figsize=(10,5), row_cluster=False)
    g.fig.subplots_adjust(right=0.7)
    g.ax_cbar.set_position((0.1, .2, .03, .4))
    loadings = pca.components_
    return(loadings)
    
def show_distance_matrix():
    distance = load_funct_states_metric()
    cols = list(range(len(load_funct_states_metric()) + 1))[1:]
    dist = pd.DataFrame(nrm_df, columns=cols, index =cols)
    sns.clustermap(dist, cmap="Blues_r", figsize=(8, 8))
    
def show_cluster(cluster_number, analysis_type):
    path = clusters_path + "/Cluster_" + f"{cluster_number:02}_" + analysis_type + "_enrichment_analysis.txt"
    clust = pd.read_csv(path, sep = "\t", index_col=0)
    return(clust)    
