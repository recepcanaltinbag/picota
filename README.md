

<p align="center">
<img src="logo/picota_logo.png" alt="Picota logo" width="150"/>
</p>

An extendible [P]()ipeline for [I]()dentification of [CO]()mposite [T]()ransposons from [A]()ssembly graphs



- [PICOTA](#picota)
  - [Installation](#installation)
  - [Downloading DBs](#downloading-dbs)
  - [Test](#test)
  - [Quick Start](#quick-start)
  - [Usage](#usage)
    - ['sra_download' module: Downloading Raw Read Data](#sra_download-module-downloading-raw-read-data)
    - ['assembly' module: Generating Assembly Graph](#assembly-module-generating-assembly-graph)
    - ['analysis' module: Analysis of Assembly Graph](#analysis-module-analysis-of-assembly-graph)
    - ['db' module: Downloading Databases](#db-module-downloading-databases)
    - ['scoring' module: Scoring Composite Transposon Candidates](#scoring-module-scoring-candidates)
    - ['all' module: All in One](#all-module-all-in-one-command)
  - [Method](#method)
    - [Composite Transposons](#composite-transposons)
    - [Pipeline of Picota](#pipeline-of-picota)
    - [Assembly Graph Analysis](#assembly-graph-analysis)
    - [Scoring Functions](#scoring-functions)
    - [Future Analysis](#future-analysis)
    - [Clustering](#clustering)
  - [FAQ](#faq)
    - [How Can I Cite Picota](#how-can-i-cite)
  - [License](#license)



# PICOTA
An extendible [P]()ipeline for [I]()dentification of [CO]()mposite [T]()ransposons from [A]()ssembly graphs

## Installation


### Requirements

- Linux (not tested on Windows and MacOS)
- [Phyton](https://www.python.org/) 3.8 or later

### Required Packages 

- [BioPython](https://biopython.org/) (biopython-1.4.0)
- [pandas](https://pandas.pydata.org/) (pandas-1.4.3)
- [requests](https://pypi.org/project/requests/) (requests-2.27.1)
- [tqdm](https://github.com/tqdm/tqdm) (tqdm-4.62.3)

### Required Tools
- [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit)
- [fastp](https://github.com/OpenGene/fastp)
- [spaDES](https://github.com/ablab/spades)
- [prodigal](https://github.com/hyattpd/Prodigal)
- [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)



### Installation of PICOTA

Although other environment management systems can also be used, we recommend to use conda:

```
conda update conda --all
conda create --name picota python=3.8
conda activate picota

```

You can use pip to install dependecies

```
pip install biopython
pip install pandas
pip install requests
pip install tqdm
```
```
wget https://github.com/recepcanaltinbag/picota/archive/refs/heads/main.zip
unzip picota-main.zip
cd picota-main
python picota.py -h 
```

## Downloading DBs

## Test

## Quick Start

## Usage

PICOTA is consist of subcommands which enable you to divide or combine them as your needs. 

### 'sra_download' module: Downloading Raw Read Data

PICOTA needs [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) for this module. Install SRA Toolkit and make it available with adding it in your PATH or you can give the path of fastq_dump with optional parameter of `--path_of_fastq_dump`. `fastq_dump` is used to obtain fastq files from sra files, so without it, the pipeline will not start. 

|__type__ |__command__ |__description__ |
| --- | --- | --- |
| required | `sra_list` | `.txt` file includes list of SRA Accessions splitted with new lines |
| required | `out_folder` | out folder for the pipeline, sra files will be downloaded in `sra_download` folder in out folder |
| optional | `--path_of_fastq_dump` | `default='fastq-dump'`, path of `fastq-dump` program, if you dont have sra-toolkit in PATH, you need to give the path of the `fastq-dump` |
| optional | `--keep_sra_file` | `default=False` to keep `sra file`, because they will be deleted after usage in default |
| optional | `--force` | `default=False` to force overwrite to output |


### 'assembly' module: Generating Assembly Graph

PICOTA needs [fastp](https://github.com/OpenGene/fastp) to filter raw reads and [SPAdes](https://github.com/ablab/spades) to generate assembly graph (`.gfa file`)


|__type__ |__command__ |__description__ |
| --- | --- | --- |
| required | `rawread_folder` | folder path of raw reads |
| required | `out_folder` | out folder for the pipeline, sra files will be downloaded in `gfa_files` and `contig_files` folder in out folder |
| optional | `--batch` | `default=False`, if given, batch analysis will be made (PICOTA assumes there are folders which have fastq file inside in the given rawread_folder.), else in default, raw read folder contains fastq files (not folders). Using this option is good when you want to analysis more than one genome |
| optional | `--threads` | `default=2`, number of threads for assembly process |
| optional | `--k_mer_list` | `default=79,99`, Longer list make analysis longer, you can add different k-mers to the list with comma. They must be odd and maximum 127. |
| optional | `--quiet` | `default=False`, activate this if you want quiet mode for assembly |
| optional | `--reads_to_process` | `default=5000000`, fastp option, you can increase if you need more reads to process |
| optional | `--fastp_q` | `default=20`, fastp option, quality threshold for the raw reads |
| optional | `--keep_temp_files` | `default=False`, to keep temp assembly files, otherwise just gfa and fasta file will be stayed |
| optional | `--skip_filtering` | `default=False`, to skip filtering process |
| optional | `--meta` | `default=''`, if you use metagenomic data |
| optional | `--path_of_spades` | `default='spades.py'`, spades.py path if you dont have spades.py in PATH |
| optional | `--path_of_fastp` | `default='fastp'`, fastp path if you dont have fastp in PATH |


### 'analysis' module: Analysis of Assembly Graph

|__type__ |__command__ |__description__ |
| --- | --- | --- |
| required | `gfa_folder` | folder path of `.gfa` files |
| required | `out_folder` | out folder for the pipeline, sra files will be downloaded in `analysis` folder in out folder |
| optional | `--find_all_path` | `default=True`, finding all possible paths is best, but sometimes it can take too long running time, I recommend to try default unless you need very high speed |
| optional | `--path_limit` | `default=15`, High number of path limit will slow the analysis, I recommend a value between 10-25, because if the path is longer than 25, it is probably not a composite transposon unless the average node lenghts are very small |
| optional | `--min_size_of_cycle` | `default=1500`, min size for composite transposons (assumption in cyclic form) |
| optional | `--max_size_of_cycle` | `default=100000`, max size for composite transposons (assumption in cyclic form) |
| optional | `--name_prefix_cycle` | `default=''`, prefix of cycles in output fasta file |
| optional | `--min_component_number` | `default=1`, min size of components |
| optional | `--max_component_number` | `default=25`, min size of components |
| optional | `--k_mer_sim` | `default=200`, splitting k_mers when looking similarities, kind of ANI |
| optional | `--threshold_sim` | `default=80`, threshold in similarity elimination |


### 'db' module: Downloading Databases

### 'scoring' module: Scoring Candidates

### 'all' module: All In One Command

## Method

### Composite Transposons

Transposons are the DNA sequences can mobilize and alter its position in genome, changing the genome size. Also, they can contribute to genome and gene evolution. Composite transposons are composed of two transposon and genetic material inside the flanking transposons. Since transposition event is somewhat stochastic, these transposons that are close to each other can act as one transposon and faciliate the transfer of the genetic material between.

We are in the brink of a "post-antibiotic" era. As more an more bacteria gain resistance to antibiotics due to misuse and overuse of antibiotics, it is more important now then ever to understand how bacteria transmit antibiotic resistance genes to each other if we don't want to go back to an era where even small scratches are deadly. Composite transposons also often carry antibiotic resistence genes, and sometimes metabolic genes such as degradation of xenobitiotics and metal resistance. Due to the nature of the transposons, these genes can be mobilized and potentially introduced to other bacteria. Therefore, it is crucial to identify and understand composite transposons if we want to understand the transfer of antibiotic resistance between different organisms and find interesting metabolic functions of genes that can be transferred.

There are many tools to find the transposons but there is a few tools related to composite transposons. One of them is [TnComp_finder](https://github.com/danillo-alvarenga/tncomp_finder). The limitation of TnComp_Finder is that it can be only used in complete genomes. For the incomplete genomes, it will probably fail because short-read sequencing (most of the time) is not able to capture composite transposons due to repetetive nature of these transposons. Because the genome databases are mostly consists of incomplete genomes, a tool that can work also in incomplete genomes is essential, so we developed PICOTA.

### Pipeline of PICOTA

To capture the composite transposons in incomplete genomes, assembly graphs built from the raw reads can be very helpful unlike using fasta sequences since the links that are between different parts of genomics sequences are lost in fasta files. Graph algorithms can be used to find transposon motifs, previously known transposons and possible novel composite transposons can be identified even *de-novo*. Then, the candidates can be investigated in existing gene databases for functional analysis. 

### Assembly Graph Analysis



### Scoring Functions
### Future Analysis
### Clustering

## FAQ

### How can I cite


## Licence
