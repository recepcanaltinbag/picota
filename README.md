

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

[see](#db-module-downloading-databases)

## Test

## Quick Start

## Usage

PICOTA is consist of subcommands which enable you to divide or combine them as your needs. 

### 'sra_download' module: Downloading Raw Read Data

PICOTA needs [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) for this module. Install SRA Toolkit and make it available with adding it in your PATH or you can give the path of fastq_dump with optional parameter of `--path_of_fastq_dump`. `fastq_dump` is used to obtain fastq files from sra files, so without it, the pipeline will not start. 

`sra_list.txt` file (list of SRA Accessions splitted with new lines) should look like this:

```
ERR9786618
SRR6242084
SRR10869808
```
Example usage of the sra_download subcommand with default parameters
```
python picota.py sra_download sra_list.txt out_folder_of_pipeline 
```

Detailed usage:


|__type__ |__command__ |__description__ |
| --- | --- | --- |
| required | `sra_list` | `.txt` file includes list of SRA Accessions splitted with new lines |
| required | `out_folder` | out folder for the pipeline, sra files will be downloaded in `sra_download` folder in out folder |
| optional | `--path_of_fastq_dump` | `default='fastq-dump'`, path of `fastq-dump` program, if you dont have sra-toolkit in PATH, you need to give the path of the `fastq-dump` |
| optional | `--keep_sra_file` | `default=False` to keep `sra file`, because they will be deleted after usage in default |
| optional | `--force` | `default=False` to force overwrite to output |


### 'assembly' module: Generating Assembly Graph

PICOTA needs [fastp](https://github.com/OpenGene/fastp) to filter raw reads and [SPAdes](https://github.com/ablab/spades) to generate assembly graph (`.gfa file`)

There is two mode in this module, so be careful when you enter the input folder to the command

If, in `default` mode, Picota assumes there is just one genome's raw read in `rawread_folder` = ERR9786618:

```
|-- ERR9786618
  |-- A1.fastq
  |-- A2.fastq
```

adding `--batch` command, so Batch is `True`, (to analyze more than one we suggested this) Picota assumes the folder structure:

```
|-- rawread_folder
  |-- ERR9786618
    |-- ERR9786618_1.fastq
    |-- ERR9786618_2.fastq
  |-- SRR6242084
    |-- SRR6242084_1.fastq
    |-- SRR6242084_2.fastq
...
```


Example usage of the sra_download subcommand with default parameters

for Default (no Batch):
```
python picota.py assembly ERR9786618 out_folder_of_pipeline 
```
Batch:
```
python picota.py assembly --batch rawread_folder out_folder_of_pipeline 
```

Detailed usage:



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

Example usage of the analysis subcommand with default parameters

```
python picota.py analysis out_folder_of_pipeline/gfa_files out_folder_of_pipeline 
```

Detailed usage:


|__type__ |__command__ |__description__ |
| --- | --- | --- |
| required | `gfa_folder` | folder path of `.gfa` files |
| required | `out_folder` | out folder for the pipeline, result will be in `analysis` folder in out folder |
| optional | `--find_all_path` | `default=True`, finding all possible paths is best, but sometimes it can take too long running time, I recommend to try default unless you need very high speed |
| optional | `--path_limit` | `default=10`, High number of path limit will slow the analysis, I recommend a value between 10-25, because if the path is longer than 25, it is probably not a composite transposon unless the average node lenghts are very small |
| optional | `--min_size_of_cycle` | `default=1500`, min size for composite transposons (assumption in cyclic form) |
| optional | `--max_size_of_cycle` | `default=100000`, max size for composite transposons (assumption in cyclic form) |
| optional | `--name_prefix_cycle` | `default=''`, prefix of cycles in output fasta file |
| optional | `--min_component_number` | `default=1`, min size of components |
| optional | `--max_component_number` | `default=25`, min size of components |
| optional | `--k_mer_sim` | `default=200`, splitting k_mers when looking similarities, kind of ANI |
| optional | `--threshold_sim` | `default=80`, threshold in similarity elimination |


### 'db' module: Downloading Databases

Example usage of the db subcommand with default parameters

To download Antibiotics DB ([CARD](https://card.mcmaster.ca/)):

```
python picota.py db antibiotics 
```

To download Kegg Xenobiotic Metabolism DB ([KEGG](https://www.genome.jp/kegg/pathway.html)):

```
python picota.py db xenobiotics 
```
`xenobiotics_pathway_list` file (a file with KEGG ko numbers splitted by newline for xenobiotics pathway, come with default but you can update with adding more ko numbers from KEGG:

```
ko00362
ko00627
ko00364
ko00625
ko00361
```


To download IS Finder Insertion Sequence DB ([ISFinder](https://isfinder.biotoul.fr/)):

```
python picota.py db insertion_sequences 
```

Detailed usage:

|__type__ |__command__ |__description__ |
| --- | --- | --- |
| required | `name` | choices=[`antibiotics`, `xenobiotics`,`insertion_sequences`], select which db you want to download |
| optional | `--db_folder` | `default='DBs'`, folder of DBs, I recommend to use default naming, because scoring module will need this folder, you have to specify if needed in its input |
| optional | `--link_of_antibiotics` | `default='https://card.mcmaster.ca/latest/data'`, link of Antibiotics db to download |
| optional | `--link_of_is` | `default='https://github.com/thanhleviet/ISfinder-sequences/raw/master/IS.fna'`, link of ISs db to download |
| optional | `--input_pathway_list` | `default='xenobiotics_pathway_list'`, a file with KEGG ko numbers splitted by newline for xenobiotics pathway |
| optional | `--kegg_db_temp_folder` | `default='kegg_db_temp'`, temp folder for KEGG files |
| optional | `--metabolism_folder` | `default='Xenobiotics'`, metabolism_folder for KEGG files |
| optional | `--max_number_of_attends` | `default=50`, maximum number of attend for KEGG download |
| optional | `--antibiotics_folder` | `default='Antibiotics'`, antibiotics out folder |
| optional | `--model_type` | `default='nucleotide_fasta_protein_homolog_model'`, antibiotics model for CARD DB |
| optional | `--is_folder` | `default='InsertionSequences'`, IS folder |
| optional | `--is_file_name` | `default='IS.database.fasta'`, IS file name |


### 'scoring' module: Scoring Candidates

PICOTA needs [prodigal](https://github.com/hyattpd/Prodigal) to find CDS regions and [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) for (`blastn`) and (`madeblastdb`) programs for sequence search.

Example usage of the analysis subcommand with default parameters

```
python picota.py scoring out_folder_of_pipeline/analysis_folder out_folder_of_pipeline 
```

Detailed usage:

|__type__ |__command__ |__description__ |
| --- | --- | --- |
| required | `analysis_folder` | Folder path of cycle fasta files, output of analysis module |
| required | `out_folder` | out folder for the pipeline, the result will be in `picota_out` |
| optional | `--path_of_prodigal` | `default='prodigal'`, prodigal path if you dont have prodigal in PATH |
| optional | `--path_of_blastn` | `default='blastn'`, blastn path if you dont have blastn in PATH |
| optional | `--path_of_makeblastdb` | `default='makeblastdb'`, makeblastdb path if you dont have makeblastdb in PATH |
| optional | `--mean_of_comptns` | `default=12312`, Mean length of Composite Transposons |
| optional | `--std_of_comptns` | `default=1871`, Standard deviation of Composite Transposons |
| optional | `--total_score_type` | `default=0`, choices=[`0`, `1`,`2`] Can be 0,1 or 2. To select the main scoring system in the alimination process |
| optional | `--threshold_final_score` | `default=50`, Lower scores than `50` will be eliminated in final results |
| optional | `--max_z` | `default=20`, Maximum Z score, it is needed to calculate inverse normalization of z-scores |
| optional | `--dist_type` | `default=1`, choices=[`0`, `1`]  dist type 1 related with negative z scores is 0,  dist type 0 is normal dist of z scores, 1 is more accurate because lower lenghts means low number of genes, so score will be lower anyway |
| optional | `--splitted_cycle_folder` | `default='splitted_cycle'`, splitted cycles, it is temp folder for prodigal |
| optional | `--prodigal_out_folder` | `default='prodigal'`, prodigal out folder, it is temp folder |
| optional | `--out_blast_folder` | `default='blast_files'`, blast folder for temp blast files and dbs |
| optional | `--picota_out_folder` | `default='picota_out'`, picota out folder for temp picota files |
| optional | `--path_to_antibiotics` | `default='DBs/Antibiotics/nucleotide_fasta_protein_homolog_model.fasta'`, path of fasta files for antibitoics |
| optional | `--path_to_xenobiotics` | `default='DBs/Xenobiotics/allXenobiotics.fasta'`, path of fast files for xenobiotics |
| optional | `--path_to_ises` | `default='DBs/InsertionSequences/IS.database.fasta'`, path of fast files for insertion sequences |




### 'all' module: All In One Command

All module starts from sra download module to the end of the PICOTA pipeline. All the optional parameters mentioned above can be used in this module.

Example usage of the analysis subcommand with default parameters

```
python picota.py all sra_list.txt out_folder_of_pipeline 
```



## Method

### Composite Transposons

Transposons are the DNA sequences can mobilize and alter its position in genome, changing the genome size. Also, they can contribute to genome and gene evolution. Composite transposons are composed of two transposon and genetic material inside the flanking transposons. Since transposition event is somewhat stochastic, these transposons that are close to each other can act as one transposon and faciliate the transfer of the genetic material between.

We are in the brink of a "post-antibiotic" era. As more an more bacteria gain resistance to antibiotics due to misuse and overuse of antibiotics, it is more important now then ever to understand how bacteria transmit antibiotic resistance genes to each other if we don't want to go back to an era where even small scratches are deadly. Composite transposons also often carry antibiotic resistence genes, and sometimes metabolic genes such as degradation of xenobitiotics and metal resistance. Due to the nature of the transposons, these genes can be mobilized and potentially introduced to other bacteria. Therefore, it is crucial to identify and understand composite transposons if we want to understand the transfer of antibiotic resistance between different organisms and find interesting metabolic functions of genes that can be transferred.

There are many tools to find the transposons but there is a few tools related to composite transposons. One of them is [TnComp_finder](https://github.com/danillo-alvarenga/tncomp_finder). The limitation of TnComp_Finder is that it can be only used in complete genomes. For the incomplete genomes, it will probably fail because short-read sequencing (most of the time) is not able to capture composite transposons due to repetetive nature of these transposons. Because the genome databases are mostly consists of incomplete genomes, a tool that can work also in incomplete genomes is essential, so we developed PICOTA.

### Pipeline of PICOTA

To capture the composite transposons in incomplete genomes, assembly graphs built from the raw reads can be very helpful unlike using fasta sequences since the links that are between different parts of genomics sequences are lost in fasta files. Graph algorithms can be used to find transposon motifs, previously known transposons and possible novel composite transposons can be identified even *de-novo*. Then, the candidates can be investigated in existing gene databases for functional analysis. 

<p align="center">
<img src="logo/1_pipeline.png" alt="Picota logo" width="850"/>
</p>


### Assembly Graph Analysis


<p align="center">
<img src="logo/2_flowchart.png" alt="Picota logo" width="850"/>
</p>


### Scoring Functions

There are three different scoring as Scoring0, Scoring1, Scoring2. If Scoring2 == 0, then there is no insertion sequence in the candidate. If Score0 > Score1, then there is more than one interesting gene (such as antibiotics, xenobiotics etc.) in the candidate. I suggest the table below for the interpretation based on Score1: 

|__score__ |__interpretation__ |
| --- | --- | 
| Score1 = 0 | not enough |
| 0 < Score1 <= 50| low |
| 50 < Score1 <= 100 | moderate |
| 100 < Score1 <= 150 | high |
| 150 < Score1 <= 300 | perfect |



### Future Analysis


### Clustering

## FAQ

### How can I cite


## Licence
