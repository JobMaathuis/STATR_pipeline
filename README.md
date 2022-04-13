# Simple Translatome Analysis Tool for Ribo-Seq (STATR) pipeline
Environmental factors can have a huge impact on the condition of a cell. This condition may be reflected by the gene expression. Using Ribo-Seq this gene expression can be accurately measured. Ribo-seq, also known as ribosome profiling, is a next-generation sequencing (NGS) technology, which can analyse the ribosome-protected mRNAs. Ribo-Seq has some major advanteges compared to similar techniques. In contrast to similar technique, Ribo-Seq analyses the messenger-RNA (mRNA) that is currently being translated on a ribosome, also known as the translatome, instead of all mRNA in the cell. This provides a more accurate overview of the protein levels in a cell and can therefore be used to compare cells under different conditions.

In order to compare the Ribo-Seq data from different cells this pipeline was created. The pipeline can have multiple samples (with repeats) as input and creates the following outputs:

* Genome-wide Ribo-Seq profile (can be visualised with free software like Integrative Genomic Viewer)
* Ribosome densities file
* Principal Component Analysis (PCA) plot
* Heatmap

Example outputs of the PCA plot and heatmap can be found in the result directory

## Pipeline overview
STATR pipeline directed acyclic graph (dag)

![STATR pipeline directed acyclic graph](https://github.com/JobMaathuis/STATR_pipeline/blob/main/images/dag.png)

The function of each rules is explained under the Rules section

## Repository structure

```
STATR_pipeline/
├── config/
│   ├── config.yaml
├── images/
│   └── dag.png
├── logs/
├── resources/
│    ├── scripts/
│    │   ├── Python_scripts/
│    │   │   ├── CheckPeriodicity.py
│    │   │ 	 ├── FormatDESeqInput.py
│    │   │   ├── GenerateProfile.py
│    │   │   └── ParseGenomeAnnotation.py
│    │   └── R_scripts/
│    │       ├── InstallPacakges.R
│    │       └── RunDESeq.R
│    └── files/
│        ├── Design_sheet.txt
│        └── RiboSeq_adapter_as.fa
├── result/    
├── workflow/
│   ├── rules/
│   │   ├── align_reads.smk
│   │   ├── create_ribosome_profile.smk
│   │   ├── decompile_alignment.smk
│   │   ├── find_degs.smk
│   │   └── trim_reads.smk
│   └── main.smk
├── LICENSE
└── README.md

```


## Rules
The following rules were used in the snakemake pipeline:

` rule trim_reads` 

Trims the input files where it removes the adapter and other illumina-specific sequences.

`rule index_genome`

Indexes the genome using Bowtie2 (which uses FM index)

`rule align_reads`

Aligns the trimmed reads to the indexed genome, which creates a .sam file for every sample

`rule convert_sam_to_bam`

Converts the .sam files to .bam files, which is a binary format and advantageous for computer programs

`rule sort_bam`

Sorts the .bam files by leftmost coordinates using samtools

`rule convert_bam_to_bed`

In order to store the genomic regions by coordinates and annotation the .bam files are converted to .bed files using bedtools

`rule parse_genome_anntation`

Parses the genome annotation of the file to a readable and clean format (GFF3)

`rule calculate_ribosome_density`

Exploits the periodicity in a way that the average ribosome density over a coding DNA sequence can be calculated

`rule generate profile`

Generates a genome-wide Ribo-Seq profile

`rule compute_coverage`

Computes the coverage of each of the given genome features

`rule format_deseq_input`

Combines all of the samples in one file in a format that can be used by DESeq2

`rule run_deseq`

Runs a DESeq2 analysis on the samples, which returns a PCA plot and a heatmap as final output

`rule move_txt_files`

Moves the less interesting .txt files to its corresponding directory

## Installation
The following installations need to be executed in order to use the pipeline:

### Snakemake
Snakemake can be installed by executing the following command:

`pip install snakemake`

### Trimmomatic
Open a terminal in the resources folder and execute the following commands:

`wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip`
`unzip Trimmomatic-0.39.zip`

### Bioconda
The following packages are installed using Bioconda, which should be first installed:

`curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

`sh Miniconda3-latest-Linux-x86_64.sh`

`conda config --add channels defaults`

`conda config --add channels bioconda`

`conda config --add channels conda-forge`

#### Bowtie2
Bowtie2 can be installed as follows:

`conda install -c bioconda bowtie2`

#### Samtools 
Samtools can be installed as follows:

`conda install -c bioconda samtools`

#### Bedtools 
Samtools can be installed as follows:

`conda install -c bioconda bedtools` 


### R packages
All of the R packages can be installed by executing the `InstallPackages.R` script

## Usage
Firstly, the files listed below are needed. Put these in the `resources/files` directory

* at least two input files (gzipped Ribo-Seq fastq files)
* reference genome, in both .fa and .gff3 format
* Ribo-Seq adapter file in .fa format
* Design sheet (edit example in the files directory)

Secondly, the config file needs to be configured. The config file can be found in `config/config.yaml`. Change the following paths in this file:

* wdir: path to your working directory
* samples: name of the different samples
* adapter: name of the adapter file 
* genome: name of the genome file

The pipeline can be used as follows:

`snakemake --snakefile {path to snakefile} --cores {amount of cores}`

## Contact
The owner of the repository can be contacted by emailing to:

j.maathuis@st.hanze.nl
