# Simple Translatome Analysis Tool for Ribo-Seq (STATR) pipeline
Enviromental factors can have a huge impact on the condition of a cell. This condition may be reflected by the gene expression. Using Ribo-Seq this gene expression can be accuratly measured. Ribo-seq, also known as ribosome profiling, is a next-generation sequencing (NGS) technology, which can analyse the ribosome-protected mRNAs. Ribo-Seq has some major advantegous compared to similar techniques. In contrast to similar technique, Ribo-Seq analyses the messenger-RNA (mRNA) that is currently being translated on a ribosome. This provides a more accurate overview of the protein levels in a cell and can therefore be used to compare cells under different conditions.

In order to compare the Ribo-Seq data from different cells this pipeline was created. The pipeline can have multiple samples (with repeats) as input and creates the following outpus:
* Dendogram
* Heatmap
* Principal Component Analysis (PCA)

## Pipeline overview
![STATR pipeline directed acyclic graph](https://github.com/JobMaathuis/STATR_pipeline/blob/main/images/dag.png)

The function of each rules is explaned under the Rules section

## Repository structre
```
eindopdracht/
├── config/
│   ├── config.yaml
├── images/
│   └── dag.png
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
└── README.md
```
## Installation
The following installations need to be executed in order to use the pipeline:
### Snakemake
Snakemake can be installed by executing the following coomand:

`pip install snakemake`
### Trimmomatic
Open a terminal in the resources folder and execute the following commands:

`wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip`
`unzip Trimmomatic-0.39.zip`
### R packages
All of the R packages can be installed by executing the `InstallPacakges.R` script

## Usage
In order te use the pipeline the config file needs to be configured. The configfile can be found in `config/config.yaml`. Change the following paths in this file:

* wdir: path to your working direcotry
* files: path to your input files (relative to the working directory)
* samples: name of the different samples
* adapter: name of the adapter file (which should be placed in the same directory as the input files)
* genome: name of the NCBI reference sequence of the used genome (which should be placed in the same directory as the input files)

The pipeline can be used as follows:

`snakemake --snakefile {path to snakefile} --cores {amount of cores}`

## Rules
The following rules were used in the snakemake pipeline:

` rule trim_reads` 

Trims the input files where it reomoves the adapter and other illuminca-specific seqeunces.

`rule index_genome`

Indexes the genome using Bowtie2 (which uses FM index)

`rule align_reads`

Aligns the trimmed reads to the indexed genome, which creates a .sam file for every sample

`rule convert_sam_to_bam`

Converts the .sam files to .bam files, which is a binary format and advantageous for computer programs

`rule sort_bam`

Sorts the .bam files by leftmost coordinates using samtools

`rule convert_bam_to_bed`

In order to store the genomic regions by coordinates and annotation the.bam files are converted to .bed files using bedtools

`rule parse_genome_anntation`

...

`rule check_periodicity`

...

`rule generate profile`

...

`rule compute_coverage`

...

`rule format_deseq_input`

...

`rule run_deseq`

...




## Contact
Nog aanvullen met tekstje

j.maathuis@st.hanze.nl
