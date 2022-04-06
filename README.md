# Simple Translatome Analysis Tool for Ribo-Seq (STATR) pipeline
Enviromental factors can have a huge impact on the condition of a cell. This condition may be reflected by the gene expression. Using Ribo-Seq this gene expression can be accuratly measured. Ribo-seq, also known as ribosome profiling, is a next-generation sequencing (NGS) technology, which can analyse the ribosome-protected mRNAs. Ribo-Seq has some major advantegous compared to similar techniques. In contrast to similar technique, Ribo-Seq analyses the messenger-RNA (mRNA) that is currently being translated on a ribosome. This provides a more accurate overview of the protein levels in a cell and can therefore be used to compare cells under different conditions.

In order to compare the Ribo-Seq data from different cells this pipeline was created. The pipeline can have multiple samples (with repeats) as input and creates the following outpus:
* Dendogram
* Heatmap
* Principal Component Analysis (PCA)

## Pipeline overview
![]()

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
The pipeline can be used as follows:

`snakemake --snakefile {path to snakefile} --cores {amount of cores}`

## Contact
Nog aanvullen met tekstje

j.maathuis@st.hanze.nl
