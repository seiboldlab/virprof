[![Unit Tests](https://github.com/seiboldlab/virprof/actions/workflows/main.yml/badge.svg)](https://github.com/seiboldlab/virprof/actions/workflows/main.yml)

Virprof - rest of title here
============================

Virprof identifies viral infections and commensal microbes in RNA-Seq (and metagenomic) NGS data. It produces full-length phylogeny-quality viral genomes for samples with sufficient viral load and 

I. Installing Virprof
---------------------

1. Install Bioconda

   If you don't already have a working Bioconda installation, see [Installing Bioconda](https://bioconda.github.io/user/install.html#install-conda) for instructions. We recommend that you immediately install `mamba` (`conda install -n base mamba`) and use it to install software packages instead of `conda`. 

2. Install Virprof

   There is no package for Virprof yet (working on it). To install the development version directly from github, create a clean directory in which you want to run your analyses and run the following commands:

   ```bash
   # Clone the repository
   mamba install --name base git # unless you already have git installed
   git clone https://github.com/seiboldlab/virprof.git
   # Install YMP
   mamba create --name ymp --file virprof/environment.yaml
   conda activate ymp
   # Copy basic configuration
   cp virprof/ymp.yml .
   ```

3. Install Reference Data

   Most reference data needed by Virprof is installed automatically. Only the BLAST NT database and the host reference need to be installed manually because they are very large and likely to already be available on your compute infrastructure somewhere.
   
   Create a folder `databases/` and populate it with the following subdirectories (which may be symbolic links pointing to where you really have those files stored):
   - `nt`: This is the BLAST NT database from NCBI. The easiest way to download it is using the `update_blastdb.pl` script that is part of the `blast` package:

     ```
     mkdir -p databases/nt
     cd databases/nt
     ymp env run blast -- update_blastdb.pl --decompress nt
     cd ../..
     ```
     
     The download will be about 70GB and about 100GB will be needed for the unpacked files.
     
   - `Homo_sapiens`: These are the human genome reference files, including the Bowtie2 index, from the Illumina iGenomes.

     ```
     mkdir -p databases
     cd databases
     wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Homo_sapiens/UCSC/hg38/Homo_sapiens_UCSC_hg38.tar.gz
     tar xf Homo_sapiens_UCSC_hg38.tar.gz
     rm Homo_sapiens_UCSC_hg38.tar.gz
     cd ..
     ```
     
     The download will be about 15GB and about 24GB needed for the unpacked files.
     
   - `grch38_snp_tran`: The index files for HISAT2:

     ```
     mkdir databases
     cd databases
     wget http://genome-idx.s3.amazonaws.com/hisat/grch38_snptran.tar.gz
     tar xf grch38_snptran.tar.gz
     rm grch38_snptran.tar.gz
     cd ..
     ```
     
     The download will be about 5GB and about 7GB needed for the unpacked files.

II. Configuring Samples
-----------------------

1. Prepare a sample sheet in CSV, TSV or Excel format listing your samples.

   Each row should have a sample name and the locations of the forward and reverse reads. The name of the column with the sample names must be `sample`, the other columns can have arbitrary names. You can add as many other columns as you like:
   
   |sample|fq1|fq2|comment|
   |-|-|-|-|
   |A|A.fwd.fastq.gz|A.rev.fastq.gz|some comment|
   |B|B.fwd.fastq.gz|B.rev.fastq.gz||
   
   If you have samples that were sequenced multiple times, prepend a column (e.g. `unit`) containing unique names to each sequencing unit. 
   
   |unit|sample|fq1|fq2|
   |-|-|-|-|
   |A|A|00000001_S7_L002_R1_001.fastq.gz|00000001_S7_L002_R2_001.fastq.gz|
   |B|B|00000002_S1_L002_R1_001.fastq.gz|00000002_S1_L002_R2_001.fastq.gz|
   |A2|A|other_run/00000123_S29_L001_R1_001.fastq.gz|other_run/00000123_S29_L001_R2_001.fastq.gz|
   
   Paths are relative to the location of the sample sheet file. If you want to move it around later, use absolute paths. 
   
   You can also list SRA run identifiers directly:
   
   |sample|srr|
   |-|-|
   |A|SRR123123123|
   |B|ERR456456456|
   
2. Create a YMP project for your samples

   To add this project to your Virprof configuration, add the following to the `ymp.yml` file:

   ```yaml
   projects:
     name_of_my_project:
        data: path_to_my_sample_sheet.csv
   ```

III. Running Virprof
--------------------

Virprof was implemented with YMP. For all the details, refer to their [documentation](https://ymp.readthedocs.io).

To run the generic workflow on a project you configured above on the local machine, just issue:

```
ymp make -j64 name_of_my_project.generic
```

(Replace `64` with the number of threads you want to use).

To run the pipeline on a local SLURM cluster, just run `ymp init cluster` once to create the basic configuration. Then run 

```
ymp submit name_of_my_project.generic
```

IV. Citing Virprof
------------------

TBD
