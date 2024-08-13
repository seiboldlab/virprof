# CHANGELOG

#### v0.9.0
 - feature: add FastA QC metrics entropy, homopolymer fraction
 - feature: extended QC filtering (fixes #60)
 - feature: make whitelist and filtering parameters configurable (fixes #55)

#### v0.8.2
 - use https to fetch VHDB instead of ftp
 - fix fastqc (update, hopefuly no more crashes)
 - feature: add umitools pipeline
 - feature: improve `virprof prepare-phylo` command

#### v0.8.1
 - fix RDS not copied to results folder
 
#### v0.8.0
 - breaking change: rename empty_samples to failed samples
 - feature: preclassify contigs with LCA
 - feature: add run information to exported reports
 - feature: add a sample sheet tab to report
 - feature: add a RDS version of the Excel report
 - fix: `num_threads` must be smaller than `num_samples`
 - fix: unit test setup
 - fix: use black formatting

#### v0.7.0
 - improve metadata(se)$fastqc: now includes the fastq_file_path for
   the respective mate as well as the remainder of the columns from
   the primary project table (for later grouping).
 - add entire sample sheet as metadata(se)$sample_sheet
 - fill assay matrix with zero or fake length for samples that failed
   salmon quantification.
 - parallelize parsing of salmon files (speed up big projecs)

#### v0.6.1
 - store coldata object in metadata
 - fix PCA-for-QC calculation breaks on zeros; now uses type=poscounts
   for estimating size factors.
 - fix pct_* variables containing range 0 to 1 (not 0 - 100)

#### v0.6.0
 - add VP S4 object carrying pathogen data
 - fix error on empty (after qc, after human deplete) fastq
 - fix numeric sample names break merge
 - add named report export to pathogen pipeline
 - fix raw depleted reads kept unnecessarily (space optimize)
 - add incorporation of QC stats into gene RDS
 - add export of metadata only RDS

#### v0.5.0
 - add workflow for extracting exon counts (with R GenomicAlignments)
 - add workflows for handling UMIs on bulk RNA-seq data
 - convert Gencode reference fasta.gz to bgzip with index
 - fix ulimit issue with large number of genomes detected

#### v0.4.0
 - changed RNA-seq pipelines to use STAR "twopassMode Basic"
 - add stage for extracing junctions from BAM with regtools
 - don't be silent while posting to Slack (to see network outage issues)

#### v0.3.1
 - add versioning to output RDS
 - add salmon checks (identical version, correction flags, reference hash/decoys)
 - update integration tests

#### v0.3.0
 - add generic unmapped read extract from rnaseq pipeline
 - add TRF stage for filtering tandem repeats (TR)
 - add TR filtering to pathogen pipeline
 - add pathogen pipeline (same as virus, works after rnaseq)
 - make reporting stage detect pipeline used internally
 - remove old hg38 reference
 - make new hg38g reference work with virus/pathogen pipeline
 - add mouse reference
 - remove one-of references (covid, hrvc, palmenberg hrv)
 - add Salmon sample metadata to output RDS (unique to metadata, different to colData)
 - use standard YMP stage for making blast index on reference genome
 - have pathogen pipeline report same way as rnaseq pipelines

#### v0.2.0
 - use salmon index with full decoy
 - add our own tximport stage

#### v0.1.0:
 - begin versioning
