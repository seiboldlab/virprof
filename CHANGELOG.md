CHANGELOG
---------
v0.5.0:
 - added pipeline for extracting exon counts (with R GenomicAlignments)

v0.4.0:
 - changed RNA-seq pipelines to use STAR "twopassMode Basic"
 - add stage for extracing junctions from BAM with regtools
 - don't be silent while posting to Slack (to see network outage issues)

v0.3.1:
 - add versioning to output RDS
 - add salmon checks (identical version, correction flags, reference hash/decoys)
 - update integration tests

v0.3.0:
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

v0.2.0:
 - use salmon index with full decoy
 - add our own tximport stage

v0.1.0:
 - begin versioning
