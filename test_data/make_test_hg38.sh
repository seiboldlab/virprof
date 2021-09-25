#!/bin/bash
set -e

ORIG_GTF=references/hg38g/ALL.gtf
ORIG_FA=references/hg38g/ALL.fasta.gz
ORIG_TX=references/hg38g/ALL.tx.fasta.gz

OUT_GTF=virprof/test_data/test.gtf.gz
OUT_FA=virprof/test_data/test.fasta.gz
OUT_TX=virprof/test_data/test.tx.fasta.gz

MAX_BP=1200000

WORK=$(mktemp -d)
trap "rm -rf $WORK" EXIT

# Truncate GTF at boundary
awk '/chr1/ {if ($5>'$MAX_BP') {exit}; print}' $ORIG_GTF |\
    gzip -c > $OUT_GTF

# Truncate FASTA at boundary
gunzip -dc $ORIG_FA |\
    awk '/^>/ {print} /^[^>]/ {print; total=total+length($0); if (total>'$MAX_BP') {exit}}' |\
    gzip -c > $OUT_FA

# Extract gene_id
gunzip -dc $OUT_GTF |\
    awk '{print $10}' |\
    sort | uniq |\
    sed 's/^"//; s/";$//'> $WORK/gene_ids.txt

# Filter TX FASTA
gunzip -dc $ORIG_TX |\
    awk '
BEGIN {
  while ((getline line < "'$WORK/gene_ids.txt'") > 0) {
    geneids[line] = 1
  }
  FS="|"
  found=0
}
/^>/ {
  found=0
  if (geneids[$2]) {
    print;
    found=1	
  }
}
found && $0 ~ /^[^>]/ {
  print
}
' |\
    gzip -c > $OUT_TX

ymp make refhg38gtest.index_rsem
ymp env run rsem -- \
    rsem-simulate-reads \
    ref_hg38gtest.index_rsem/ALL.rsem \
    virprof/test_data/HR5366.model \
    virprof/test_data/HR5366.isoforms.results \
    0.5 1000 \
    virprof/test_data/sim




head $WORK/gene_ids.txt
tail $WORK/gene_ids.txt



# Extract gene_id for beginning of genome
#awk '/chr1/ {if ($5>'$MAX_BP') {exit; }; print $10}' $ORIG_GTF | uniq |  > $WORK/gene_ids.txt



ls -latr $WORK
