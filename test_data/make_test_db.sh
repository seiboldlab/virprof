#!/bin/bash

SRCSTAGES=test.qc.ref_hg38.deplete.assemble_meta.coverage.ref_hg38.annotate_blastE2MegaBest.blastfilter_vpU200.ref_NT.annotate_blast*
zgrep -v "^#" $SRCSTAGES/*.blast7.gz | cut -f2 -d $'\t' | sort| uniq > expected_accs.txt

blastdbcmd -db databases/nt/nt -entry_batch expected_accs.txt -out expected_subjects.fasta

blastdbcmd -db databases/nt/nt -entry_batch expected_accs.txt -out expected_subjects.taxid_map -outfmt "%a %T"

makeblastdb -in expected_subjects.fasta -parse_seqids -blastdb_version 5 -taxid_map expected_subjects.taxid_map -title "NT subset for testing"  -dbtype nucl -out databases/nt_test/nt
