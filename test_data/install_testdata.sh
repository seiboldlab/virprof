#!/bin/bash
set -e

# Create basic ymp.yml
cat >ymp.yml <<EOF
include:
  - virprof/virprof.yml
resource_limits:
  mem:
    max: 6G
EOF
echo "Created ymp.yml:"
echo "---"
cat ymp.yml
echo "---"
echo


unpack() {
    src=virprof/test_data/$1
    dst=databases/$2
    mkdir -p $dst
    echo "Unpacking $src to $dst"
    tar xfv $src -C $dst
    echo "DONE"
    echo
}

unpack Homo_sapiens_UCSC_hg38_test.tar.bz2
unpack grch38_snp_tran_test.tar.bz2
unpack taxdump_test.tar.bz2 NCBI_taxonomy
unpack nt_test.tar.bz2 nt
