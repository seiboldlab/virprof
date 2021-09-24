#!/bin/bash
set -e

cp virprof/test_data/ymp.yml .

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
unpack nt_test.tar.bz2 nt
