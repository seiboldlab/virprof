#!/bin/bash
set -e

SRC=${BASH_SOURCE[0]}
while [ -L "$SRC" ]; do
  DIR=$(cd -P "$(dirname "$SRC")" >/dev/null 2>&1 && pwd)
  SOURCE=$(readlink "$SRC")
  [[ $SRE != /* ]] && SRC=$DIR/$SRC
done
DIR=$(cd -P "$(dirname "$SRC")" >/dev/null 2>&1 && pwd)

cp "$DIR/ymp.yml" .

unpack() {
    src="$DIR/$1"
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
