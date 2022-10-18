#!/bin/bash
set -e

# Find path to self (symlink acrobatics)
SRC=${BASH_SOURCE[0]}
while [ -L "$SRC" ]; do
  DIR=$(cd -P "$(dirname "$SRC")" >/dev/null 2>&1 && pwd)
  SOURCE=$(readlink "$SRC")
  [[ $SRE != /* ]] && SRC=$DIR/$SRC
done
# Self is in test_data directory
TESTDATA=$(cd -P "$(dirname "$SRC")" >/dev/null 2>&1 && pwd)
# which is in base directory
BASE=$(cd -P "$(dirname "$TESTDATA")" >/dev/null 2>&1 && pwd)
# get relative path to that
RELBASE=$(python -c "import os.path; print(os.path.relpath('$BASE'))")

sed "s%virprof/virprof.yml%$RELBASE/virprof.yml%" "$TESTDATA/ymp.yml"  > ymp.yml

unpack() {
    src="$TESTDATA/$1"
    dst=databases/$2
    mkdir -p $dst
    echo "Unpacking $src to $dst" >&2
    tar xfv $src -C $dst
    echo "DONE" >&2
    echo >&2
}

unpack Homo_sapiens_UCSC_hg38_test.tar.bz2
unpack grch38_snp_tran_test.tar.bz2
unpack nt_test.tar.bz2 nt

mkdir -p test_data
tail -n +2 $TESTDATA/test.csv |
    while IFS=, read unit sample fq1 fq2 virus; do
	[ -z "$unit" ] && continue
	for n in 1 2; do
	    fq=fq$n
	    echo making test_data/${!fq} >&2
	    (
		echo "+ $TESTDATA/sim${unit}_$n.fq.gz" >&2
		gzip -dc $TESTDATA/sim${unit}_$n.fq.gz |\
		    sed 's/^@/@sim'${unit}'_/'
		echo -n $virus+ | while read -d + vir; do
		    [ -z "$vir" ] && continue
		    echo "+ $TESTDATA/test_${vir}.R.$n.fq.gz" >&2
		    gzip -dc $TESTDATA/test_$vir.R$n.fq.gz |\
			sed 's/^@/@vir_'${vir}'_/'
		done
	    ) | gzip -c > test_data/${!fq} 
	done 
    done
cp $TESTDATA/test.csv test_data
