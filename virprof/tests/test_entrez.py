"""Test cases for features module"""

from ..entrez import FeatureTableParser

# FIXME: multiple db_xref entries probably possible

feature_table_raw = "\n".join(
    (
        ">Feature ref|NC_045512.2|",
        # feature with no meta data
        "1	265	5'UTR",
        "266	21555	gene",
        "			gene	ORF1ab",
        "			locus_tag	GU280_gp01",
        "			db_xref	GeneID:43740578",
        # joined region:
        "266	13468	CDS",
        "13468	21555",
        "			product	ORF1ab polyprotein",
        "			protein_id	ref|YP_009724389.1|",
        "			note	pp1ab; translated by -1 ribosomal frameshift",
        "			ribosomal_slippage",
        "			exception	ribosomal slippage",
        "266	13483	CDS",
        "			product	ORF",
        "29675	29903	3'UTR",
        # second feature:
        ">Feature other_source|accession|comment",
        # feature with no meta data and extending before sequence
        "<501	600	5'UTR",
        # gene on reverse strand
        "800	700	gene",
        "			gene	reverse_strand",
        # feature extending beyond end
        "900	>1000	3'UTR",
    )
)

feature_table_expected = {
    "NC_045512.2": [
        (1, 265, "5'UTR", {}),
        (
            266,
            21555,
            "gene",
            {
                "gene": "ORF1ab",
                "locus_tag": "GU280_gp01",
                "db_xref": "GeneID:43740578",
            },
        ),
        (
            266,
            13468,
            "CDS",
            {
                # "_multipart": True,
                "product": "ORF1ab polyprotein",
                "protein_id": "ref|YP_009724389.1|",
                "note": "pp1ab; translated by -1 ribosomal frameshift",
                "ribosomal_slippage": "",
                "exception": "ribosomal slippage",
            },
        ),
        (
            13468,
            21555,
            "CDS",
            {
                # "_multipart": True,
                "product": "ORF1ab polyprotein",
                "protein_id": "ref|YP_009724389.1|",
                "note": "pp1ab; translated by -1 ribosomal frameshift",
                "ribosomal_slippage": "",
                "exception": "ribosomal slippage",
            },
        ),
        (266, 13483, "CDS", {"product": "ORF"}),
        (29675, 29903, "3'UTR", {}),
    ],
    "accession": [
        (501, 600, "5'UTR", {"_left_open": True}),
        (700, 800, "gene", {"_reversed": True, "gene": "reverse_strand"}),
        (900, 1000, "3'UTR", {"_right_open": True}),
    ],
}


def test_feature_parser():
    parser = FeatureTableParser()
    result = parser.parse(feature_table_raw)
    assert result.keys() == feature_table_expected.keys()
    for res, exp in zip(result.values(), feature_table_expected.values()):
        assert res == exp

    assert result == feature_table_expected
