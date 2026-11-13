"""
tests/test_sense_intron_precedence.py
=====================================
A base inside an intron of a gene on the read's own strand, which is also
inside an exon of a gene on the other strand, counts as intronic.

GTF layout (chr1)
-----------------
HOST    +  1000-9000   exons 1000-2000 and 8000-9000, intron 2001-7999
NESTED  -  4000-5000   single exon, inside the HOST intron
ALONE   -  20000-21000 single exon, no gene on the other strand
"""

from __future__ import annotations

import textwrap

import pysam
import pytest

from scnoisemeter.constants import ReadCategory
from scnoisemeter.modules.annotation import build_annotation_index
from scnoisemeter.modules.classifier import ReadClassifier

_GTF = textwrap.dedent("""\
    chr1\tHAVANA\tgene\t1000\t9000\t.\t+\t.\tgene_id "ENSG00000000001"; gene_name "HOST"; gene_type "protein_coding";
    chr1\tHAVANA\ttranscript\t1000\t9000\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001"; gene_type "protein_coding";
    chr1\tHAVANA\texon\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "1"; gene_type "protein_coding";
    chr1\tHAVANA\texon\t8000\t9000\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "2"; gene_type "protein_coding";
    chr1\tHAVANA\tgene\t4000\t5000\t.\t-\t.\tgene_id "ENSG00000000002"; gene_name "NESTED"; gene_type "lncRNA";
    chr1\tHAVANA\ttranscript\t4000\t5000\t.\t-\t.\tgene_id "ENSG00000000002"; transcript_id "ENST00000000002"; gene_type "lncRNA";
    chr1\tHAVANA\texon\t4000\t5000\t.\t-\t.\tgene_id "ENSG00000000002"; transcript_id "ENST00000000002"; exon_number "1"; gene_type "lncRNA";
    chr1\tHAVANA\tgene\t20000\t21000\t.\t-\t.\tgene_id "ENSG00000000003"; gene_name "ALONE"; gene_type "lncRNA";
    chr1\tHAVANA\ttranscript\t20000\t21000\t.\t-\t.\tgene_id "ENSG00000000003"; transcript_id "ENST00000000003"; gene_type "lncRNA";
    chr1\tHAVANA\texon\t20000\t21000\t.\t-\t.\tgene_id "ENSG00000000003"; transcript_id "ENST00000000003"; exon_number "1"; gene_type "lncRNA";
""")

_HEADER = pysam.AlignmentHeader.from_dict({
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 248_956_422}],
})


@pytest.fixture(scope="module")
def classifier(tmp_path_factory):
    gtf = tmp_path_factory.mktemp("nested") / "nested.gtf"
    gtf.write_text(_GTF)
    return ReadClassifier(build_annotation_index(gtf, cache=False), reference=None)


def _read(pos: int, *, reverse: bool) -> pysam.AlignedSegment:
    seg = pysam.AlignedSegment(_HEADER)
    seg.query_name = f"read_{pos}_{'rev' if reverse else 'fwd'}"
    seg.query_sequence = "ACGT" * 25
    seg.flag = 16 if reverse else 0
    seg.reference_id = 0
    seg.reference_start = pos
    seg.mapping_quality = 60
    seg.cigar = [(0, 100)]
    seg.set_tag("NH", 1)
    return seg


def test_sense_intron_wins_over_antisense_exon(classifier):
    result = classifier.classify(_read(4400, reverse=False))
    assert result.category is ReadCategory.INTRONIC_PURE
    assert result.base_counts == {ReadCategory.INTRONIC_PURE: 100}


def test_read_entering_the_nested_exon_stays_intronic(classifier):
    # 49 bases of HOST intron alone, then 51 inside the NESTED exon as well.
    result = classifier.classify(_read(3950, reverse=False))
    assert result.base_counts == {ReadCategory.INTRONIC_PURE: 100}


def test_nested_gene_keeps_its_own_sense_reads(classifier):
    result = classifier.classify(_read(4400, reverse=True))
    assert result.category is ReadCategory.EXONIC_SENSE


def test_antisense_exon_without_a_sense_gene_is_still_antisense(classifier):
    result = classifier.classify(_read(20400, reverse=False))
    assert result.category is ReadCategory.EXONIC_ANTISENSE
