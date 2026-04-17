"""
tests/conftest.py
=================
Shared session-scoped fixtures for the scNoiseMeter test suite.

The tiny_gtf and tiny_bam fixtures create minimal synthetic data once per
test session. test_integration.py uses them to exercise the full
annotation-build → pipeline → metrics path without any external files.

GTF layout (chr1 only; chrM is excluded from annotation by design)
-------------------------------------------------------------------
GENE_A  chr1  +  1000–6000 (protein_coding)
  exon1:  1000–2000
  intron1: 2001–2999   (implicit; gene body minus exons)
  exon2:  3000–4000
  intron2: 4001–4999
  exon3:  5000–6000

GENE_B  chr1  −  15000–25000 (lncRNA)
  exon1:  15000–17000
  intron1: 17001–19999
  exon2:  20000–22000

Intergenic on chr1: positions < 999, 6000–14999, 25000+

BAM read inventory (all 100 bp, CIGAR = 100M)
----------------------------------------------
All positions are 0-based reference_start.

read_exonic_sense      chr1:1100  +   CB=cell1  NH=1  → EXONIC_SENSE
read_exonic_antisense  chr1:1300  −   CB=cell1  NH=1  → EXONIC_ANTISENSE
read_multimapper       chr1:1500  +   CB=cell1  NH=3  → MULTIMAPPER
read_chimeric          chr1:1700  +   CB=cell1  NH=1  SA chr1:50000 → CHIMERIC
read_intronic_boundary chr1:1980  +   CB=cell2  NH=1  → INTRONIC_BOUNDARY
  (20 exonic bases in exon1, 80 intronic in intron1; no N in CIGAR)
read_intronic_pure     chr1:2100  +   CB=cell2  NH=1  → INTRONIC_PURE
read_exonic_sense_neg  chr1:15100 −   CB=cell2  NH=1  → EXONIC_SENSE (GENE_B)
read_intergenic        chr1:30000 +   CB=cell2  NH=1  → INTERGENIC_*
read_mito              chrM:200   +   CB=cell1  NH=1  → MITOCHONDRIAL
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pysam
import pytest

BARCODE_CELL1 = "AAACCTGAGAAACCAT"
BARCODE_CELL2 = "AAACCTGAGAAACCGG"

_TINY_GTF = textwrap.dedent("""\
    chr1\tHAVANA\tgene\t1000\t6000\t.\t+\t.\tgene_id "ENSG00000000001"; gene_name "GENE_A"; gene_type "protein_coding";
    chr1\tHAVANA\ttranscript\t1000\t6000\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001"; gene_type "protein_coding";
    chr1\tHAVANA\texon\t1000\t2000\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "1"; gene_type "protein_coding";
    chr1\tHAVANA\texon\t3000\t4000\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "2"; gene_type "protein_coding";
    chr1\tHAVANA\texon\t5000\t6000\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST00000000001"; exon_number "3"; gene_type "protein_coding";
    chr1\tHAVANA\tgene\t15000\t25000\t.\t-\t.\tgene_id "ENSG00000000002"; gene_name "GENE_B"; gene_type "lncRNA";
    chr1\tHAVANA\ttranscript\t15000\t25000\t.\t-\t.\tgene_id "ENSG00000000002"; transcript_id "ENST00000000002"; gene_type "lncRNA";
    chr1\tHAVANA\texon\t15000\t17000\t.\t-\t.\tgene_id "ENSG00000000002"; transcript_id "ENST00000000002"; exon_number "1"; gene_type "lncRNA";
    chr1\tHAVANA\texon\t20000\t22000\t.\t-\t.\tgene_id "ENSG00000000002"; transcript_id "ENST00000000002"; exon_number "2"; gene_type "lncRNA";
""")


@pytest.fixture(scope="session")
def tiny_gtf(tmp_path_factory) -> Path:
    """Plain-text GTF with two genes on chr1."""
    tmp = tmp_path_factory.mktemp("annotation")
    p = tmp / "tiny.gtf"
    p.write_text(_TINY_GTF)
    return p


def _make_seg(header, name, flag, ref_id, pos, cb, ub, nh, *, sa=None):
    seg = pysam.AlignedSegment(header)
    seg.query_name = name
    seg.query_sequence = "ACGT" * 25          # 100 bp
    seg.flag = flag
    seg.reference_id = ref_id
    seg.reference_start = pos
    seg.mapping_quality = 60
    seg.cigar = [(0, 100)]                     # 100M
    seg.set_tag("CB", cb)
    seg.set_tag("UB", ub)
    seg.set_tag("NH", nh)
    if sa:
        seg.set_tag("SA", sa)
    return seg


@pytest.fixture(scope="session")
def tiny_bam(tmp_path_factory) -> Path:
    """
    Coordinate-sorted, indexed BAM with one read per major category.
    Returns the path to the sorted BAM (.bai index written alongside it).
    """
    tmp = tmp_path_factory.mktemp("bam")
    raw  = tmp / "tiny_unsorted.bam"
    out  = tmp / "tiny.bam"

    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            # Full GRCh38 chr1 length so _select_contigs does not filter it out
            {"SN": "chr1", "LN": 248_956_422},
            {"SN": "chrM", "LN": 16_569},
        ],
        # minimap2 @PG triggers ONT platform detection in inspect_bam
        "PG": [{"ID": "minimap2", "PN": "minimap2", "VN": "2.26"}],
    })

    FWD, REV = 0, 16  # SAM flag bits

    # Reads in coordinate order (chr1 then chrM; within chr1, ascending position).
    #
    # Cell-read budget:
    #   cell1: exonic_sense ×7, exonic_antisense, multimapper, chimeric, mito = 11 reads
    #   cell2: intronic_boundary, intronic_pure, exonic_sense_neg ×7, intergenic = 10 reads
    # Both cells reach the 10-read threshold required for the per-cell table.
    reads = [
        # --- categorical reads ---
        _make_seg(header, "read_exonic_sense",      FWD, 0, 1100,  BARCODE_CELL1, "AAAAAAAAAAAA", 1),
        _make_seg(header, "read_exonic_antisense",  REV, 0, 1300,  BARCODE_CELL1, "CCCCCCCCCCCC", 1),
        _make_seg(header, "read_multimapper",       FWD, 0, 1500,  BARCODE_CELL1, "GGGGGGGGGGGG", 3),
        _make_seg(header, "read_chimeric",          FWD, 0, 1700,  BARCODE_CELL1, "TTTTTTTTTTTT", 1,
                  sa="chr1,50000,+,100M,60,0;"),
        _make_seg(header, "read_intronic_boundary", FWD, 0, 1980,  BARCODE_CELL2, "ACACACACACAC", 1),
        _make_seg(header, "read_intronic_pure",     FWD, 0, 2100,  BARCODE_CELL2, "TGTGTGTGTGTG", 1),
        _make_seg(header, "read_exonic_sense_neg",  REV, 0, 15100, BARCODE_CELL2, "GCGCGCGCGCGC", 1),
        _make_seg(header, "read_intergenic",        FWD, 0, 30000, BARCODE_CELL2, "ATATATATATAT", 1),
        _make_seg(header, "read_mito",              FWD, 1, 200,   BARCODE_CELL1, "CGCGCGCGCGCG", 1),
        # --- padding reads to reach the 10-read/cell threshold ---
        # cell1: 6 more exonic_sense reads (total cell1 = 11)
        _make_seg(header, "read_c1_pad1", FWD, 0, 1101, BARCODE_CELL1, "AAACAAACAAAC", 1),
        _make_seg(header, "read_c1_pad2", FWD, 0, 1102, BARCODE_CELL1, "AAAGAAAGAAAG", 1),
        _make_seg(header, "read_c1_pad3", FWD, 0, 1103, BARCODE_CELL1, "AACCAACCAACC", 1),
        _make_seg(header, "read_c1_pad4", FWD, 0, 1104, BARCODE_CELL1, "AAGCAAGCAAGC", 1),
        _make_seg(header, "read_c1_pad5", FWD, 0, 1105, BARCODE_CELL1, "AATCAATCAATC", 1),
        _make_seg(header, "read_c1_pad6", FWD, 0, 1106, BARCODE_CELL1, "ACACACACAAAC", 1),
        # cell2: 6 more exonic_sense reads on GENE_B -strand (total cell2 = 10)
        _make_seg(header, "read_c2_pad1", REV, 0, 15101, BARCODE_CELL2, "CAAACAAACAAA", 1),
        _make_seg(header, "read_c2_pad2", REV, 0, 15102, BARCODE_CELL2, "CAAGCAAGCAAG", 1),
        _make_seg(header, "read_c2_pad3", REV, 0, 15103, BARCODE_CELL2, "CAACCAACCAAC", 1),
        _make_seg(header, "read_c2_pad4", REV, 0, 15104, BARCODE_CELL2, "CAATCAATCAAT", 1),
        _make_seg(header, "read_c2_pad5", REV, 0, 15105, BARCODE_CELL2, "CACACACACAAC", 1),
        _make_seg(header, "read_c2_pad6", REV, 0, 15106, BARCODE_CELL2, "CGAACGAACGAA", 1),
    ]

    with pysam.AlignmentFile(str(raw), "wb", header=header) as bam:
        for r in reads:
            bam.write(r)

    pysam.sort("-o", str(out), str(raw))
    pysam.index(str(out))
    return out
