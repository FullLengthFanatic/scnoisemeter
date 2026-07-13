"""Regression tests for explicit multi-mapper evidence.

MAPQ 0 is not, by itself, proof of multiple alignments.  In particular, some
long-read BAM producers emit MAPQ 0 for all primary records.  Treating it as a
multi-mapper criterion collapses every genomic category into MULTIMAPPER.
"""

from types import SimpleNamespace

from scnoisemeter.modules.classifier import ReadClassifier


class _Read(SimpleNamespace):
    def has_tag(self, tag):
        return tag in self.tags

    def get_tag(self, tag):
        return self.tags[tag]


def _read(*, mapq=60, **tags):
    return _Read(mapping_quality=mapq, tags=tags)


def test_mapq_zero_without_multi_hit_tag_is_not_multimapper():
    """Protect PacBio/Kinnex-style BAMs whose primary records use MAPQ 0."""
    assert not ReadClassifier._is_multimapper(_read(mapq=0))


def test_low_nonzero_mapq_without_multi_hit_tag_is_not_multimapper():
    assert not ReadClassifier._is_multimapper(_read(mapq=1))


def test_nh_greater_than_one_is_multimapper():
    assert ReadClassifier._is_multimapper(_read(NH=2))


def test_nh_one_is_not_multimapper_even_at_mapq_zero():
    assert not ReadClassifier._is_multimapper(_read(mapq=0, NH=1))


def test_nonempty_xa_is_multimapper_when_nh_is_absent():
    assert ReadClassifier._is_multimapper(_read(XA="chr2,+100,100M,1;"))


def test_empty_xa_is_not_multimapper():
    assert not ReadClassifier._is_multimapper(_read(mapq=0, XA=""))
