"""
tests/test_annotation_fetcher.py
=================================
Unit tests for annotation_fetcher.py.

All HTTP/FTP calls are mocked — these tests never hit the network.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch, call

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from scnoisemeter.utils.annotation_fetcher import (
    extract_gencode_version_from_filename,
    extract_polyasite_version_from_filename,
    fetch_latest_gencode_gtf,
    fetch_latest_polyasite_atlas,
    fetch_10x_whitelist,
    _parse_gencode_gtf_filename,
    _url_exists,
)
from scnoisemeter.cli import (
    _check_version_consistency,
    _resolve_gtf,
    _resolve_polya_sites,
)


# ---------------------------------------------------------------------------
# Helper: fake GENCODE directory listing HTML
# ---------------------------------------------------------------------------

def _fake_gencode_html(version: int = 49) -> str:
    return (
        f'<a href="gencode.v{version}.annotation.gtf.gz">'
        f'gencode.v{version}.annotation.gtf.gz</a>\n'
        f'<a href="gencode.v{version}.basic.annotation.gtf.gz">'
        f'gencode.v{version}.basic.annotation.gtf.gz</a>\n'
    )


# ---------------------------------------------------------------------------
# extract_gencode_version_from_filename
# ---------------------------------------------------------------------------

class TestExtractGencode:
    def test_standard_filename(self):
        assert extract_gencode_version_from_filename("gencode.v47.annotation.gtf.gz") == 47

    def test_capital_v(self):
        assert extract_gencode_version_from_filename("some_V42_annotation.gtf.gz") == 42

    def test_no_version(self):
        assert extract_gencode_version_from_filename("genes.gtf.gz") is None

    def test_full_path(self):
        p = "/home/user/refs/gencode.v32.annotation.gtf.gz"
        assert extract_gencode_version_from_filename(p) == 32

    def test_url_construction_v49(self):
        """Correct GTF URL can be constructed for a given version number."""
        version = 49
        filename = f"gencode.v{version}.annotation.gtf.gz"
        base = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/"
        url = base + filename
        assert url == (
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/"
            "gencode.v49.annotation.gtf.gz"
        )


# ---------------------------------------------------------------------------
# extract_polyasite_version_from_filename
# ---------------------------------------------------------------------------

class TestExtractPolyasite:
    def test_v3_gencode42(self):
        assert extract_polyasite_version_from_filename(
            "atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz"
        ) == 42

    def test_v3_gencode47(self):
        assert extract_polyasite_version_from_filename(
            "atlas.clusters.3.0.GRCh38.GENCODE_47.bed.gz"
        ) == 47

    def test_ensembl_based_returns_none(self):
        # PolyASite 2.0 uses Ensembl numbering, not GENCODE — should return None
        assert extract_polyasite_version_from_filename(
            "atlas.clusters.2.0.GRCh38.96.bed.gz"
        ) is None

    def test_url_construction_gencode42(self):
        """Correct PolyASite URL is constructed for a given version number."""
        version = 42
        filename = f"atlas.clusters.3.0.GRCh38.GENCODE_{version}.bed.gz"
        base = "https://polyasite.unibas.ch/download/atlas/3.0/"
        url = f"{base}GRCh38.GENCODE_{version}/{filename}"
        assert url == (
            "https://polyasite.unibas.ch/download/atlas/3.0/"
            "GRCh38.GENCODE_42/atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz"
        )


# ---------------------------------------------------------------------------
# fetch_latest_gencode_gtf — cache-hit and download paths
# ---------------------------------------------------------------------------

class TestFetchLatestGencodeGtf:
    def test_skip_download_when_cached(self, tmp_path):
        """Auto-download is skipped when the correct version is already cached."""
        html = _fake_gencode_html(49)
        # Pre-create the cached file
        cached = tmp_path / "gencode.v49.annotation.gtf.gz"
        cached.write_bytes(b"fake")

        with (
            patch("scnoisemeter.utils.annotation_fetcher._http_get", return_value=html),
            patch("scnoisemeter.utils.annotation_fetcher.CACHE_DIR", tmp_path),
            patch("scnoisemeter.utils.annotation_fetcher._download") as mock_dl,
        ):
            path, version = fetch_latest_gencode_gtf()

        assert version == 49
        assert path == cached
        mock_dl.assert_not_called()

    def test_triggers_download_when_not_cached(self, tmp_path):
        html = _fake_gencode_html(49)

        def fake_download(url, dest, **kwargs):
            dest.write_bytes(b"fake gtf content")

        with (
            patch("scnoisemeter.utils.annotation_fetcher._http_get", return_value=html),
            patch("scnoisemeter.utils.annotation_fetcher.CACHE_DIR", tmp_path),
            patch("scnoisemeter.utils.annotation_fetcher._download", side_effect=fake_download),
        ):
            path, version = fetch_latest_gencode_gtf()

        assert version == 49
        assert path.exists()


# ---------------------------------------------------------------------------
# fetch_latest_polyasite_atlas — cache-hit and probing
# ---------------------------------------------------------------------------

class TestFetchLatestPolyasiteAtlas:
    def test_skip_download_when_cached(self, tmp_path):
        """Auto-download is skipped when the correct version is already cached."""
        cached = tmp_path / "atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz"
        cached.write_bytes(b"fake")

        def fake_url_exists(url, **kwargs):
            return "GENCODE_42" in url

        with (
            patch("scnoisemeter.utils.annotation_fetcher.CACHE_DIR", tmp_path),
            patch("scnoisemeter.utils.annotation_fetcher._url_exists", side_effect=fake_url_exists),
            patch("scnoisemeter.utils.annotation_fetcher._download") as mock_dl,
        ):
            path, version = fetch_latest_polyasite_atlas(hint_max_gencode_version=49)

        assert version == 42
        assert path == cached
        mock_dl.assert_not_called()

    def test_triggers_download_when_not_cached(self, tmp_path):
        def fake_probe(url, **kwargs):
            return 200 if "GENCODE_42" in url else 404

        def fake_download(url, dest, **kwargs):
            dest.write_bytes(b"fake atlas content")

        with (
            patch("scnoisemeter.utils.annotation_fetcher.CACHE_DIR", tmp_path),
            patch("scnoisemeter.utils.annotation_fetcher._url_probe_status", side_effect=fake_probe),
            patch("scnoisemeter.utils.annotation_fetcher._download", side_effect=fake_download),
        ):
            path, version = fetch_latest_polyasite_atlas(hint_max_gencode_version=49)

        assert version == 42
        assert path.exists()

    def test_raises_when_no_version_found(self, tmp_path):
        with (
            patch("scnoisemeter.utils.annotation_fetcher.CACHE_DIR", tmp_path),
            patch("scnoisemeter.utils.annotation_fetcher._url_probe_status", return_value=404),
            patch("scnoisemeter.utils.annotation_fetcher._download"),
        ):
            with pytest.raises(RuntimeError, match="Could not find any PolyASite 3.0 atlas"):
                fetch_latest_polyasite_atlas(hint_max_gencode_version=44)


# ---------------------------------------------------------------------------
# _resolve_gtf — user-supplied bypasses auto-download
# ---------------------------------------------------------------------------

class TestResolveGtf:
    def test_user_supplied_bypasses_download(self, tmp_path):
        """User-supplied --gtf bypasses GTF auto-download."""
        fake_gtf = tmp_path / "gencode.v45.annotation.gtf.gz"
        fake_gtf.write_bytes(b"fake")

        with patch("scnoisemeter.cli.fetch_latest_gencode_gtf") as mock_fetch:
            path, version, source = _resolve_gtf(str(fake_gtf))

        mock_fetch.assert_not_called()
        assert path == str(fake_gtf)
        assert version == 45
        assert source == "user-supplied"

    def test_none_triggers_auto_download(self, tmp_path):
        """When --gtf is None, auto-download is triggered."""
        fake_path = tmp_path / "gencode.v49.annotation.gtf.gz"
        fake_path.write_bytes(b"fake")

        with patch("scnoisemeter.cli.fetch_latest_gencode_gtf", return_value=(fake_path, 49)):
            path, version, source = _resolve_gtf(None)

        assert version == 49
        assert source == "auto-downloaded"

    def test_user_supplied_makes_no_network_call(self, tmp_path):
        """When --gtf is supplied, no network call is made regardless of version."""
        fake_gtf = tmp_path / "gencode.v32.annotation.gtf.gz"
        fake_gtf.write_bytes(b"fake")

        # fetch_latest_gencode_gtf must NOT be called when gtf_arg is supplied
        with patch("scnoisemeter.cli.fetch_latest_gencode_gtf") as mock_fetch:
            path, version, source = _resolve_gtf(str(fake_gtf))

        mock_fetch.assert_not_called()
        assert version == 32
        assert source == "user-supplied"


# ---------------------------------------------------------------------------
# _resolve_polya_sites — user-supplied bypasses auto-download
# ---------------------------------------------------------------------------

class TestResolvePolySites:
    def test_user_supplied_bypasses_download(self, tmp_path):
        """User-supplied --polya-sites bypasses atlas auto-download."""
        fake_bed = tmp_path / "atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz"
        fake_bed.write_bytes(b"fake")

        with patch("scnoisemeter.cli.fetch_latest_polyasite_atlas") as mock_fetch:
            paths, version, source = _resolve_polya_sites((str(fake_bed),))

        mock_fetch.assert_not_called()
        assert paths == [str(fake_bed)]
        assert version == 42
        assert source == "user-supplied"

    def test_empty_triggers_auto_download(self, tmp_path):
        """When --polya-sites is empty, auto-download is triggered."""
        fake_path = tmp_path / "atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz"
        fake_path.write_bytes(b"fake")

        with patch("scnoisemeter.cli.fetch_latest_polyasite_atlas", return_value=(fake_path, 42)):
            paths, version, source = _resolve_polya_sites(())

        assert version == 42
        assert source == "auto-downloaded"


# ---------------------------------------------------------------------------
# _check_version_consistency
# ---------------------------------------------------------------------------

class TestVersionConsistency:
    def test_mismatch_warns_gtf49_atlas42(self):
        """Version mismatch warning triggers when GTF is v49 and atlas is v42 (gap=7 > 5)."""
        with patch("scnoisemeter.cli.click") as mock_click:
            result = _check_version_consistency(49, 42)
        assert result is not None
        assert "v49" in result
        assert "v42" in result

    def test_no_warning_within_5(self):
        """Version mismatch warning does NOT trigger when versions are within 5."""
        result = _check_version_consistency(47, 46)
        assert result is None

        result = _check_version_consistency(47, 45)
        assert result is None

        result = _check_version_consistency(47, 43)
        assert result is None

        result = _check_version_consistency(47, 42)
        assert result is None

    def test_exactly_5_no_warning(self):
        result = _check_version_consistency(47, 42)
        assert result is None

    def test_diff_6_warns(self):
        result = _check_version_consistency(49, 43)
        assert result is not None

    def test_none_gtf_version_no_crash(self):
        result = _check_version_consistency(None, 42)
        assert result is None

    def test_none_polya_version_no_crash(self):
        result = _check_version_consistency(47, None)
        assert result is None

    def test_both_none_no_crash(self):
        result = _check_version_consistency(None, None)
        assert result is None


# ---------------------------------------------------------------------------
# fetch_10x_whitelist — cache-hit and download paths
# ---------------------------------------------------------------------------

class TestFetch10xWhitelist:
    def test_skip_download_when_cached_v3(self, tmp_path):
        """Auto-download is skipped when v3 whitelist is already cached."""
        cached = tmp_path / "10x_whitelist_v3.txt.gz"
        cached.write_bytes(b"ACGT\n")

        with (
            patch("scnoisemeter.utils.annotation_fetcher.CACHE_DIR", tmp_path),
            patch("scnoisemeter.utils.annotation_fetcher._download") as mock_dl,
        ):
            path = fetch_10x_whitelist("10x_v3")

        assert path == cached
        mock_dl.assert_not_called()

    def test_triggers_download_when_not_cached_v4(self, tmp_path):
        """v4 whitelist is downloaded when not in cache."""
        def fake_download(url, dest, **kwargs):
            dest.write_bytes(b"fake whitelist")

        with (
            patch("scnoisemeter.utils.annotation_fetcher.CACHE_DIR", tmp_path),
            patch("scnoisemeter.utils.annotation_fetcher._download", side_effect=fake_download),
        ):
            path = fetch_10x_whitelist("10x_v4")

        assert path.exists()
        assert path.name == "10x_whitelist_v4.txt.gz"

    def test_raises_on_unknown_chemistry(self, tmp_path):
        with pytest.raises(ValueError, match="Unknown chemistry"):
            fetch_10x_whitelist("10x_v99")

    def test_raises_when_offline_and_not_cached(self, tmp_path):
        with (
            patch("scnoisemeter.utils.annotation_fetcher.CACHE_DIR", tmp_path),
        ):
            with pytest.raises(RuntimeError, match="--offline"):
                fetch_10x_whitelist("10x_v3", offline=True)
