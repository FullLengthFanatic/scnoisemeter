"""
annotation_fetcher.py
=====================
Auto-downloads GENCODE GTF and PolyASite atlas to the user's cache directory
so that different labs always use a matched, versioned reference pair without
having to supply --gtf and --polya-sites manually.

GENCODE GTF
-----------
The current latest release is discovered by parsing the GENCODE FTP directory
listing at:
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/

Cache: ~/.cache/scnoisemeter/gencode.vXX.annotation.gtf.gz

Network calls are only made when the cache is empty.  If a cached file exists
it is returned immediately without any network access.

PolyASite atlas
---------------
PolyASite 3.0 (single-cell atlas) uses GENCODE version numbers in its path.
We probe candidate URLs from the current GENCODE version downward to find the
highest available atlas release.  Directory listing is disabled on the server,
so probing is the only option.

Cache: ~/.cache/scnoisemeter/atlas.clusters.3.0.GRCh38.GENCODE_XX.bed.gz

Network calls are only made when the cache is empty.
"""

from __future__ import annotations

import io
import logging
import re
import sys
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from typing import Optional, Tuple

from scnoisemeter import __version__

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CACHE_DIR = Path.home() / ".cache" / "scnoisemeter"

_GENCODE_LATEST_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/"
)
# URL template for a specific GENCODE release, e.g. release_42 → v42
_GENCODE_RELEASE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_{version}/gencode.v{version}.annotation.gtf.gz"
)
_POLYASITE3_BASE_URL = "https://polyasite.unibas.ch/download/atlas/3.0/"

# Earliest PolyASite 3.0 GRCh38 atlas we are aware of (GENCODE 42).
# Probing never goes below this.
_POLYASITE_MIN_GENCODE = 42

# PolyA_DB v4 — bulk+long-read validated PAS atlas (Wistar Institute, 2026).
# Contains two files: hg38.PAS.main.tsv (long-read validated, ~281 k sites)
# and hg38.PAS.max.tsv (all identified sites, ~1.43 M).
# We use main.tsv by default (higher confidence).
_POLYADB4_ZIP_URL = (
    "https://exon.apps.wistar.org/polya_db/v4/download/4.1/HumanPas.zip"
)
_POLYADB4_MAIN_ENTRY = "HumanPas/hg38.PAS.main.tsv"
_POLYADB4_CACHE_NAME = "polyadb4.hg38.PAS.main.bed"

# FANTOM5 — CAGE-based TSS atlas, phase 1+2, hg38 (RIKEN, 2017).
# BED6 format: chrom, start, end, name, score, strand (0-based).
# _load_tss_sites() uses columns 0-2 (midpoint), so BED6 works as-is.
_FANTOM5_CAGE_URL = (
    "https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/"
    "extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz"
)
_FANTOM5_CACHE_NAME = "fantom5.hg38.CAGE_peaks.bed.gz"

_USER_AGENT = f"scnoisemeter/{__version__} (https://github.com/scnoisemeter)"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fetch_latest_gencode_gtf(offline: bool = False) -> Tuple[Path, int]:
    """
    Return the path to a GENCODE human annotation GTF.

    Cache-first: if a cached ``gencode.vXX.annotation.gtf.gz`` file exists,
    it is returned immediately — **no network call is made**.  A network call
    is only made when the cache is empty (no matching file found).

    Parameters
    ----------
    offline : bool
        If True, never make network calls.  Raises RuntimeError if the cache
        is empty.

    Returns
    -------
    path : Path
        Path to the cached (possibly freshly downloaded) GTF.
    version : int
        GENCODE release number (e.g. 49).
    """
    cache_dir = _ensure_cache_dir()

    # Use cached file if available (no network call)
    cached = _find_cached_gencode_gtf(cache_dir)
    if cached is not None:
        path, version = cached
        logger.info("Using cached GENCODE v%d GTF: %s", version, path)
        print(f"Using GENCODE v{version} (from ~/.cache/scnoisemeter/)")
        return path, version

    # Cache is empty
    if offline:
        raise RuntimeError(
            "No cached GENCODE GTF found and --offline mode is active. "
            "Run once with network access to populate the cache, "
            "or supply --gtf explicitly."
        )

    # Make network call only now (cache was empty)
    html = _http_get(_GENCODE_LATEST_URL)
    filename = _parse_gencode_gtf_filename(html)
    version = _extract_gencode_version_int(filename)
    cached_path = cache_dir / filename

    url = _GENCODE_LATEST_URL + filename
    print(f"Downloading GENCODE v{version} annotation GTF …")
    _download(url, cached_path)
    print(f"Using GENCODE v{version} (downloaded to ~/.cache/scnoisemeter/)")

    return cached_path, version


def fetch_gencode_gtf_version(version: int, offline: bool = False) -> Tuple[Path, int]:
    """
    Return the path to a specific GENCODE human annotation GTF release.

    Cache-first: if ``gencode.v{version}.annotation.gtf.gz`` is already in
    the cache directory, it is returned immediately.

    Parameters
    ----------
    version : int
        GENCODE release number (e.g. 42).
    offline : bool
        If True, raise RuntimeError when the cache is empty.

    Returns
    -------
    path : Path
        Path to the cached (possibly freshly downloaded) GTF.
    version : int
        The requested GENCODE release number.
    """
    cache_dir = _ensure_cache_dir()
    filename = f"gencode.v{version}.annotation.gtf.gz"
    cached_path = cache_dir / filename

    if cached_path.exists() and cached_path.stat().st_size > 0:
        logger.info("Using cached GENCODE v%d GTF: %s", version, cached_path)
        print(f"Using GENCODE v{version} (from ~/.cache/scnoisemeter/)")
        return cached_path, version

    if offline:
        raise RuntimeError(
            f"No cached GENCODE v{version} GTF found and --offline mode is active. "
            f"Run once with network access to populate the cache, "
            f"or supply --gtf explicitly."
        )

    url = _GENCODE_RELEASE_URL.format(version=version)
    status = _url_probe_status(url)
    if status != 200:
        raise RuntimeError(
            f"GENCODE v{version} annotation GTF not found (HTTP {status}). "
            f"URL: {url}. "
            f"Check the release number or supply --gtf with a local file."
        )

    print(f"Downloading GENCODE v{version} annotation GTF …")
    _download(url, cached_path)
    print(f"Using GENCODE v{version} (downloaded to ~/.cache/scnoisemeter/)")
    return cached_path, version


def fetch_latest_polyasite_atlas(
    hint_max_gencode_version: Optional[int] = None,
    offline: bool = False,
) -> Tuple[Path, int]:
    """
    Find and return the path to the highest PolyASite 3.0 atlas for GRCh38.

    Cache-first: if a cached ``atlas.clusters.3.0.GRCh38.GENCODE_XX.bed.gz``
    file exists, it is returned immediately — **no network call is made**.
    A network call is only made when the cache is empty.

    Parameters
    ----------
    hint_max_gencode_version : int, optional
        Start probing from this GENCODE version (used only when a download is
        actually needed).  If None and the cache is empty, the GENCODE FTP is
        queried to discover the current latest version.
    offline : bool
        If True, never make network calls.  Raises RuntimeError if the cache
        is empty.

    Returns
    -------
    path : Path
        Path to the cached BED.gz atlas.
    version : int
        GENCODE version the atlas was built against (e.g. 42).
    """
    cache_dir = _ensure_cache_dir()

    # Use cached file if available (no network call)
    cached = _find_cached_polyasite_atlas(cache_dir)
    if cached is not None:
        path, version = cached
        logger.info(
            "Using cached PolyASite atlas GENCODE v%d: %s", version, path
        )
        print(
            f"Using PolyASite atlas v3.0 / GENCODE v{version} "
            f"(from ~/.cache/scnoisemeter/)"
        )
        return path, version

    # Cache is empty
    if offline:
        raise RuntimeError(
            "No cached PolyASite atlas found and --offline mode is active. "
            "Run once with network access to populate the cache, "
            "or supply --polya-sites explicitly."
        )

    # Make network call(s) only now (cache was empty)
    if hint_max_gencode_version is None:
        html = _http_get(_GENCODE_LATEST_URL)
        fname = _parse_gencode_gtf_filename(html)
        hint_max_gencode_version = _extract_gencode_version_int(fname)

    found_version: Optional[int] = None
    found_filename: Optional[str] = None
    all_network_errors: bool = True

    for v in range(hint_max_gencode_version, _POLYASITE_MIN_GENCODE - 1, -1):
        filename = f"atlas.clusters.3.0.GRCh38.GENCODE_{v}.bed.gz"
        url = f"{_POLYASITE3_BASE_URL}GRCh38.GENCODE_{v}/{filename}"
        status = _url_probe_status(url)
        if status is None:
            # Network error on this probe; keep trying but track it
            logger.debug("Network error probing PolyASite URL: %s", url)
        else:
            all_network_errors = False
            if status == 200:
                found_version = v
                found_filename = filename
                break

    if found_version is None or found_filename is None:
        if all_network_errors:
            raise RuntimeError(
                f"All PolyASite URL probes failed with network errors. "
                f"The server at {_POLYASITE3_BASE_URL} may be unreachable or "
                f"the URL scheme may have changed. "
                f"Supply --polya-sites explicitly to bypass auto-discovery."
            )
        raise RuntimeError(
            f"Could not find any PolyASite 3.0 atlas for GRCh38 at GENCODE "
            f"versions {_POLYASITE_MIN_GENCODE}–{hint_max_gencode_version}. "
            f"The server responded but no atlas file was found for these versions. "
            f"Check {_POLYASITE3_BASE_URL} or supply --polya-sites explicitly."
        )

    cached_path = cache_dir / found_filename
    url = (
        f"{_POLYASITE3_BASE_URL}GRCh38.GENCODE_{found_version}/{found_filename}"
    )
    print(f"Downloading PolyASite atlas GENCODE v{found_version} …")
    _download(url, cached_path)
    print(
        f"Using PolyASite atlas v3.0 / GENCODE v{found_version} "
        f"(downloaded to ~/.cache/scnoisemeter/)"
    )

    return cached_path, found_version


def fetch_polyadb4_atlas(offline: bool = False) -> Path:
    """
    Return the path to the PolyA_DB v4 GRCh38 PAS atlas (main collection).

    Cache-first: if ``polyadb4.hg38.PAS.main.bed`` exists in the cache
    directory, it is returned immediately — no network call is made.

    The original ZIP contains a tab-delimited TSV with a ``PAS_ID`` column
    in ``chr:strand:position`` format (1-based).  On first download, the
    TSV is converted to BED3 (0-based, chrom/start/end) and that BED is
    cached; subsequent runs read the BED directly.

    Parameters
    ----------
    offline : bool
        If True, raise RuntimeError when the cache is empty.

    Returns
    -------
    path : Path
        Path to the cached BED file.
    """
    cache_dir = _ensure_cache_dir()

    cached = _find_cached_polyadb4(cache_dir)
    if cached is not None:
        logger.info("Using cached PolyA_DB v4 atlas: %s", cached)
        print("Using PolyA_DB v4 PAS atlas (from ~/.cache/scnoisemeter/)")
        return cached

    if offline:
        raise RuntimeError(
            "No cached PolyA_DB v4 atlas found and --offline mode is active. "
            "Run once with network access to populate the cache, "
            "or supply --polya-sites explicitly."
        )

    dest = cache_dir / _POLYADB4_CACHE_NAME
    print("Downloading PolyA_DB v4 atlas (main collection) …")
    _download_polyadb4_to_bed(dest)
    print("Using PolyA_DB v4 PAS atlas (downloaded to ~/.cache/scnoisemeter/)")
    return dest


def fetch_fantom5_cage_peaks(offline: bool = False) -> Path:
    """
    Return the path to the FANTOM5 hg38 CAGE peak atlas (phase 1+2).

    Cache-first: if ``fantom5.hg38.CAGE_peaks.bed.gz`` exists in the cache
    directory, it is returned immediately — no network call is made.

    The file is BED6 (chrom, start, end, name, score, strand, 0-based).
    It is stored as-is; no format conversion is performed.

    Parameters
    ----------
    offline : bool
        If True, raise RuntimeError when the cache is empty.

    Returns
    -------
    path : Path
        Path to the cached BED.gz file.
    """
    cache_dir = _ensure_cache_dir()

    cached = _find_cached_fantom5(cache_dir)
    if cached is not None:
        logger.info("Using cached FANTOM5 CAGE atlas: %s", cached)
        print("Using FANTOM5 CAGE peak atlas (from ~/.cache/scnoisemeter/)")
        return cached

    if offline:
        raise RuntimeError(
            "No cached FANTOM5 CAGE atlas found and --offline mode is active. "
            "Run once with network access to populate the cache, "
            "or supply --tss-sites explicitly."
        )

    dest = cache_dir / _FANTOM5_CACHE_NAME
    print("Downloading FANTOM5 hg38 CAGE peak atlas …")
    _download(_FANTOM5_CAGE_URL, dest)
    print("Using FANTOM5 CAGE peak atlas (downloaded to ~/.cache/scnoisemeter/)")
    return dest


def extract_gencode_version_from_filename(filename: str) -> Optional[int]:
    """
    Extract a GENCODE release number from a filename, or return None.

    Handles patterns such as:
      gencode.v32.annotation.gtf.gz  →  32
      some_file_v38.gtf              →  38
      custom_V42_annotation.gtf.gz  →  42
    """
    match = re.search(r"[vV](\d+)", Path(filename).name)
    if match:
        return int(match.group(1))
    return None


def extract_polyasite_version_from_filename(filename: str) -> Optional[int]:
    """
    Extract the GENCODE version embedded in a PolyASite atlas filename.

    Handles patterns such as:
      atlas.clusters.3.0.GRCh38.GENCODE_42.bed.gz  →  42
      atlas.clusters.2.0.GRCh38.96.bed.gz            →  None  (Ensembl, not GENCODE)
    """
    match = re.search(r"GENCODE[_.](\d+)", Path(filename).name)
    if match:
        return int(match.group(1))
    return None


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _ensure_cache_dir() -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR


def _find_cached_gencode_gtf(cache_dir: Path) -> Optional[Tuple[Path, int]]:
    """
    Return the highest-version cached GENCODE GTF in *cache_dir*, or None.

    Looks for files matching ``gencode.v*.annotation.gtf.gz``.
    """
    best_version: Optional[int] = None
    best_path: Optional[Path] = None
    for f in cache_dir.glob("gencode.v*.annotation.gtf.gz"):
        v = extract_gencode_version_from_filename(f.name)
        if v is not None and (best_version is None or v > best_version):
            best_version = v
            best_path = f
    if best_path is not None:
        return best_path, best_version  # type: ignore[return-value]
    return None


def _find_cached_polyasite_atlas(cache_dir: Path) -> Optional[Tuple[Path, int]]:
    """
    Return the highest-version cached PolyASite 3.0 atlas in *cache_dir*, or None.

    Looks for files matching ``atlas.clusters.3.0.GRCh38.GENCODE_*.bed.gz``.
    """
    best_version: Optional[int] = None
    best_path: Optional[Path] = None
    for f in cache_dir.glob("atlas.clusters.3.0.GRCh38.GENCODE_*.bed.gz"):
        v = extract_polyasite_version_from_filename(f.name)
        if v is not None and (best_version is None or v > best_version):
            best_version = v
            best_path = f
    if best_path is not None:
        return best_path, best_version  # type: ignore[return-value]
    return None


def _find_cached_polyadb4(cache_dir: Path) -> Optional[Path]:
    """Return the cached PolyA_DB v4 BED if it exists and is non-empty."""
    p = cache_dir / _POLYADB4_CACHE_NAME
    if p.exists() and p.stat().st_size > 0:
        return p
    return None


def _find_cached_fantom5(cache_dir: Path) -> Optional[Path]:
    """Return the cached FANTOM5 CAGE BED.gz if it exists and is non-empty."""
    p = cache_dir / _FANTOM5_CACHE_NAME
    if p.exists() and p.stat().st_size > 0:
        return p
    return None


def _download_polyadb4_to_bed(dest: Path) -> None:
    """
    Download the PolyA_DB v4 ZIP, extract the main TSV, convert to BED3,
    and write atomically to *dest*.

    Conversion: PAS_ID column has format ``chr1:+:940182`` (1-based position).
    Each site becomes a BED3 row: ``chrom\\tstart\\tend`` where start = pos-1
    and end = pos (0-based half-open, one nucleotide wide).

    The ZIP is held in memory (≈55 MB); the output BED is streamed to a
    .tmp file and renamed on success.
    """
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    req = urllib.request.Request(_POLYADB4_ZIP_URL, headers={"User-Agent": _USER_AGENT})
    try:
        with urllib.request.urlopen(req) as resp:
            total = int(resp.getheader("Content-Length", 0) or 0)
            downloaded = 0
            buf = io.BytesIO()
            chunk_size = 1 << 20
            while True:
                chunk = resp.read(chunk_size)
                if not chunk:
                    break
                buf.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = downloaded / total * 100
                    sys.stderr.write(
                        f"\r  {downloaded / 1e6:.1f} / {total / 1e6:.1f} MB ({pct:.0f}%)"
                    )
                    sys.stderr.flush()
        if total:
            sys.stderr.write("\n")
            sys.stderr.flush()
        if total and downloaded != total:
            raise RuntimeError(
                f"Download truncated: expected {total} bytes, got {downloaded}. "
                f"URL: {_POLYADB4_ZIP_URL}"
            )
        buf.seek(0)
        n_written = 0
        with zipfile.ZipFile(buf) as zf:
            if _POLYADB4_MAIN_ENTRY not in zf.namelist():
                raise RuntimeError(
                    f"Expected {_POLYADB4_MAIN_ENTRY!r} inside the downloaded ZIP "
                    f"but it was not found.  The PolyA_DB v4 archive structure may "
                    f"have changed.  Supply --polya-sites explicitly as a workaround."
                )
            with zf.open(_POLYADB4_MAIN_ENTRY) as tsv_fh, open(tmp, "wb") as bed_fh:
                header_skipped = False
                for raw_line in tsv_fh:
                    line = raw_line.decode("utf-8", errors="replace").rstrip("\r\n")
                    if not header_skipped:
                        header_skipped = True
                        continue  # skip column header row
                    if not line:
                        continue
                    pas_id = line.split("\t", 1)[0]
                    parts = pas_id.split(":")
                    if len(parts) != 3:
                        continue
                    chrom, _strand, pos_str = parts
                    try:
                        pos_1based = int(pos_str)
                    except ValueError:
                        continue
                    start = pos_1based - 1  # convert to 0-based
                    end = pos_1based        # half-open: [start, end)
                    bed_fh.write(f"{chrom}\t{start}\t{end}\n".encode())
                    n_written += 1
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise
    tmp.rename(dest)
    logger.info("PolyA_DB v4: wrote %d BED3 sites to %s", n_written, dest)


def _parse_gencode_gtf_filename(html: str) -> str:
    """
    Find the primary (full) GENCODE annotation GTF filename in a directory
    listing HTML page.  The primary file is named:
      gencode.vXX.annotation.gtf.gz
    (not the basic/primary_assembly/chr_patch variants).
    """
    # Match the canonical annotation GTF — not basic, not primary_assembly, etc.
    # The href must end immediately after .gz (no extra suffix).
    for pattern in [
        r"(gencode\.v\d+\.annotation\.gtf\.gz)(?=[\"' <>])",
        r"(gencode\.v\d+\.annotation\.gtf\.gz)",
    ]:
        match = re.search(pattern, html)
        if match:
            return match.group(1)
    raise RuntimeError(
        "Could not find gencode.vXX.annotation.gtf.gz in GENCODE directory listing. "
        f"URL: {_GENCODE_LATEST_URL}"
    )


def _extract_gencode_version_int(filename: str) -> int:
    v = extract_gencode_version_from_filename(filename)
    if v is None:
        raise RuntimeError(
            f"Could not extract GENCODE version number from filename: {filename!r}"
        )
    return v


def _http_get(url: str, timeout: int = 30) -> str:
    """Fetch *url* and return the response body as a UTF-8 string."""
    req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def _url_exists(url: str, timeout: int = 15) -> bool:
    """Return True if *url* responds with HTTP 2xx, False otherwise."""
    return _url_probe_status(url, timeout=timeout) == 200


def _url_probe_status(url: str, timeout: int = 15) -> Optional[int]:
    """
    HEAD *url* and return the HTTP status code, or None on network error.

    Returning None (vs. a non-200 status code) lets callers distinguish
    "server reachable but file absent" from "server unreachable entirely".
    """
    try:
        req = urllib.request.Request(
            url, method="HEAD", headers={"User-Agent": _USER_AGENT}
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status
    except urllib.error.HTTPError as exc:
        return exc.code
    except (urllib.error.URLError, OSError):
        return None


def _download(url: str, dest: Path, chunk_size: int = 1 << 20) -> None:
    """Stream-download *url* to *dest* atomically.

    Writes to a .tmp sibling file and renames it to *dest* only on success.
    A partial download (interrupted connection, size mismatch) is deleted so
    the cache never contains a corrupt file.
    """
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    try:
        with urllib.request.urlopen(req) as resp:
            total = int(resp.getheader("Content-Length", 0) or 0)
            downloaded = 0
            with open(tmp, "wb") as fh:
                while True:
                    chunk = resp.read(chunk_size)
                    if not chunk:
                        break
                    fh.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        pct = downloaded / total * 100
                        sys.stderr.write(
                            f"\r  {downloaded / 1e6:.1f} / {total / 1e6:.1f} MB ({pct:.0f}%)"
                        )
                        sys.stderr.flush()
        if total:
            sys.stderr.write("\n")
            sys.stderr.flush()
        if total and downloaded != total:
            raise RuntimeError(
                f"Download truncated: expected {total} bytes, got {downloaded}. "
                f"URL: {url}"
            )
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise
    tmp.rename(dest)
