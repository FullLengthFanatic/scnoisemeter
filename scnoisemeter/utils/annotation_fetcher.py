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

import logging
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional, Tuple

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CACHE_DIR = Path.home() / ".cache" / "scnoisemeter"

_GENCODE_LATEST_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/"
)
_POLYASITE3_BASE_URL = "https://polyasite.unibas.ch/download/atlas/3.0/"

# Earliest PolyASite 3.0 GRCh38 atlas we are aware of (GENCODE 42).
# Probing never goes below this.
_POLYASITE_MIN_GENCODE = 42

_USER_AGENT = "scnoisemeter/0.1.6 (https://github.com/scnoisemeter)"


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

    for v in range(hint_max_gencode_version, _POLYASITE_MIN_GENCODE - 1, -1):
        filename = f"atlas.clusters.3.0.GRCh38.GENCODE_{v}.bed.gz"
        url = f"{_POLYASITE3_BASE_URL}GRCh38.GENCODE_{v}/{filename}"
        if _url_exists(url):
            found_version = v
            found_filename = filename
            break

    if found_version is None or found_filename is None:
        raise RuntimeError(
            f"Could not find any PolyASite 3.0 atlas for GRCh38 at GENCODE "
            f"versions {_POLYASITE_MIN_GENCODE}–{hint_max_gencode_version}. "
            f"Check {_POLYASITE3_BASE_URL}"
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
    try:
        req = urllib.request.Request(
            url, method="HEAD", headers={"User-Agent": _USER_AGENT}
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return 200 <= resp.status < 300
    except (urllib.error.HTTPError, urllib.error.URLError, OSError):
        return False


def _download(url: str, dest: Path, chunk_size: int = 1 << 20) -> None:
    """Stream-download *url* to *dest*, printing a simple progress indicator."""
    req = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    with urllib.request.urlopen(req) as resp:
        total = int(resp.getheader("Content-Length", 0) or 0)
        downloaded = 0
        with open(dest, "wb") as fh:
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
