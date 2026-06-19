"""One-call OCR / text-layer helper for the Corpus Document Policy (task #400).

This is the *code layer* of the OCR-first rule in
``docs/notes/corpus-document-policy.md`` §1. A document in the private
``cyclers_pdf/papers/`` corpus is only "processed" once it is text-searchable;
this module makes that step a single idempotent call any future digest agent
reuses, instead of vision-reading hundreds of page-images through the Read tool
(the cost trap the policy was adopted to close — Szebehely 1967, a 661-page
image-only scan, burned large token cost that way).

Policy §1 flow this module implements:

1. **Probe for an existing text layer** (``pdftotext`` + page count). Most
   modern arXiv/journal PDFs already carry one — return the original path
   unchanged, no OCR.
2. **OCR image-only PDFs once** with ``ocrmypdf --skip-text`` (Tesseract under
   the hood — a cheap deterministic CPU step), caching the result.
3. Return the path to the text-layer PDF; callers then ``pdftotext`` it.

Hybrid caveat (policy §1, "Hybrid for precision")
-------------------------------------------------
Tesseract garbles math / subscripts / tables on old scans. The OCR'd layer
produced here is for **navigation** — finding the right chapter/section/page.
Precision sourced values (equations, table cells) must still come from Claude
vision (the Read tool on the page image) on the 2-3 precision-critical pages,
or from the original typeset PDF. Do not quote a sourced numeric constant from
this OCR text.

Caching convention (mirrors ``verify/spice_kernels.py``)
--------------------------------------------------------
The OCR'd PDF is written into a ``cyclerfinder_ocr`` subdirectory of the
astropy cache dir — the same on-demand-produce-if-missing / cache / return-path
pattern the SPICE leapseconds kernel uses. The generated PDF is large, binary,
and regenerable, so like the cached BSP/LSK it is **never committed to the
repo**.

The ``ocr`` extra (``ocrmypdf`` + the ``tesseract-ocr`` / ``ghostscript``
system binaries) is optional. This module imports nothing from it at import
time and degrades with a clear message if it is absent, so the core never
depends on it.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

# Probe threshold (documented constant; see policy §1). At or above this many
# pdftotext characters per page ⇒ usable text layer (no OCR); below ⇒ image-only.
#
# Measured on the corpus (2026-06-19): a genuine text layer yields ~3700
# chars/page (braik-ross-2026 arXiv: 215409 chars / 58 pages); an image-only
# scan yields EXACTLY ~1 char/page — pdftotext emits one form-feed (\f) byte per
# page and nothing else (Szebehely 1967: 342 chars / 342 pages = 1.000). The two
# populations are ~3 orders of magnitude apart, so the threshold sits anywhere in
# the gap. A floor of 1.0 is WRONG: it coincides exactly with the image-only
# value and misclassifies image-only scans as text-bearing (the bug the slow
# OCR test caught). 10 chars/page is ~10x above the form-feed noise floor and
# ~370x below real text — comfortably inside the gap.
#
# Edge case (honest): a scan carrying a *thin pre-existing* OCR layer (a few
# words per page) could read just above this floor and wrongly skip a re-OCR.
# That is rare; for such a document call with a higher min_chars_per_page or
# force OCR explicitly. We do not over-engineer per-page distribution analysis
# here — the bimodal corpus does not need it.
DEFAULT_MIN_CHARS_PER_PAGE = 10.0

# ocrmypdf flags (policy §1): --skip-text leaves any page that already has text
# untouched (so a partly-OCR'd PDF is not double-burned) and OCRs only the
# image pages; --quiet keeps the subprocess output clean for programmatic use.
_OCRMYPDF_FLAGS = ("--skip-text", "--quiet")

_OCR_EXTRA_HINT = (
    "the 'ocr' extra is required for this step (ocrmypdf + the tesseract-ocr "
    "and ghostscript system binaries). Install with `uv sync --extra ocr` and "
    "ensure `tesseract` and `gs` are on PATH (apt: tesseract-ocr ghostscript)."
)


def _require_binary(name: str) -> str:
    """Return the path to system binary ``name`` or raise a clear ImportError.

    The OCR step shells out to external tools rather than importing a Python
    package, so the "optional dependency missing" failure is surfaced the same
    ImportError-style way an absent ``import`` would be — pointing at the ``ocr``
    extra (policy §1 tooling note).
    """
    path = shutil.which(name)
    if path is None:
        raise ImportError(f"{name!r} not found on PATH: {_OCR_EXTRA_HINT}")
    return path


def _page_count(pdf_path: Path) -> int:
    """Return the PDF page count.

    Prefers ``pdfinfo`` (poppler, the same toolchain as ``pdftotext``); falls
    back to ``pypdfium2`` (already an ocrmypdf dependency) if pdfinfo is absent.
    """
    pdfinfo = shutil.which("pdfinfo")
    if pdfinfo is not None:
        out = subprocess.check_output([pdfinfo, str(pdf_path)], text=True)
        for line in out.splitlines():
            if line.startswith("Pages:"):
                return int(line.split(":", 1)[1].strip())
        raise RuntimeError(f"pdfinfo gave no Pages: line for {pdf_path}")
    try:
        import pypdfium2 as pdfium  # type: ignore[import-untyped]
    except ImportError as exc:  # pragma: no cover - defensive
        raise ImportError(
            f"neither pdfinfo nor pypdfium2 available to count pages: {_OCR_EXTRA_HINT}"
        ) from exc
    doc = pdfium.PdfDocument(str(pdf_path))
    try:
        return len(doc)
    finally:
        doc.close()


def _pdftotext(pdf_path: Path) -> str:
    """Extract the PDF's current text layer via ``pdftotext`` (poppler).

    ``-layout`` preserves the physical layout (policy §1 uses ``pdftotext
    -layout``); ``-`` streams to stdout so nothing is written to disk.
    """
    pdftotext = _require_binary("pdftotext")
    return subprocess.check_output([pdftotext, "-layout", str(pdf_path), "-"], text=True)


def has_text_layer(
    pdf_path: str | os.PathLike[str],
    min_chars_per_page: float = DEFAULT_MIN_CHARS_PER_PAGE,
) -> bool:
    """Return True if the PDF already carries a usable text layer.

    Probes with ``pdftotext`` and divides the extracted character count by the
    page count (policy §1's "text-layer probe"). At or above
    ``min_chars_per_page`` ⇒ usable text layer (no OCR needed); below ⇒
    image-only.
    """
    pdf = Path(pdf_path)
    n_chars = len(_pdftotext(pdf))
    n_pages = max(_page_count(pdf), 1)
    return (n_chars / n_pages) >= min_chars_per_page


def ensure_text_layer(
    pdf_path: str | os.PathLike[str],
    *,
    cache_dir: str | os.PathLike[str] | None = None,
    min_chars_per_page: float = DEFAULT_MIN_CHARS_PER_PAGE,
) -> Path:
    """Return a path to a text-searchable version of ``pdf_path`` (policy §1).

    If the PDF already has a usable text layer (``has_text_layer``), the
    original path is returned unchanged. Otherwise the PDF is OCR'd once with
    ``ocrmypdf --skip-text --quiet`` into ``cache_dir`` and the cached path is
    returned.

    Caching mirrors ``verify/spice_kernels.ensure_leapseconds_kernel``:
    ``cache_dir`` defaults to a ``cyclerfinder_ocr`` subdirectory of the astropy
    cache dir, and the OCR step is idempotent — a cached OCR'd file that is
    newer than the source is reused (no re-OCR). The OCR'd PDF is large/binary/
    regenerable and is never committed to the repo.

    Hybrid caveat: the returned PDF's OCR text is for navigation; precision
    sourced values still come from vision on the page image or the original (see
    the module docstring).
    """
    pdf = Path(pdf_path)
    if not pdf.exists():
        raise FileNotFoundError(f"PDF not found: {pdf}")

    if has_text_layer(pdf, min_chars_per_page=min_chars_per_page):
        return pdf

    if cache_dir is None:
        from astropy.config.paths import get_cache_dir

        cache_dir = Path(get_cache_dir()) / "cyclerfinder_ocr"
    cache_path = Path(cache_dir)
    cache_path.mkdir(parents=True, exist_ok=True)
    out_path = cache_path / f"{pdf.stem}.ocr.pdf"

    # Idempotent: reuse a cached OCR'd file that is newer than the source.
    if out_path.exists() and out_path.stat().st_mtime >= pdf.stat().st_mtime:
        return out_path

    ocrmypdf = _require_binary("ocrmypdf")
    _require_binary("tesseract")  # surfaced early with the ocr-extra hint
    subprocess.run(
        [ocrmypdf, *_OCRMYPDF_FLAGS, str(pdf), str(out_path)],
        check=True,
    )
    return out_path


def extract_text(
    pdf_path: str | os.PathLike[str],
    *,
    cache_dir: str | os.PathLike[str] | None = None,
    min_chars_per_page: float = DEFAULT_MIN_CHARS_PER_PAGE,
) -> str:
    """Return the text of ``pdf_path``, OCR-ing first if it is image-only.

    Convenience wrapper: ``ensure_text_layer`` then ``pdftotext -layout``. This
    is the one call a digest agent makes to get navigable text from any corpus
    document regardless of whether it arrived with a text layer (policy §1).
    Precision values still come from vision on the page (module docstring).
    """
    text_pdf = ensure_text_layer(
        pdf_path, cache_dir=cache_dir, min_chars_per_page=min_chars_per_page
    )
    return _pdftotext(text_pdf)
