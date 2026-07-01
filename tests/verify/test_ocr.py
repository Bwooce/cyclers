"""Tests for the Corpus Document Policy OCR helper (task #400).

Everything here SKIPS cleanly (never fails) when the optional ``ocr`` extra /
system binaries (tesseract, ghostscript, ocrmypdf, poppler) or the private
paper corpus is absent — mirroring how the SPICE/rebound
cross-check tests skip when their data/extra is missing. CI without the corpus
or the OCR toolchain therefore reports skips, not failures.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from cyclerfinder.verify.ocr import (
    ensure_text_layer,
    extract_text,
    has_text_layer,
)

# Private corpus root (see project memory: reference for the private corpus).
# Never linked from the public repo; this test only reads it if present locally.
_PAPERS = Path(os.environ.get("CYCLERS_CORPUS_DIR", "/home/bruce/dev/corpus/papers"))

# A modern, text-layer-abundant arXiv paper in the corpus (chosen for plentiful
# embedded text — the has_text_layer-true case).
_TEXT_LAYER_PDF = (
    _PAPERS / "braik-ross-2026-orbital-networks-three-body-problem-arxiv-2605.31543.pdf"
)

# The one confirmed image-only corpus book (policy §1 motivating example).
_SZEBEHELY = _PAPERS / "szebehely-1967-theory-of-orbits-restricted-problem-three-bodies-book.pdf"


def _have(*binaries: str) -> bool:
    return all(shutil.which(b) is not None for b in binaries)


def _first_text_layer_pdf() -> Path | None:
    """Return a corpus PDF that already has a usable text layer, if any.

    Prefers the named modern paper; otherwise scans papers/ for the first PDF
    that probes as text-layer, so the test does not hinge on one exact filename.
    """
    if not _have("pdftotext", "pdfinfo"):
        return None
    if _TEXT_LAYER_PDF.exists() and has_text_layer(_TEXT_LAYER_PDF):
        return _TEXT_LAYER_PDF
    if not _PAPERS.is_dir():
        return None
    for pdf in sorted(_PAPERS.glob("*.pdf")):
        try:
            if has_text_layer(pdf):
                return pdf
        except (subprocess.SubprocessError, OSError):
            continue
    return None


def test_has_text_layer_true() -> None:
    """A modern text-layer PDF probes as having a usable text layer (no OCR)."""
    if not _have("pdftotext", "pdfinfo"):
        pytest.skip("poppler (pdftotext/pdfinfo) not installed")
    pdf = _first_text_layer_pdf()
    if pdf is None:
        pytest.skip("no text-layer PDF available in the private corpus")
    assert has_text_layer(pdf) is True
    # ensure_text_layer must return the ORIGINAL path unchanged (no OCR step).
    assert ensure_text_layer(pdf) == pdf
    # And the convenience wrapper yields the real text.
    assert len(extract_text(pdf)) > 1000


@pytest.mark.slow
def test_ensure_text_layer_ocrs_image_only(tmp_path: Path) -> None:
    """OCR adds a text layer to an image-only scan (Szebehely 2-page sample).

    The full 661-page book is too slow for a unit test, so a 2-page sample is
    carved out with ghostscript and OCR'd; we assert post-OCR char count >>
    pre-OCR (which is ~0 for the image-only sample).
    """
    if not _have("pdftotext", "pdfinfo", "ocrmypdf", "tesseract", "gs"):
        pytest.skip("OCR toolchain (poppler/ocrmypdf/tesseract/ghostscript) not installed")
    if not _SZEBEHELY.exists():
        pytest.skip("Szebehely 1967 not present in the private corpus")

    # Carve a 2-page sample with ghostscript (pages 30-31: body text, not the
    # cover/blank front matter).
    sample = tmp_path / "szebehely_sample.pdf"
    subprocess.run(
        [
            "gs",
            "-q",
            "-dNOPAUSE",
            "-dBATCH",
            "-sDEVICE=pdfwrite",
            "-dFirstPage=30",
            "-dLastPage=31",
            f"-sOutputFile={sample}",
            str(_SZEBEHELY),
        ],
        check=True,
    )
    assert sample.exists()

    # Pre-OCR: image-only sample probes as NOT having a text layer.
    assert has_text_layer(sample) is False
    pre = subprocess.check_output(["pdftotext", "-layout", str(sample), "-"], text=True)
    pre_chars = len(pre.strip())

    # OCR into the tmp cache dir (NOT the default astropy cache — keep the test
    # self-contained and disposable).
    out = ensure_text_layer(sample, cache_dir=tmp_path / "ocr_cache")
    assert out != sample  # an OCR'd copy was produced
    assert out.exists()

    post = extract_text(sample, cache_dir=tmp_path / "ocr_cache")
    post_chars = len(post.strip())
    assert post_chars > pre_chars
    assert post_chars > 200  # a couple of real text pages yields plenty


def test_missing_pdf_raises(tmp_path: Path) -> None:
    """ensure_text_layer raises FileNotFoundError for a non-existent PDF."""
    if not _have("pdftotext", "pdfinfo"):
        pytest.skip("poppler (pdftotext/pdfinfo) not installed")
    with pytest.raises(FileNotFoundError):
        ensure_text_layer(tmp_path / "nope.pdf")
