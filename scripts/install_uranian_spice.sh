#!/usr/bin/env bash
# Install JPL NAIF SPICE kernels for the Uranian system into the local GMAT
# install (Phase 4.1 / #335 Part A).
#
# What this script installs
# -------------------------
# * ura111.bsp -- the legacy "URA111" Uranian satellite ephemeris.
#   Sources: ARIEL (701), UMBRIEL (702), TITANIA (703), OBERON (704),
#   MIRANDA (705), URANUS (799), plus EARTH-BARYCENTER (3) and SUN (10)
#   used as parent-frame anchors.
#   Time span: 1900-01-01 to 2099-12-24 ET.
#   Size: 162 MB (the smallest URA kernel covering all five regular moons).
#   URL: https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/ura111.bsp
#
# Why ura111 (not the larger ura116xl / ura184 / ura111xl-* releases)
# -------------------------------------------------------------------
# * The XL ("extended-length") releases are extreme-precision multi-GB
#   files (e.g. ura116xl.bsp ~ 692 MB; ura111xl-701.bsp ~ 2 GB) needed for
#   spacecraft navigation against the Uranus 1986 Voyager-2 / 2020s
#   stellar-occultation campaigns; not warranted for cycler-class V4-strict
#   gauntlets where the moon eccentricity / inclination corrections
#   dominate the perturbation budget.
# * The 184 series covers only the planetary barycenters and Earth; it has
#   no satellite content.
# * ura111 is the standard JPL/NAIF entry-level Uranian satellite kernel:
#   the one that GMAT R2022a is documented against.
#
# GMAT already ships
# ------------------
# * Leap-seconds: ~/GMAT/R2022a/data/time/SPICELeapSecondKernel.tls
# * Planetary constants: ~/GMAT/R2022a/data/planetary_coeff/SPICEPlanetaryConstantsKernel.tpc
# * DE405 / DE421 / DE424 planetary ephemerides:
#     ~/GMAT/R2022a/data/planetary_ephem/spk/DE405AllPlanets.bsp etc.
#
# Verification
# ------------
# After install, run:
#   uv run python scripts/verify_uranian_spice.py
# which loads the kernel via spiceypy, queries Umbriel + Oberon states at
# J2000+0d, J2000+1d, J2000+7d, and prints the SMA / eccentricity / inclination
# the kernel implies. Expected eccentricities (Murray-Dermott Table A.7):
#   e_Miranda  = 0.0013
#   e_Ariel    = 0.0012
#   e_Umbriel  = 0.0039
#   e_Titania  = 0.0011
#   e_Oberon   = 0.0014
#
# Run as
# ------
#   bash scripts/install_uranian_spice.sh
#
# The download is idempotent (curl -C - resumes); re-running this script
# after a partial or completed download is safe.

set -euo pipefail

GMAT_ROOT="${GMAT_ROOT:-$HOME/GMAT/R2022a}"
URA_DIR="${GMAT_ROOT}/data/planetary_ephem/spk/uranian"
URA_KERNEL_URL="https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/a_old_versions/ura111.bsp"
URA_KERNEL_NAME="ura111.bsp"
URA_KERNEL_EXPECTED_BYTES="169928704"

if [[ ! -d "${GMAT_ROOT}" ]]; then
    echo "error: GMAT install not found at ${GMAT_ROOT}" >&2
    echo "Set GMAT_ROOT to the install path and re-run." >&2
    exit 1
fi

mkdir -p "${URA_DIR}"

target="${URA_DIR}/${URA_KERNEL_NAME}"
if [[ -f "${target}" ]]; then
    have_bytes=$(stat -c %s "${target}")
    if [[ "${have_bytes}" == "${URA_KERNEL_EXPECTED_BYTES}" ]]; then
        echo "[ura111] already installed at ${target} (${have_bytes} bytes); skipping download."
    else
        echo "[ura111] partial file at ${target} (${have_bytes} of ${URA_KERNEL_EXPECTED_BYTES} bytes); resuming."
        curl -C - -L --fail -o "${target}" "${URA_KERNEL_URL}"
    fi
else
    echo "[ura111] downloading ${URA_KERNEL_URL} -> ${target}"
    curl -L --fail -o "${target}" "${URA_KERNEL_URL}"
fi

# Quick verification: ensure size matches.
actual_bytes=$(stat -c %s "${target}")
if [[ "${actual_bytes}" != "${URA_KERNEL_EXPECTED_BYTES}" ]]; then
    echo "warning: ura111.bsp size ${actual_bytes} != expected ${URA_KERNEL_EXPECTED_BYTES}" >&2
    echo "The download may be incomplete or the kernel may have been updated upstream." >&2
fi

# Sanity-print the first 4 bytes (NAIF DAF header magic).
header=$(head -c 8 "${target}")
if [[ "${header}" =~ ^DAF/SPK ]]; then
    echo "[ura111] DAF/SPK header OK"
else
    echo "warning: unexpected header (got '${header}'); file may be corrupted" >&2
fi

echo
echo "URA SPICE install complete:"
echo "  kernel: ${target}"
echo "  size:   ${actual_bytes} bytes"
echo
echo "Next step: uv run python scripts/verify_uranian_spice.py"
