"""ML guard rails for the discovery program (#256).

Currently houses the **false-positive flagger** (a non-blocking probability
that a SILVER closure resembles one of our own past-bug signatures). The
flagger NEVER auto-rejects a candidate; it only routes high-probability
results to human re-check. See ``falsepos_flagger.FalsePosFlagger``.
"""

from cyclerfinder.ml.falsepos_features import (
    FEATURE_NAMES,
    SILVER_RECORD_KEYS,
    extract_features,
)
from cyclerfinder.ml.falsepos_flagger import FalsePosFlagger
from cyclerfinder.ml.falsepos_labels import (
    KNOWN_FALSE_POSITIVES,
    KNOWN_TRUE_REPRODUCTIONS,
    build_training_set,
)

__all__ = [
    "FEATURE_NAMES",
    "KNOWN_FALSE_POSITIVES",
    "KNOWN_TRUE_REPRODUCTIONS",
    "SILVER_RECORD_KEYS",
    "FalsePosFlagger",
    "build_training_set",
    "extract_features",
]
