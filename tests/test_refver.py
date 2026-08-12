#!/usr/bin/env python3
"""Unit tests for reference-version aliasing.

Last update: 2026-08-06

Covers:
- folding: every alias resolves to its canonical version
- uniqueness: no alias is claimed by two versions
- matching: case-insensitive and stripped, exact not substring
- unknown: folds to lowercase and is reported as unknown
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "config"))
from const import (  # noqa: E402
    REFVERS,
    REFVERS_ALIAS,
    canonical_refver,
    is_known_refver,
)


def test_every_alias_key_is_canonical():
    """REFVERS_ALIAS is keyed by the canonical versions, and covers all of them."""
    assert sorted(REFVERS_ALIAS) == sorted(REFVERS)


def test_canonical_versions_map_to_themselves():
    for refver in REFVERS:
        assert canonical_refver(refver) == refver


@pytest.mark.parametrize(
    "canon,alias", [(c, a) for c, al in REFVERS_ALIAS.items() for a in al]
)
def test_every_alias_resolves_to_its_canonical(canon, alias):
    assert canonical_refver(alias) == canon
    assert is_known_refver(alias)


def test_no_alias_is_claimed_twice():
    """An alias listed under two versions would make the fold order-dependent."""
    seen = [a.lower() for aliases in REFVERS_ALIAS.values() for a in aliases]
    assert len(seen) == len(set(seen))


@pytest.mark.parametrize("spelling", ["GRCh38", "grch38", "GRCH38", "  GRCh38  "])
def test_match_is_case_insensitive_and_stripped(spelling):
    assert canonical_refver(spelling) == "hg38"


def test_unknown_folds_to_lowercase_and_is_not_known():
    """An unsupported reference still yields a stable token, so it can select records."""
    assert canonical_refver("GRCh38-GIABv3") == "grch38-giabv3"
    assert not is_known_refver("GRCh38-GIABv3")


def test_matching_is_exact_not_substring():
    """A sub-flavor must not fold into its base build; that would mix two FASTAs."""
    assert canonical_refver("GRCh38-GIABv3") != "hg38"
    assert canonical_refver("hg38_alt") != "hg38"


def test_non_str_input_is_coerced():
    """A JSON number reaching the lookup must not raise."""
    assert canonical_refver(38) == "38"
