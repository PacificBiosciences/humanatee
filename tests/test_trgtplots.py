"""Tests for the trgtpop module."""

import pytest
from humanatee import trgtplots

# turn off docstring checking for this file
# ruff: noqa: D101 D102 D103

# test_resources = resources.files('tests')
# popdb = test_resources.joinpath('test_data', 'hprc.trgt.vcf')


@pytest.mark.parametrize(
    'pop_al',
    [
        ('{"92,92": 709}'),
        (None),
    ],
)
def test_pophist_popal(tmp_path, pop_al):
    trgtplots.PopHist(
        title='Test title',
        pop_al=pop_al,
        short_allele_cutoff=92,
        long_allele_cutoff=92,
        filename=tmp_path / 'test.png',
    )


@pytest.mark.parametrize(
    'short_allele, long_allele',
    [
        (None, None),
        (42, 42),
    ],
)
def test_pophist_cutoffs(tmp_path, short_allele, long_allele):
    trgtplots.PopHist(
        title='Test title',
        pop_al='{"96,96": 405, "92,96": 112, "92,92": 10}',
        short_allele_cutoff=short_allele,
        long_allele_cutoff=long_allele,
        filename=tmp_path / 'test.png',
    )
