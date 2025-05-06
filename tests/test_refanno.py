"""Tests for the refanno module."""

import os
from importlib import resources

import pytest
from humanatee import refanno

# turn off docstring checking for this file
# ruff: noqa: D101 D102 D103
test_resources = resources.files('tests')
pop_vcf = test_resources.joinpath('test_data', 'trgt.hprc.hg.sorted.vcf')
pop_db = test_resources.joinpath('test_data', 'trgt.hprc.hg.db')
sorted_vcf = test_resources.joinpath('test_data', 'trgt.sorted.vcf')
gff = test_resources.joinpath('test_data', 'genes.gff3.gz')
bed = test_resources.joinpath('test_data', 'annotation.bed')


@pytest.mark.parametrize(
    'annotation',
    [
        ([gff, 'csq', 1000]),
        ([gff, 'csq', 1000, '--sort']),
        ([bed, 'csq', 1000]),
        ([bed, 'oddregions', 100, '--sort']),
        ([gff, 'csq', 1000, '--db', pop_db]),
    ],
)
def test_refanno_main(annotation, tmp_path):
    cmdargs = [
        '--vcf',
        pop_vcf,
        '--source',
        'trgt',
        '--prefix',
        tmp_path / 'test.trgtref',
        '--annotation',
    ] + annotation
    refanno.refanno_main([str(arg) for arg in cmdargs])
    # TODO: instead of checking for existence, check for content
    assert os.path.exists(tmp_path / 'test.trgtref.db')


@pytest.mark.parametrize(
    'bed_start, bed_end, strand, expected_result',
    [
        (5, 8, '+', -2),
        (5, 8, '-', 2),
        (5, 10, '+', 0),
        (5, 10, '-', 0),
        (5, 13, '+', 0),
        (5, 13, '-', 0),
        (5, 20, '+', 0),
        (5, 20, '-', 0),
        (5, 25, '+', 0),
        (5, 25, '-', 0),
        (10, 20, '+', 0),
        (10, 20, '-', 0),
        (18, 20, '+', 0),
        (18, 20, '-', 0),
        (20, 22, '+', 0),
        (20, 22, '-', 0),
        (27, 29, '+', 7),
        (27, 29, '-', -7),
    ],
)
def test_get_bp_distance(bed_start, bed_end, strand, expected_result):
    assert refanno.get_bp_distance(10, 20, bed_start, bed_end, strand) == expected_result
