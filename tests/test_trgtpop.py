"""Tests for the trgtpop module."""

import os
from importlib import resources

from humanatee import trgtpop

# turn off docstring checking for this file
# ruff: noqa: D101 D102 D103

test_resources = resources.files('tests')
pop_vcf = test_resources.joinpath('test_data', 'trgt.hprc.vcf')
pop_db1 = test_resources.joinpath('test_data', 'trgt.hprc.hg.db')
pop_db2 = test_resources.joinpath('test_data', 'trgt.hprc.na.db')


def test_trgtpop_main(tmp_path):
    cmdargs = ['--vcf', str(pop_vcf), '--cohort-id', 'hprc', '--prefix', tmp_path / 'test.trgtpop']
    trgtpop.trgtpop_main(str(_) for _ in cmdargs)
    assert os.path.exists(tmp_path / 'test.trgtpop.db')


def tests_trgtpop_main_merge(tmp_path):
    cmdargs = ['--dbs', pop_db1, pop_db2, '--prefix', tmp_path / 'test.trgtpop']
    trgtpop.trgtpop_main(str(_) for _ in cmdargs)
    assert os.path.exists(tmp_path / 'test.trgtpop.db')


def tests_trgtpop_main_merge_dupdbs(tmp_path):
    cmdargs = ['--dbs', pop_db2, pop_db2, '--prefix', tmp_path / 'test.trgtpop']
    trgtpop.trgtpop_main(str(_) for _ in cmdargs)
    assert os.path.exists(tmp_path / 'test.trgtpop.db')
