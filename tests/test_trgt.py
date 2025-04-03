"""Tests for the trgt module."""

import os
from importlib import resources

import pytest
from humanatee import trgt

# turn off docstring checking for this file
# ruff: noqa: D101 D102 D103

test_resources = resources.files('tests')
vcf = test_resources.joinpath('test_data', 'trgt.vcf')
pli_lookup = test_resources.joinpath('test_data', 'pli.lookup')
oelof_lookup = test_resources.joinpath('test_data', 'oelof.lookup')
clinvar_lookup = test_resources.joinpath('test_data', 'clinvar_gene_desc.txt')
phenotype = test_resources.joinpath('test_data', 'phrank.tsv')
popdb = test_resources.joinpath('test_data', 'trgt.hprc.hg.db')
denovo_duo = test_resources.joinpath('test_data', 'trgt.denovo.duo.tsv')
denovo_trio = test_resources.joinpath('test_data', 'trgt.denovo.trio.tsv')
pathogenic = test_resources.joinpath('test_data', 'trgt.pathogenic.tsv')
dropouts = test_resources.joinpath('test_data', 'trgt.dropouts.tsv')


class TestAnnotateTrgt:
    @pytest.fixture()
    def trgt_db(self, tmp_path):
        return trgt.TRGTdb(
            vcf=vcf,
            sample_id='SAMPLE1',
            gene_lookups=[pli_lookup, oelof_lookup, clinvar_lookup],
            prefix=tmp_path / 'SAMPLE1.trgt.humanatee',
        )

    def test_trgtdb_consequence(self, trgt_db):
        trgt_db.cursor.execute(
            """
            SELECT ranked_consequence
            FROM Consequence
            WHERE variant_id = 'chr1_146192773_146192858_A'
            AND gene = 'NBPF10'
            AND allele_index = 1
            """
        )
        results = trgt_db.cursor.fetchall()
        assert results == [('28_intron_variant',)]

    def test_trgtdb_sampleallele(self, trgt_db):
        trgt_db.cursor.execute(
            """
            SELECT MC
            FROM SampleAllele
            WHERE variant_id = 'chr1_146192773_146192858_A'
            AND sample_id = 'SAMPLE1'
            AND genotype_index = 0
            """
        )
        results = trgt_db.cursor.fetchall()
        assert results == [('27',)]

    def test_gene_lookups(self, trgt_db):
        trgt_db.cursor.execute(
            """
            SELECT oelof
            FROM Gene
            WHERE gene = 'ACTN2'
            """
        )
        results = trgt_db.cursor.fetchall()
        assert results == [('0.12109',)]

    def test_allr_substring(self, trgt_db):
        trgt_db.cursor.execute(
            """
            SELECT SUBSTR(ALLR, INSTR(ALLR,'-'),-20) FROM SampleAllele s
            WHERE s.variant_id = 'HMNR7_VWA1'
            AND genotype_index = 0
            """
        )
        results = trgt_db.cursor.fetchall()
        assert results == [('19',)]


@pytest.mark.parametrize(
    'optional_params, filename',
    [
        ([], 'SAMPLE1.trgt.humanatee.db'),
        (['--gene-lookup', pli_lookup], 'SAMPLE1.trgt.humanatee.db'),
        (['--gene-lookup', oelof_lookup, '--gene-lookup', clinvar_lookup], 'SAMPLE1.trgt.humanatee.db'),
        (['--phenotype', phenotype], 'SAMPLE1.trgt.humanatee.db'),
        (['--gene-lookup', pli_lookup, '--phenotype', phenotype], 'SAMPLE1.trgt.humanatee.db'),
        (['--denovo', denovo_duo, '--mother'], 'SAMPLE1.trgt.humanatee.db'),
        (['--denovo', denovo_trio, '--mother', '--father'], 'SAMPLE1.trgt.humanatee.db'),
        (['--pathogenic-tsv', pathogenic, '--filter', '--tsv'], 'SAMPLE1.trgt.humanatee.filtered.tsv'),
        (['--dropouts', dropouts], 'SAMPLE1.trgt.humanatee.db'),
        (['--xlsx'], 'SAMPLE1.trgt.humanatee.filtered.xlsx'),
        (['--tsv'], 'SAMPLE1.trgt.humanatee.tsv'),
        (['--population-db', popdb], 'SAMPLE1.trgt.humanatee.db'),
        (['--population-db', popdb, '--plot'], 'SAMPLE1.trgt.humanatee.db'),
        (['--population-db', popdb, '--denovo', denovo_duo, '--mother'], 'SAMPLE1.trgt.humanatee.db'),
        (['--population-db', popdb, '--denovo', denovo_trio, '--mother', '--father'], 'SAMPLE1.trgt.humanatee.db'),
        (
            ['--population-db', popdb, '--plot', '--pathogenic-tsv', pathogenic],
            'SAMPLE1.trgt.humanatee_plots/SAMPLE1.trgt.humanatee.allele_length_hist.FECD3_TCF4.png',
        ),
    ],
)
def test_trgt_main_success(tmp_path, optional_params, filename):
    required_params = ['--vcf', vcf, '--sample-id', 'SAMPLE1', '--prefix', tmp_path / 'SAMPLE1.trgt.humanatee']
    trgt.trgt_main([str(_) for _ in required_params + optional_params])
    assert os.path.exists(tmp_path / filename)
