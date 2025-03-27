"""Tests for the trgt module."""

import os

import pytest
from humanatee import trgt

# turn off docstring checking for this file
# ruff: noqa: D101 D102 D103


class TestAnnotateTrgt:
    cur_dir = os.path.dirname(__file__)
    vcf = os.path.join(cur_dir, 'test_data/trgt.vcf')
    pli_lookup = os.path.join(cur_dir, 'test_data/pli.lookup')
    oelof_lookup = os.path.join(cur_dir, 'test_data/oelof.lookup')
    clinvar_lookup = os.path.join(cur_dir, 'test_data/clinvar_gene_desc.txt')

    @pytest.fixture()
    def trgt_db(self, tmp_path):
        return trgt.TRGTdb(
            vcf=self.vcf,
            sample_id='SAMPLE1',
            gene_lookups=[self.pli_lookup, self.oelof_lookup, self.clinvar_lookup],
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
