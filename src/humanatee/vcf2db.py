"""Module for converting VCF files to a SQLite database."""

import logging
import os

import pandas as pd

from humanatee import dbutils, utils


class VCFdb:
    """Class for managing VCF to database operations."""

    def __init__(
        self, cursor, vcf, sample_id, source, table_columns, sub_source='', sequence=False, keep_filters=[], csq=True
    ):
        self.cursor = cursor
        self.vcf = vcf
        self.sample_id = sample_id
        self.source = source
        self.sub_source = sub_source
        self.record_source = (
            f'{source}_{sub_source}' if sub_source else source
        )  # avoids duplicate variant_ids from 2+ paraphase VCFs
        self.table_columns = utils.flatten(
            [[self.TableColumn(table, *col) for col in cols] for table, cols in table_columns.items()]
        )
        self.csq_format = []
        self.first_record = None

        filename = os.path.join(os.path.dirname(__file__), 'data/ensembl_consequences.tsv')
        self.consequence_order = list(pd.read_csv(filename, sep='\t', header=0)['so_term'])
        self.impact_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']

        if source not in ['trgt', 'sawfish', 'deepvariant', 'paraphase']:
            raise ValueError(f'VCF type {source} not recognized')

        self.samples = list(vcf.header.samples)
        if sample_id not in self.samples and source != 'paraphase':
            raise ValueError(f'Sample {sample_id} not found in VCF')
        if source == 'paraphase':
            logging.info(f'Haplotypes: {self.samples}')
        else:
            logging.info(f'Samples: {self.samples}')

        self.add_columns_to_tables()
        if 'CSQ' in self.vcf.header.info.keys() and csq:
            self.csq_format = vcf.header.info['CSQ'].description.split('Format: ')[1].split('|')
            self.prep_csq()

        # iterate through variant records
        for record in self.vcf:
            if not self.first_record:
                self.first_record = record
            if keep_filters and (record.filter.keys() == {} or record.filter.keys()[0] not in keep_filters):
                continue

            self.add_record(record)

    class TableColumn:
        """Class for adding custom columns to tables based on VCF fields."""

        def __init__(self, table, sql_label, sql_type, vcf_column, field_name, index):
            self.table = table
            self.sql_label = sql_label
            self.sql_type = sql_type
            self.vcf_column = vcf_column
            self.field_name = field_name
            self.index = index

    def get_variant_id(self, record):
        """Get variant ID based on source."""
        match self.source:
            case 'trgt':
                variant_id = record.info['TRID']
            case 'paraphase':
                variant_id = f'{record.chrom}-{record.pos}-{record.ref}-{record.alts[0]}'
            case _:
                variant_id = record.id
        if isinstance(variant_id, tuple):
            variant_id = ','.join(variant_id)
        return variant_id

    def add_columns_to_tables(self):
        """Add custom columns to tables."""
        for col in self.table_columns:
            dbutils.add_columns(self.cursor, col.table, {col.sql_label: col.sql_type})

        # add haplotype_id if source is paraphase
        if self.source == 'paraphase':
            dbutils.add_columns(self.cursor, 'SampleAllele', {'haplotype_id': 'TEXT'})

    def prep_csq(self):
        """Prepare CSQ format for parsing."""
        expected_csq_fields = ['ALLELE_NUM', 'SYMBOL', 'Gene', 'BIOTYPE', 'Consequence', 'IMPACT']
        for col in self.table_columns:
            if col.vcf_column == 'INFO' and col.field_name == 'CSQ':
                expected_csq_fields.append(col.index)
                # replace index label string with index integer
                col.index = self.csq_format.index(col.index) if col.index in self.csq_format else None
        if not all([_ in self.csq_format for _ in expected_csq_fields]):
            raise ValueError(f'Not all required fields found in CSQ: {expected_csq_fields}')

    def add_record(self, record):
        """Add a record to the database."""
        variant_id = self.get_variant_id(record)
        logging.debug(f'Variant ID: {variant_id}')
        self.add_vcf_item(
            self.cursor,
            'Variant',
            record,
            self.sample_id,
            [_ for _ in self.table_columns if _.table == 'Variant'],
            {
                'source': self.record_source,
                'variant_id': variant_id,
                'chrom': record.chrom,
                'start': record.pos,
                'end': record.stop,
            },
        )
        if self.csq_format:
            self.add_consequence(
                self.cursor,
                variant_id,
                self.record_source,
                record,
                self.consequence_order,
                self.impact_order,
                self.csq_format,
            )
        for sample in self.samples:
            for genotype_index, allele_index in enumerate(record.samples[sample]['GT']):
                if self.source == 'paraphase':
                    fields = {
                        'variant_id': variant_id,
                        'genotype_index': self.samples.index(sample),
                        'allele_index': allele_index,
                        'sample_id': self.sample_id,
                        'source': self.record_source,
                        'haplotype_id': sample,
                    }
                else:
                    fields = {
                        'variant_id': variant_id,
                        'genotype_index': genotype_index,
                        'allele_index': allele_index,
                        'sample_id': sample,
                        'source': self.record_source,
                    }

                self.add_vcf_item(
                    self.cursor,
                    'SampleAllele',
                    record,
                    sample,
                    [_ for _ in self.table_columns if _.table == 'SampleAllele'],
                    fields,
                )

    def add_consequence(self, cursor, variant_id, source, record, consequence_order, impact_order, csq_format):
        """Add consequence information to Consequence table."""
        if 'CSQ' not in record.info.keys():
            return
        csq_list = [csq_string.split('|') for csq_string in record.info['CSQ']]
        # relevant_alleles = record.samples[sample_id]["GT"]
        # csq_list = [_ for _ in csq_list if int(_[0]) in relevant_alleles]
        if not csq_list:
            return
        for csq in csq_list:
            # ALLELE_NUM,SYMBOL,Gene,BIOTYPE,Consequence,IMPACT,HGVSc,HGVSp,CLIN_SIG,gnomADg_AF
            csq_dict = dict(zip(csq_format, csq))
            gene_symbol = None if csq_dict['SYMBOL'] == '' else csq_dict['SYMBOL']
            gene = gene_symbol if gene_symbol else csq_dict['Gene']
            csq_out = '/'.join(
                [gene] + [csq_dict[_] for _ in ['BIOTYPE', 'Consequence', 'IMPACT', 'HGVSc', 'HGVSp', 'HGVS_OFFSET']]
            )
            consequence_rank = min([consequence_order.index(_) for _ in csq_dict['Consequence'].split('&')]) + 1
            ranked_consequence = f'{str(consequence_rank).zfill(2)}_{csq_dict["Consequence"]}'
            impact_rank = impact_order.index(csq_dict['IMPACT']) + 1
            ranked_impact = f'{str(impact_rank)}_{csq_dict["IMPACT"]}'
            values = {
                'variant_id': variant_id,
                'source': source,
                'allele_index': csq_dict['ALLELE_NUM'],
                'gene': gene_symbol,
                'ensembl_id': csq_dict['Gene'],
                'csq': csq_out,
                'ranked_consequence': ranked_consequence,
                'ranked_impact': ranked_impact,
            }
            try:
                dbutils.add_table_row(cursor, 'Consequence', values)
            except Exception as e:
                if 'FOREIGN KEY constraint failed' in str(e):
                    dbutils.add_table_row(cursor, 'Gene', {'gene': gene_symbol})
                    dbutils.add_table_row(cursor, 'Consequence', values)
                    return
                else:
                    raise e

    def add_vcf_item(self, cursor, table, record, record_sample, record_fields=[], fields={}):
        """Add a row to a SQL table with at least some fields defined in VFC record."""
        record_dict = {'INFO': record.info, 'FORMAT': record.samples[record_sample]}
        column_dict = {}
        if ('allele_index' not in fields.keys() or fields['allele_index'] is None) and table == 'SampleAllele':
            return
        for col in record_fields:
            # get value from record
            if col.vcf_column in record_dict.keys():
                try:
                    value = record_dict[col.vcf_column][col.field_name]
                except KeyError:
                    value = None
            else:
                column_dict[col.sql_label] = col.vcf_column
                continue

            if col.field_name == 'CSQ':
                value = value[0].split('|') if value else None
            index = col.index
            if isinstance(index, str):
                index = fields[index]  # for genotype_index or allele_index
            if isinstance(value, list) or isinstance(value, tuple):
                if index is not None:
                    try:
                        value = value[index]
                    except IndexError:
                        value = None
                else:
                    value = ','.join(map(str, value))
            column_dict[col.sql_label] = value
        column_dict.update(fields)
        dbutils.add_table_row(cursor, table, column_dict)


def add_gene_lookups(cursor, gene_lookups, phenotype=None):
    """Add gene lookup tables to database. Best to perform before reading VCF to avoid unnecessary updates."""
    logging.info('Adding gene lookups to db')
    df = pd.DataFrame(columns=['gene'])
    for gene_lookup in gene_lookups:
        label = os.path.basename(gene_lookup.name).split('.')[0].replace('-', '_').replace('.', '_')
        new_df = pd.read_csv(gene_lookup, sep='\t', names=['gene', label], dtype=str)
        new_df = new_df.astype(str).groupby('gene').agg(lambda x: '|'.join(x)).reset_index()
        df = df.merge(new_df, on='gene', how='outer')
    if phenotype:
        df = df.merge(pd.read_csv(phenotype, sep='\t', names=['gene', 'phenotype']), on='gene', how='outer')
    placeholder = tuple(df.to_dict(orient='records'))
    for label in df.columns[1:]:
        sql_type = 'REAL' if label == 'phenotype' else 'TEXT'
        dbutils.add_columns(cursor, 'Gene', {label: sql_type})
    values = f':{", :".join(df.columns)}'
    dbutils.smart_execute(
        cursor,
        f"""
        INSERT INTO Gene
        VALUES({values})
        """,
        (placeholder),
        many=True,
    )


def create_consequence_summary_view(cursor):
    """Create a view that summarizes consequence information for each variant."""
    dbutils.smart_execute(
        cursor,
        """
        CREATE VIEW consequence_summary AS
        SELECT
            variant_id,
            source,
            gene,
            MIN(ranked_consequence) AS highest_consequence,
            MIN(ranked_impact) AS highest_impact,
            csq
        FROM Consequence
        GROUP BY variant_id, ensembl_id
        """,
    )


def create_gene_summary_view(cursor, gene_lookup_labels):
    """Create a view that summarizes gene information for each variant. Requires existence of consequence summary view."""
    group_lookups = [f', GROUP_CONCAT(g.{label}, ";") AS {label}' for label in gene_lookup_labels]
    coalesce_lookups = [f', COALESCE({label}, "") AS {label}' for label in gene_lookup_labels]
    dbutils.smart_execute(
        cursor,
        f"""
        CREATE VIEW gene_summary AS
        SELECT
            c.variant_id,
            c.source,
            MIN(highest_consequence) AS highest_consequence,
            MIN(highest_impact) AS highest_impact,
            GROUP_CONCAT(c.gene, ';') AS gene,
            GROUP_CONCAT(c.csq, ';') AS consequence {' '.join(group_lookups)}
        FROM consequence_summary AS c
        LEFT JOIN (SELECT gene{''.join(coalesce_lookups)} FROM Gene) AS g
        ON c.gene = g.gene
        GROUP BY c.variant_id, c.source
        """,
    )


def create_variant_summary_view(cursor):
    """Create a view that combines variant and gene information. Requires existence of gene summary view."""
    # get names of columns without duplicating variant_id
    gene_columns = [
        f'g.{_[1]}'
        for _ in dbutils.smart_execute(cursor, 'PRAGMA table_info(gene_summary);').fetchall()
        if _[1] != 'variant_id'
    ]
    gene_columns = ', '.join(gene_columns)
    dbutils.smart_execute(
        cursor,
        f"""
        CREATE VIEW variant_summary AS
        SELECT
            v.*, {gene_columns}
        FROM Variant AS v
        LEFT JOIN gene_summary AS g
        ON v.variant_id = g.variant_id
        """,
    )


def create_sample_summary_view(cursor, sample_id, no_concat_columns=[]):
    """Create a view that summarizes sample information for each variant."""
    columns = dbutils.smart_execute(cursor, 'PRAGMA table_info(SampleAllele);').fetchall()
    columns = [_[1] for _ in columns]
    concat_columns = [
        _ for _ in columns if _ not in ['variant_id', 'source', 'sample_id'] and _ not in no_concat_columns
    ]

    coalesce = [f', COALESCE({col}, "") AS {col}' for col in concat_columns]
    group_concat = [f', NULLIF(GROUP_CONCAT({col}, ","), ",") AS {col}' for col in concat_columns]
    no_concat = [f', {col}' for col in no_concat_columns if col in columns]

    dbutils.smart_execute(
        cursor,
        f"""
        CREATE VIEW sample_summary AS
        SELECT
            variant_id, source, sample_id {' '.join(group_concat)} {' '.join(no_concat)}
        FROM (SELECT variant_id, source, sample_id {' '.join(coalesce)} {' '.join(no_concat)} FROM SampleAllele)
        WHERE sample_id = '{sample_id}'
        GROUP BY variant_id, sample_id
        """,
    )


def create_summary_views(cursor, gene_lookup_labels, sample_id, no_concat_columns=[]):
    """Create all summary views."""
    create_consequence_summary_view(cursor)
    create_sample_summary_view(cursor, sample_id, no_concat_columns)
    create_gene_summary_view(cursor, gene_lookup_labels)
    create_variant_summary_view(cursor)
