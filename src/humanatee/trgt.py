"""Annotation of TRGT VCF."""

import argparse
import csv
import json
import logging
import os
import re
import sqlite3
import sys
from collections import Counter
from importlib import resources
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

import humanatee
from humanatee import dbutils, trgtplots, utils, vcf2db


class TRGTdb:
    """Annotate TRGT VCF based on user-provided files."""

    def __init__(self, sample_id, vcf, prefix, gene_lookups=[], phenotype=None):
        """Annotate VCF based on optional inputs."""
        self.sample_id = sample_id
        self.gene_lookup_labels = [os.path.basename(file.name).split('.')[0] for file in gene_lookups]
        if phenotype:
            self.gene_lookup_labels.insert(0, 'phenotype')
        self.prefix = prefix
        self.resources = resources.files('humanatee')

        logging.debug('Initializing output database, removing if it already exists')
        db_path = f'{prefix}.db'
        # TODO: add in-memory only option
        humanatee.silent_remove(db_path)
        if os.path.dirname(db_path) != '':
            os.makedirs(os.path.dirname(db_path), exist_ok=True)  # make parent directories if they don't exist
        self.conn = sqlite3.connect(db_path)
        self.cursor = self.conn.cursor()
        dbutils.smart_execute(self.cursor, 'PRAGMA foreign_keys = ON;')
        dbutils.smart_execute(self.cursor, 'PRAGMA synchronous=OFF;')

        logging.debug('Loading database schema for VCF')
        vcf_schema = self.resources.joinpath('schemas', 'vcf.sql').read_text()
        self.cursor.executescript(vcf_schema)

        if gene_lookups or phenotype:
            vcf2db.add_gene_lookups(self.cursor, gene_lookups, phenotype)

        logging.info('Reading sample TRGT VCF')
        v = pysam.VariantFile(vcf)
        variant_columns = [
            ('MOTIFS', 'TEXT', 'INFO', 'MOTIFS', None),
            ('STRUC', 'TEXT', 'INFO', 'STRUC', None),
            ('PS', 'INTEGER', 'FORMAT', 'PS', None),
        ]
        sample_allele_columns = [
            ('AL', 'INTEGER', 'FORMAT', 'AL', 'genotype_index'),
            ('ALLR', 'TEXT', 'FORMAT', 'ALLR', 'genotype_index'),
            ('SD', 'INTEGER', 'FORMAT', 'SD', 'genotype_index'),
            ('MC', 'TEXT', 'FORMAT', 'MC', 'genotype_index'),
            ('AM', 'REAL', 'FORMAT', 'AM', 'genotype_index'),
        ]
        vcf2db.VCFdb(
            self.cursor,
            v,
            self.sample_id,
            'trgt',
            table_columns={'Variant': variant_columns, 'SampleAllele': sample_allele_columns},
        )

    def commit_and_close(self):
        """Commit and close the database."""
        logging.info('Committing and closing database')
        self.conn.commit()
        self.cursor.close()
        self.conn.close()

    def write_tabular(self, filter):
        """Write TRGT data and filtered data to TSVs."""
        logging.info('Writing tabular output')
        results = dbutils.smart_execute(self.cursor, 'SELECT * FROM TRGT_main').fetchall()
        with open(f'{self.prefix}.tsv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow([_[0] for _ in self.cursor.description])  # header
            writer.writerows(results)
        if filter:
            df = pd.read_sql_query('SELECT * FROM TRGT_filtered', self.conn)
            df.to_csv(f'{self.prefix}.filtered.tsv', index=False, sep='\t')

    def write_xlsx(self):
        """Write TRGT filtered data to Excel."""
        logging.info('Writing Excel output')
        df = pd.read_sql_query('SELECT * FROM TRGT_filtered', self.conn)
        with pd.ExcelWriter(f'{self.prefix}.filtered.xlsx') as writer:
            df.to_excel(writer, sheet_name='repeats', index=False)

    def __load_pathogenic_info__(self, anno_tsv):
        """Load pathogenic annotations into Variant table."""
        logging.info('Loading pathogenic annotations')
        # add columns to Variant table for pathogenic info
        # TODO: require minimum columns but add extra columns if they exist
        columns = {
            'pathogenic_motifs': 'TEXT',
            'normal_min': 'INTEGER',
            'normal_max': 'INTEGER',
            'pathogenic_expansion_min': 'INTEGER',
            'pathogenic_contraction_max': 'INTEGER',
            'disease_gene': 'TEXT',
            'disease_id': 'TEXT',
            'disease_name': 'TEXT',
            'inheritance': 'TEXT',
            'notes': 'TEXT',
        }
        for label, type in columns.items():
            dbutils.add_column(self.cursor, 'Variant', label, type)

        # load pathogenic csv
        records = csv.DictReader(anno_tsv, delimiter='\t')
        for row in records:
            trid, motifs, struc = tuple([_.split('=')[1] for _ in row['info'].split(';')])
            start = int(row['start']) - 1  # bed to vcf
            items = {
                k: row[k] for k in row.keys() if row[k] and k in columns
            }  # dictreader returns empty strings for missing values
            items = ', '.join([f"{k}='{v}'" for k, v in items.items()])
            dbutils.smart_execute(
                self.cursor,
                f"""
                UPDATE Variant
                SET {items}
                WHERE variant_id = '{trid}'
                AND chrom = '{row["chrom"]}'
                AND start = '{start}'
                AND end = '{row["end"]}'
                AND motifs = '{motifs}'
                AND struc = '{struc}';
                """,
            )

    def __update_variant_pathogenicity__(
        self, sample_id, trid, idx, motifs, pathogenic_motifs, mc_string, normal_min, normal_max, expansion, contraction
    ):
        #  ('R01_89_SF_FAM861929_2997404D01_PRO1', 'FRAXE_AFF2', 0, 'GCC', 'GCC', '21', 4, 39, 200, None)
        pathogenic_motif_count = sum_mc(mc_string, motifs, pathogenic_motifs)
        pathogenicity = None
        if pathogenic_motif_count is not None:
            breaks = [contraction, normal_min, normal_max, expansion]
            ranges = [[0, contraction], [normal_min, normal_max], [expansion, float('inf')]]
            labels = ['Pathogenic', 'Normal', 'Pathogenic']
            if not contraction:
                breaks, ranges, labels = breaks[1:], ranges[1:], labels[1:]
                ranges[0][0] = 0
            if breaks != sorted(breaks) or (len(set(breaks)) != len(breaks) and normal_min != normal_max):
                raise ValueError(f'Normal or pathogenic thresholds provided for {trid} are not valid')
            pathogenicity = 'Premutation'
            for i, range in enumerate(ranges):
                if utils.in_range(pathogenic_motif_count, range):
                    pathogenicity = labels[i]
                    break

            x = f", pathogenic = '{pathogenicity}'" if pathogenicity else ''
            dbutils.smart_execute(
                self.cursor,
                f"""
                UPDATE SampleAllele
                SET pathogenic_MC = '{pathogenic_motif_count}'{x}
                WHERE variant_id = '{trid}'
                AND genotype_index = '{idx}'
                AND sample_id = '{sample_id}';
                """,
            )

    def identify_pathogenic_expansions(self, pathogenic_tsv):
        """Annotate pathogenic variants based on pathogenic annotations."""
        self.__load_pathogenic_info__(pathogenic_tsv)
        logging.info('Annotating pathogenic expansions')

        # add columns to SampleAllele table for pathogenic info
        columns = {'pathogenic_MC': 'INTEGER', 'pathogenic': 'TEXT'}
        for label, type in columns.items():
            dbutils.add_column(self.cursor, 'SampleAllele', label, type)

        # test pathogenicity for loci that are in Pathogenic table
        dbutils.smart_execute(
            self.cursor,
            """
            SELECT
                SampleAllele.sample_id,
                SampleAllele.variant_id,
                SampleAllele.genotype_index,
                Variant.motifs,
                Variant.pathogenic_motifs,
                SampleAllele.MC,
                Variant.normal_min,
                Variant.normal_max,
                Variant.pathogenic_expansion_min,
                Variant.pathogenic_contraction_max
            FROM Variant
            JOIN SampleAllele
            ON (Variant.variant_id = SampleAllele.variant_id)
            WHERE Variant.normal_min IS NOT NULL;
            """,
        )
        sample_alleles = self.cursor.fetchall()
        for x in sample_alleles:
            self.__update_variant_pathogenicity__(*x)

    def identify_population_expansions(self, pop_db):
        """Annotate expanded alleles based on population cutoffs."""
        logging.info('Annotating expansions based on population cutoffs')

        # add columns to Variant table for population info
        columns = {
            'pop_upper': 'INTEGER',
            'pop_cutoff_n': 'INTEGER',
            'pop_recessive_upper': 'INTEGER',
            'pop_recessive_cutoff_n': 'INTEGER',
            'pop_allele_length_genotypes': 'TEXT',
        }
        for label, type in columns.items():
            dbutils.add_column(self.cursor, 'Variant', label, type)

        # add columns to SampleAllele table for population info
        dbutils.add_column(self.cursor, 'SampleAllele', 'expanded', 'TEXT')

        # attach population database and update Variant table
        dbutils.smart_execute(self.cursor, f'ATTACH DATABASE "{pop_db}" AS pop')
        dbutils.smart_execute(
            self.cursor,
            """
                            UPDATE Variant
                            SET
                                pop_upper = p.upper,
                                pop_cutoff_n = p.n_alleles,
                                pop_recessive_upper = p.recessive_upper,
                                pop_recessive_cutoff_n = p.recessive_n_alleles,
                                pop_allele_length_genotypes = p.allele_length_genotypes
                            FROM pop.Locus AS p
                            WHERE Variant.variant_id = p.trid
                            AND Variant.chrom = p.chrom
                            AND Variant.start = p.start
                            AND Variant.end = p.end
                            AND Variant.motifs = p.motifs
                            AND Variant.struc = p.struc
                            """,
        )

        # update SampleAllele table with expanded status based on population cutoffs
        dbutils.smart_execute(
            self.cursor,
            """
            UPDATE SampleAllele as s
            SET expanded = CASE
                WHEN SUBSTR(s.ALLR, INSTR(s.ALLR,'-'),-20) > v.pop_upper THEN 'Expanded'
                WHEN SUBSTR(s.ALLR, INSTR(s.ALLR,'-'),-20) <= v.pop_recessive_upper THEN 'Normal'
                WHEN SUBSTR(s.ALLR, INSTR(s.ALLR,'-'),-20) <= v.pop_upper
                AND SUBSTR(s.ALLR, INSTR(s.ALLR,'-'),-20) > v.pop_recessive_upper THEN 'Expanded*'
                END
            FROM Variant as v
            WHERE s.variant_id = v.variant_id;
            """,
        )
        # only keep genotypes that are Expanded/Expanded, Expanded/Expanded*, Normal/Expanded, Expanded, Expanded*
        dbutils.smart_execute(
            self.cursor,
            """
            UPDATE SampleAllele
            SET expanded = NULL
            WHERE variant_id IN
            (SELECT variant_id
            FROM (SELECT variant_id, GROUP_CONCAT(expanded) AS expanded FROM SampleAllele GROUP BY variant_id)
            WHERE expanded = 'Expanded*,Normal'
            OR expanded = 'Normal,Expanded*'
            OR expanded = 'Expanded*,Expanded*'
            OR expanded = 'Normal,Normal'
            OR expanded = 'Normal')
            """,
        )

    def annotate_denovos(self, denovo, mother, father):
        """Annotate de novo variants based on denovo file."""
        logging.info('Annotating de novo variants')
        records = csv.DictReader(denovo, delimiter='\t')
        header = records.fieldnames

        if 'a_coverage' in header:
            if mother and father:
                raise ValueError('TRGT denovo was run in duo mode but both mother and father were provided')
            else:
                header = [re.sub(r'^a_', 'child_', h) for h in header]
                header = [re.sub(r'_a$', '_child', h) for h in header]
                x = 'mother' if mother else 'father'
                header = [re.sub(r'^b_', f'{x}_', h) for h in header]
                header = [re.sub(r'_b$', f'_{x}', h) for h in header]
                records.fieldnames = header
        else:
            if not mother and not father:
                raise ValueError('TRGT denovo was run in duo mode but either mother or father was not provided')

        # replace child_dropout with dropout
        index = header.index('child_dropout')
        records.fieldnames[index] = 'coverage_dropout'

        sample_allele_columns = {
            'denovo_coverage': 'INTEGER',
            'allele_coverage': 'INTEGER',
            'allele_ratio': 'REAL',
            'child_ratio': 'REAL',
            'allele_origin': 'TEXT',
            'denovo_status': 'TEXT',
        }
        sample_allele_columns = {k: v for k, v in sample_allele_columns.items() if k in header}
        variant_columns = {
            'father_dropout': 'TEXT',  # TODO: make Rust dropout tool that take multi-sample VCF and has columns for each sample
            'mother_dropout': 'TEXT',
            'mother_AL': 'TEXT',  # TODO: get parent AL/MC from merged VCF and then ignore these columns
            'father_AL': 'TEXT',
            'mother_MC': 'TEXT',
            'father_MC': 'TEXT',
        }
        variant_columns = {k: v for k, v in variant_columns.items() if k in header}
        for column, sql_type in sample_allele_columns.items():
            dbutils.add_column(self.cursor, 'SampleAllele', column, sql_type)
        for column, sql_type in variant_columns.items():
            dbutils.add_column(self.cursor, 'Variant', column, sql_type)

        trid = None
        for row in records:
            if row['denovo_status'] == '.':
                del row['denovo_status'], row['allele_ratio'], row['allele_coverage']
            sample_allele_items = ', '.join([f"{k}='{v}'" for k, v in row.items() if k in sample_allele_columns.keys()])
            dbutils.smart_execute(
                self.cursor,
                f"""
                UPDATE SampleAllele
                SET {sample_allele_items}
                WHERE
                    variant_id = '{row["trid"]}'
                    AND sample_id = '{self.sample_id}'
                    AND allele_index = '{row["genotype"]}'
                    AND genotype_index = '{row["index"]}';
                """,
            )
            variant_items = ', '.join(
                [f"{k}='{v}'" for k, v in row.items() if k in variant_columns.keys() and v != '.']
            )
            if row['trid'] != trid and len(variant_items) > 0:
                # update variant table
                dbutils.smart_execute(
                    self.cursor,
                    f"""
                    UPDATE Variant
                    SET {variant_items}
                    WHERE variant_id = '{row["trid"]}'
                    """,
                )
            trid = row['trid']

    def annotate_coverage_dropouts(self, dropouts):
        """Annotate coverage dropouts based on dropouts file."""
        logging.info('Annotating coverage dropouts')

        dbutils.add_column(self.cursor, 'Variant', 'coverage_dropout', 'TEXT')
        records = csv.DictReader(dropouts, delimiter='\t', fieldnames=['chrom', 'start', 'end', 'info', 'dropout'])
        for row in records:
            trid, motifs, struc = tuple([_.split('=')[1] for _ in row['info'].split(';')])
            start = int(row['start']) - 1  # bed to vcf
            dropout = 'FD' if row['dropout'] == 'FullDropout' else 'HD'
            dbutils.smart_execute(
                self.cursor,
                f"""
                UPDATE Variant
                SET coverage_dropout = '{dropout}'
                WHERE variant_id = '{trid}'
                AND chrom = '{row["chrom"]}'
                AND start = '{start}'
                AND end = '{row["end"]}'
                AND motifs = '{motifs}'
                AND struc = '{struc}';
                """,
            )

    def create_summary_views(self):
        """Create main TRGT views in DB."""
        vcf2db.create_summary_views(self.cursor, self.gene_lookup_labels, self.sample_id)

        # get sort order for columns in main TRGT views
        q = """SELECT v.*, s.*
            FROM variant_summary AS v
            LEFT JOIN sample_summary AS s
            ON v.variant_id = s.variant_id
            ORDER BY v.priority_rank
            """
        dbutils.smart_execute(
            self.cursor,
            f"""
            CREATE TEMP VIEW TRGT_temp AS
            {q}
            """,
        )
        columns = [_[1] for _ in dbutils.smart_execute(self.cursor, 'PRAGMA table_info(TRGT_temp);').fetchall()]
        column_order = []
        trgt_columns = self.resources.joinpath('data', 'trgt_columns.tsv').read_text()
        for line in trgt_columns.split('\n'):
            column_order.append(line.rstrip().split('\t')[0])
        if len(self.gene_lookup_labels) > 0:
            idx = column_order.index('gene') + 1
            column_order[idx:idx] = self.gene_lookup_labels

        columns = sorted([col for col in columns if col in column_order], key=column_order.index)

        # create main TRGT views
        dbutils.smart_execute(
            self.cursor,
            f"""
            CREATE VIEW TRGT_main AS
            SELECT {','.join(columns)} FROM
            ({q});
            """,
        )

        conditions = []
        if 'expanded' in columns:
            conditions += ["expanded LIKE '%Expanded%'"]
        if 'pathogenic' in columns:
            conditions += ["pathogenic LIKE '%Pathogenic%'", "pathogenic LIKE '%Premutation%'"]

        conditions = 'WHERE ' + ' OR '.join(conditions) if len(conditions) > 0 else ''
        dbutils.smart_execute(
            self.cursor,
            f"""
            CREATE VIEW TRGT_filtered AS
            SELECT * FROM TRGT_main {conditions};
            """,
        )
        if 'pathogenic' in columns:
            dbutils.smart_execute(
                self.cursor,
                f"""
                CREATE VIEW TRGT_pathogenic AS
                SELECT {','.join(columns)} FROM
                (SELECT v.*,s.*
                FROM variant_summary AS v
                LEFT JOIN sample_summary AS s
                ON v.variant_id = s.variant_id
                WHERE v.normal_min IS NOT NULL
                ORDER BY v.chrom, v.start);
                """,
            )

    def prioritize_and_plot(self, plot=True, priority_rank_max=2, allele_metric='AL'):
        """Add flags to variants based on pathogenicity, population cutoffs, denovo status, impact, etc."""
        logging.info(f"Prioritizing{' and plotting' if plot else ''} variants")

        directory = f'{self.prefix}_plots'
        plot_prefix = f"{directory}/{self.prefix.split('/')[-1]}"
        Path(directory).mkdir(parents=True, exist_ok=True)

        if 'phenotype' in self.gene_lookup_labels:
            sorted_phenotypes = dbutils.smart_execute(self.cursor, 'SELECT phenotype FROM Gene').fetchall()
            sorted_phenotypes = sorted([_[0] for _ in sorted_phenotypes if _[0] is not None and _[0] > 0], reverse=True)
            cutoff50 = np.percentile(sorted_phenotypes, 50)
        else:
            sorted_phenotypes = []
            cutoff50 = None

        dbutils.smart_execute(
            self.cursor,
            """
                SELECT
                    v.*, m.*
                FROM Variant as v
                JOIN
                TRGT_filtered as m
                ON v.variant_id = m.variant_id;
                """,
        )
        col_names = [_[0] for _ in self.cursor.description]
        # get plot settings based on if parents genotypes are provided
        plot_settings = [(allele_metric, 'Child', 'tab:olive')]
        if (mother := f'mother_{allele_metric}') in col_names:
            plot_settings.append((mother, 'Mom', 'tab:pink'))
        if (father := f'father_{allele_metric}') in col_names:
            plot_settings.append((father, 'Dad', 'tab:blue'))
        results = self.cursor.fetchall()

        for row in results:
            row_dict = {col: row[col_names.index(col)] for col in col_names}
            variant_id = row_dict['variant_id']
            parsed_row = FlagFilter(row_dict, sorted_phenotypes, cutoff50)

            dbutils.smart_execute(
                self.cursor,
                f"""
                UPDATE Variant
                SET priority_rank = '{parsed_row.rank}' {parsed_row.flags} {parsed_row.filters}
                WHERE variant_id = '{variant_id}'
                """,
            )

            if (
                plot
                and 'pop_upper' in col_names
                and parsed_row.rank <= priority_rank_max
                and 'population_monomorphic' not in parsed_row.flags
            ):
                # TODO: if pathogenic, maybe use MC instead of AL?
                trgtplots.PopHist(
                    title=f'TRID: {variant_id}',
                    pop_al=row_dict['pop_allele_length_genotypes'],
                    short_allele_cutoff=row_dict['pop_recessive_upper'],
                    long_allele_cutoff=row_dict['pop_upper'],
                    filename=f'{plot_prefix}.allele_length_hist.{variant_id.replace(",","_")[:45]}.png',
                    sample_genotypes=[row_dict[sample[0]] for sample in plot_settings],
                    sample_labels=[sample[1] for sample in plot_settings],
                    sample_colors=[sample[2] for sample in plot_settings],
                )

    def create_filtered_database(self, pathogenic=False):
        """Create a filtered database from the main database."""
        logging.info('Creating filtered database')
        # TODO: requires views TRGT_filtered and TRGT_pathogenic
        # TODO: add in-memory only option

        filtered_db = f'{self.prefix}.filtered.db'
        humanatee.silent_remove(filtered_db)
        dbutils.smart_execute(self.cursor, f'ATTACH DATABASE "{filtered_db}" AS filtered')
        dbutils.smart_execute(self.cursor, 'SELECT sql FROM sqlite_master')
        sql = self.cursor.fetchall()
        for line in sql:
            if line[0]:
                dbutils.smart_execute(
                    self.cursor,
                    line[0]
                    .replace('CREATE TABLE ', 'CREATE TABLE filtered.')
                    .replace('CREATE VIEW ', 'CREATE VIEW filtered.'),
                )

        for table in ['Gene', 'Variant', 'Consequence', 'SampleAllele']:
            logging.info(f'Populating filtered table: {table}')
            column_names = [_[1] for _ in dbutils.smart_execute(self.cursor, f'PRAGMA table_info({table});').fetchall()]
            pathogenic_ids = 'OR variant_id IN (SELECT variant_id FROM TRGT_pathogenic)' if pathogenic else ''
            where = (
                f'WHERE variant_id IN (SELECT variant_id FROM TRGT_filtered) {pathogenic_ids}'
                if table != 'Gene'
                else ''
            )
            dbutils.smart_execute(
                self.cursor,
                f"""
                INSERT INTO filtered.{table} ({', '.join(column_names)})
                SELECT * FROM {table}
                {where};
                """,
            )


def sum_mc(mc_string, motifs_string=None, target_motifs_string=None):
    """Sum motif counts from MC string."""
    mc_list = list(map(int, mc_string.split('_')))
    motifs = motifs_string.split(',') if motifs_string else None
    target_motifs = target_motifs_string.split(',') if target_motifs_string else None
    if motifs and target_motifs:
        if not all([motif in motifs for motif in target_motifs]):
            raise ValueError('Not all target motifs are in the motifs list.')
        else:
            indices = [motifs.index(motif) for motif in target_motifs if motif in motifs]
    else:
        indices = range(0, len(mc_list))
    return sum([mc_list[idx] for idx in indices])


class FlagFilter:
    """Filter and flag variants based on pathogenicity, population cutoffs, denovo status, impact, etc."""

    def __init__(self, row, sorted_phenotypes, phenotype_cutoff):
        self.flags = []
        self.filters = []

        self.filter_SD(row['SD'], row['AL'])
        if 'pathogenic' in row.keys():
            self.flag_pathogenic(row['pathogenic'])
        if 'expanded' in row.keys():
            self.flag_zygosity(row['expanded'])
            self.flag_allele_lengths(
                row['AL'],
                row['pop_allele_length_genotypes'],
                row['pop_upper'],
            )
        if 'denovo_status' in row.keys():
            self.flag_denovo(row['denovo_status'])
            self.filter_denovo(
                row['denovo_status'],
                row['allele_ratio'],
                row['child_ratio'],
            )
            if 'mother_dropout' in row.keys():
                self.filter_dropouts(row['mother_dropout'], 'mother_dropout')
            if 'father_dropout' in row.keys():
                self.filter_dropouts(row['father_dropout'], 'father_dropout')
        if 'coverage_dropout' in row.keys():
            self.filter_dropouts(row['coverage_dropout'], 'coverage_dropout')
        if 'highest_impact' in row.keys():
            self.flag_impact(row['highest_impact'])
        if 'phenotype' in row.keys():
            self.flag_phenotype(row['phenotype'], sorted_phenotypes, phenotype_cutoff)
        self.flags = ';'.join(self.flags)
        self.flags = f", flags = '{self.flags}'" if self.flags else ''
        self.filters = ';'.join(self.filters)
        self.filters = f", filters = '{self.filters}'" if self.filters else ''

        self.rank_priority(self.flags, self.filters)

    def flag_pathogenic(self, pathogenic_genotype):
        """Flag pathogenic variants."""
        if not pathogenic_genotype:
            return
        if 'Pathogenic' in pathogenic_genotype:
            self.flags.append('pathogenic')
        elif 'Premutation' in pathogenic_genotype:
            self.flags.append('premutation')

    def flag_zygosity(self, expanded_genotype):
        """Flag biallelic and hemizygous variants."""
        if not expanded_genotype:
            return
        if 'Normal' not in expanded_genotype:
            if ',' in expanded_genotype:
                self.flags.append('mode=biallelic')
            else:
                self.flags.append('mode=hemizygous')

    def flag_allele_lengths(self, sample_al_string, pop_allele_lengths, long_allele_cutoff):
        """Flag large expansions and population monomorphic variants."""
        if not all([sample_al_string, pop_allele_lengths, long_allele_cutoff]):
            return
        sample_alleles = [int(_) for _ in sample_al_string.split(',')]
        unique_alleles = set(
            humanatee.flatten(
                [genotype.split(',') for genotype in Counter(json.loads(pop_allele_lengths)).keys() if genotype != '']
            )
        )
        n_alleles = len(unique_alleles)
        if n_alleles <= 1:
            self.flags.append('misc=population_monomorphic')
        # is the expansion "big" (≥ 100% larger) compared to cutoff
        if (max(sample_alleles) - long_allele_cutoff) / (long_allele_cutoff + 1e-10) >= 1:
            self.flags.append('misc=large_expansion')

    def flag_denovo(self, denovo_status):
        """Flag de novo variants."""
        if not denovo_status:
            return
        if 'Y' in denovo_status:
            self.flags.append('mode=denovo')

    def flag_impact(self, impact):
        """Flag high and moderate impact variants."""
        if not impact:
            return
        if 'HIGH' in impact:
            self.flags.append('impact=high')
        elif 'MODERATE' in impact:
            self.flags.append('impact=moderate')

    def flag_phenotype(self, phenotype_string, sorted_phenotypes, cutoff50):
        """Flag variants with high phenotype scores."""
        if not phenotype_string or len(sorted_phenotypes) == 0:
            return
        phenotype_list = [float(_) for _ in phenotype_string.split(';') if _]
        if len(phenotype_list) == 0:
            return
        max_phenotype = max(phenotype_list)
        if max_phenotype >= cutoff50:
            self.flags.append(f'phenotype={sorted_phenotypes.index(max_phenotype) + 1}/{len(sorted_phenotypes)}')

    def filter_SD(self, SD, AL):
        """Filter variants based on spanning depth."""
        if not SD or not AL:
            return []
        SD = [int(_) for _ in SD.split(',')]
        AL = [int(_) for _ in AL.split(',')]
        # guidelines from Egor:
        # if AL <= 1000 bps, SD must be >=5
        # else SD must be >= 2 because some longer repeats are harder to sequence through so we expect to get fewer reads
        # SD ratio is useful for shorter repeats. For long expansions we expect very imbalanced coverage between short and expanded allele
        for sd, al in zip(SD, AL):
            if (al <= 1000 and sd < 5) or (al > 1000 and sd < 2):
                self.filters.append('SD')
                return

    def filter_dropouts(self, dropout, label):
        """Filter variants based on coverage dropouts."""
        if not dropout:
            return
        if 'FD' in dropout or 'HD' in dropout:
            self.filters.append(label)

    def filter_denovo(self, denovo_status, allele_ratio, child_ratio):
        """Filter de novo variants based on allele and child ratios."""
        if not all([denovo_status, allele_ratio, child_ratio]):
            return
        for denovo_status, allele_ratio, child_ratio in zip(
            denovo_status.split(','), allele_ratio.split(','), child_ratio.split(',')
        ):
            # guidelines from Tom:
            # allele_ratio should be >0.7: The closer this ratio comes to 1.0 the better.
            # Assuming that reads are partitioned correctly to their respective alleles,
            # then at a true de novo allele all reads should have de novo evidence (1.0 ratio).
            # Of course this ratio can be lower for many reasons: mosaicism, poor partitioning, read noise, etc...
            # child_ratio should be between 0.3 and 0.7: the range is chosen to balance the need to
            # detect de novo events while accounting for the possibility of allelic dropout or partial mosaicism.
            # Ratios below 0.3 or above 0.7 could indicate sample anomalies or errors in read assignment.
            if 'Y' in denovo_status:
                if float(allele_ratio) <= 0.7:
                    self.filters.append('allele_ratio')
                if float(child_ratio) < 0.3 or float(child_ratio) > 0.7:
                    self.filters.append('child_ratio')

    def rank_priority(self, flags, filters):
        """Rank variants based on flags and filters."""
        self.rank = 4
        if flags:
            if 'mode=' in flags:
                modes = 0
                if 'mode=denovo' in flags and not any(
                    [
                        filter in filters
                        for filter in ['allele_ratio', 'child_ratio', 'father_dropout', 'mother_dropout']
                    ]
                ):
                    modes = 1
                if 'mode=biallelic' in flags or 'mode=hemizygous' in flags:
                    modes = 1
                self.rank -= modes
            if 'impact=' in flags or 'misc=' in flags:
                self.rank -= 1
            if 'phenotype=' in flags:
                self.rank -= 1
            if 'pathogenic' in flags or 'premutation' in flags:
                self.rank = 1
            elif 'SD' in filters or 'coverage_dropout' in filters:
                self.rank = 4


def parse_args(cmdargs):
    """Pull the command line parameters."""
    parser = argparse.ArgumentParser(
        prog='humanatee trgt',
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # required named arguments
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument(
        '--vcf', metavar='VCF', required=True, type=argparse.FileType('r'), help='VCF file generated by TRGT'
    )
    requiredNamed.add_argument('--sample-id', metavar='ID', required=True, type=str, help='Sample ID in VCF')
    requiredNamed.add_argument('--prefix', metavar='PREFIX', required=True, type=str, help='Prefix for files')
    # optional named arguments
    parser.add_argument('--pathogenic-tsv', metavar='TSV', type=argparse.FileType('r'), help='Pathogenic annotations')
    parser.add_argument('--population-db', metavar='DB', type=str, help='Population cutoffs')
    parser.add_argument('--dropouts', type=argparse.FileType('r'), metavar='TSV', help='Coverage dropouts')
    parser.add_argument(
        '--gene-lookup',
        type=argparse.FileType('r'),
        default=[],
        metavar=('TSV'),
        action='append',
        help='A tab delimited, headerless file with gene lookup information',
    )
    parser.add_argument(
        '--phenotype',
        type=argparse.FileType('r'),
        metavar=('TSV'),
        help='A tab delimited, headerless file with 1:1 gene to phenotype information',
    )

    parser.add_argument('--denovo', type=argparse.FileType('r'), metavar='TSV', help='Output of trgt-denovo')
    parser.add_argument('--mother', action='store_true', default=False, help='Mother used for denovo')
    parser.add_argument('--father', action='store_true', default=False, help='Father used for denovo')
    parser.add_argument('--filter', action='store_true', default=False, help='Create filtered database')
    parser.add_argument('--tsv', action='store_true', default=False, help='Write tabular output')
    parser.add_argument('--xlsx', action='store_true', default=False, help='Write Excel output')
    parser.add_argument('--plot', action='store_true', default=False, help='Plot allele length histograms')
    parser.add_argument('--verbose', action='store_true', default=False, help='Verbose logging')

    parser._action_groups.reverse()
    args = parser.parse_args(cmdargs)

    if args.denovo and not any([args.mother, args.father]):
        parser.error('--denovo requires --mother and/or --father.')
    if args.plot and not args.population_db:
        parser.error('--plot requires --population-db.')
    return args


def trgt_main(cmdargs):
    """Run from command line."""
    args = parse_args(cmdargs)

    humanatee.setup_logging(args.verbose, sys.stderr, show_version=True)

    trgt_db = TRGTdb(
        vcf=args.vcf,
        sample_id=args.sample_id,
        gene_lookups=args.gene_lookup,
        phenotype=args.phenotype,
        prefix=args.prefix,
    )

    if args.denovo:
        trgt_db.annotate_denovos(args.denovo, args.mother, args.father)
    if args.dropouts:
        trgt_db.annotate_coverage_dropouts(args.dropouts)
    if args.pathogenic_tsv:
        trgt_db.identify_pathogenic_expansions(args.pathogenic_tsv)
    if args.population_db:
        trgt_db.identify_population_expansions(args.population_db)
    trgt_db.create_summary_views()
    trgt_db.prioritize_and_plot(args.plot)
    if args.filter or args.xlsx:
        trgt_db.create_filtered_database(args.pathogenic_tsv is not None)
    if args.tsv:
        trgt_db.write_tabular(args.filter)
    if args.xlsx:
        trgt_db.write_xlsx()
    trgt_db.commit_and_close()
