"""Get population summary statistics from multi-sample TRGT VCF to use for annotation."""

import argparse
import json
import logging
import sqlite3
import sys
from collections import Counter

import numpy as np
import pysam

from humanatee import dbutils, utils


class LocusStats:
    """Calculate summary statistics for a locus."""

    def __init__(self, al_genotypes, sd_genotypes=None):
        # sort alleles in each genotype by length with shortest allele first
        # if sd is provided (i.e. not merging databases), filter out genotypes where at least one allele doesn't meet SD criteria
        if sd_genotypes is None:
            sorted_genotypes = [sorted(_) for _ in al_genotypes]
        else:
            sorted_genotypes = [
                self.__filter_sd__(al_tuple, sd_tuple) for al_tuple, sd_tuple in zip(al_genotypes, sd_genotypes)
            ]
            sorted_genotypes = [sorted(_) for _ in sorted_genotypes if _ is not None]

        all_alleles = utils.flatten(sorted_genotypes)
        self.upper = self.__get_cutoffs__(all_alleles)
        self.cutoff_n = len(all_alleles)

        short_alleles = [next(iter(sample), '') for sample in sorted_genotypes]
        short_alleles = [_ for _ in short_alleles if _ != '']
        self.recessive_upper = self.__get_cutoffs__(short_alleles)
        self.recessive_cutoff_n = len(short_alleles)

        self.al_genotype_dict = Counter(
            [','.join([str(allele) for allele in genotype]) for genotype in sorted_genotypes]
        )

    def __filter_sd__(self, al_tuple, sd_tuple):
        """Filter genotypes where at least one allele doesn't meet SD criteria."""
        for al, sd in zip(al_tuple, sd_tuple):
            if not all([al, sd]):
                return None
            if al <= 1000 and sd < 5:
                return None
            elif al > 1000 and sd < 2:
                return None
        return list(al_tuple)

    def __get_cutoffs__(self, counts):
        """Get outlier cutoff based on quantiles."""
        if len(counts) > 0:
            upper = np.percentile(counts, 99)
            upper = int(np.ceil(upper))
        else:
            upper = 0
        return upper


def vcf_to_db(cursor, vcf):
    """Convert multi-sample VCF to DB with aggregate statistics."""
    v = pysam.VariantFile(vcf)
    samples = list(v.header.samples)
    logging.info(f'Samples: {samples}')

    # add records to table
    for record in v:
        trid = record.info['TRID']
        if isinstance(trid, tuple):
            trid = ','.join(record.info['TRID'])
        logging.debug(f'TRID: {trid}')
        al_genotypes = [record.samples[sample]['AL'] for sample in samples]
        sd_genotypes = [record.samples[sample]['SD'] for sample in samples]
        locus_stats = LocusStats(al_genotypes, sd_genotypes)
        row_dict = {
            'variant_id': trid,
            'source': 'trgt',
            'chrom': record.chrom,
            'start': record.pos,  # TRGT BED.start and VCF.pos are the same
            'end': record.stop,
            'motifs': ','.join(record.info['MOTIFS']),
            'struc': record.info['STRUC'],
            'pop_upper': locus_stats.upper,
            'pop_n_alleles': locus_stats.cutoff_n,
            'pop_recessive_upper': locus_stats.recessive_upper,
            'pop_recessive_n_alleles': locus_stats.recessive_cutoff_n,
            'pop_allele_length_genotypes': json.dumps(locus_stats.al_genotype_dict),
        }
        dbutils.add_table_row(cursor, 'Variant', row_dict, replace=False)

    v.close()


def merge_dbs(conn, cursor, db_list):
    """Merge list of popDBs into a single DB."""
    logging.info(f'Merging databases: {db_list}')

    def merge_dbs_helper(cursor, db2, required_tables):
        """Merge second database into cursor."""
        db2_conn = sqlite3.connect(db2)
        db2_cursor = db2_conn.cursor()
        dbutils.validate_db_tables(db2_cursor, required_tables)

        # update Locus table
        db2_cursor.execute(f'SELECT {",".join(required_tables["Variant"])} FROM Variant')
        for row in db2_cursor.fetchall():
            row = {k: v for k, v in zip(required_tables['Variant'], row)}
            try:
                dbutils.add_table_row(cursor, 'Variant', row, replace=False)
            except sqlite3.IntegrityError:
                # handle primary key conflict by replace row
                cursor.execute(
                    'SELECT pop_allele_length_genotypes FROM Variant where variant_id = ?', (row['variant_id'],)
                )
                original_al_genotypes = cursor.fetchone()[0]
                new_al_genotypes = utils.reverse_counter(row['pop_allele_length_genotypes']) + utils.reverse_counter(
                    original_al_genotypes
                )
                locus_stats = LocusStats(new_al_genotypes)
                row_dict = {
                    'variant_id': row['variant_id'],
                    'source': 'trgt',
                    'chrom': row['chrom'],
                    'start': row['start'],
                    'end': row['end'],
                    'motifs': row['motifs'],
                    'struc': row['struc'],
                    'pop_upper': locus_stats.upper,
                    'pop_n_alleles': locus_stats.cutoff_n,
                    'pop_recessive_upper': locus_stats.recessive_upper,
                    'pop_recessive_n_alleles': locus_stats.recessive_cutoff_n,
                    'pop_allele_length_genotypes': json.dumps(locus_stats.al_genotype_dict),
                }
                dbutils.add_table_row(cursor, 'Variant', row_dict, replace=True)
        db2_conn.close()

    # what table and columns do we need to merge
    columns = cursor.execute('PRAGMA table_info(Variant);').fetchall()
    required_tables = {'Variant': [column[1] for column in columns]}
    # merge databases one by one
    for db in db_list:
        merge_dbs_helper(cursor, db, required_tables)
        conn.commit()


def trgtpop_main(cmdargs):
    """Run from command line."""
    parser = argparse.ArgumentParser(description=__doc__, prog='humanatee trgtpop')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('--prefix', metavar='PREFIX', required=True, type=str, help='Prefix for output')

    mainInput = parser.add_mutually_exclusive_group(required=True)
    mainInput.add_argument('--vcf', metavar='VCF', type=str, help='Multi-sample TRGT VCF of population')
    mainInput.add_argument('--dbs', metavar='DB,DB,...', type=str, nargs='+', help='Databases to merge')
    parser.add_argument('--verbose', action='store_true', default=False, help='Verbose logging')
    parser.add_argument('--logfile', metavar='FILENAME', type=str, help='Write log to FILENAME as well as stdout')
    parser._action_groups.reverse()
    args = parser.parse_args(cmdargs)

    utils.setup_logging(
        args.verbose,
        sys.stderr if args.logfile is None else utils.LogFileStderr(args.logfile),
        show_version=True,
    )

    conn, cursor = dbutils.initialize_db(args.prefix, ['variant'])
    new_columns = {
        'motifs': 'TEXT',
        'struc': 'TEXT',
        'pop_upper': 'INTEGER',
        'pop_n_alleles': 'INTEGER',
        'pop_recessive_upper': 'INTEGER',
        'pop_recessive_n_alleles': 'INTEGER',
        'pop_allele_length_genotypes': 'TEXT',
    }
    dbutils.add_columns(cursor, 'Variant', new_columns)

    if args.dbs:
        # this will drop all columns except those relevant for merging, i.e. csq annotation will disappear
        merge_dbs(conn, cursor, args.dbs)
    else:
        vcf_to_db(cursor, args.vcf)

    conn.commit()
    conn.close()
