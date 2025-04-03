"""Get population summary statistics from multi-sample TRGT VCF to use for annotation."""

import argparse
import json
import logging
import sqlite3
import sys
from collections import Counter
from importlib import resources

import numpy as np
import pysam

import humanatee


# TODO: print some summary stats to log
# TODO: pull MC from population, get “upper” cutoffs for each motif and summed motifs (don’t use this now, but useful to have)
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

        all_alleles = humanatee.flatten(sorted_genotypes)
        self.lower, self.upper = self.__get_cutoffs__(all_alleles)
        self.cutoff_n = len(all_alleles)

        short_alleles = [next(iter(sample), '') for sample in sorted_genotypes]
        short_alleles = [_ for _ in short_alleles if _ != '']
        self.recessive_lower, self.recessive_upper = self.__get_cutoffs__(short_alleles)
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
            lower, upper = np.percentile(counts, [1, 99])
            lower = int(np.floor(lower))
            upper = int(np.ceil(upper))
        else:
            lower, upper = 0, 0
        return lower, upper


def insert_locus_as_row(cursor, trid, chrom, start, stop, motifs, struc, locus_stats, replace=False):
    """Insert locus into SQLite database."""
    command = 'INSERT OR REPLACE' if replace else 'INSERT'
    cursor.execute(
        f"""
    {command} INTO Locus (
        trid, chrom, start, end, motifs, struc, lower, upper, n_alleles,
        recessive_lower, recessive_upper, recessive_n_alleles, allele_length_genotypes
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """,
        (
            trid,
            chrom,
            start,
            stop,
            motifs,
            struc,
            locus_stats.lower,
            locus_stats.upper,
            locus_stats.cutoff_n,
            locus_stats.recessive_lower,
            locus_stats.recessive_upper,
            locus_stats.recessive_cutoff_n,
            json.dumps(locus_stats.al_genotype_dict),
        ),
    )


def vcf_to_db(vcf, cohort, prefix):
    """Convert multi-sample VCF to DB with aggregate statistics."""
    v = pysam.VariantFile(vcf)
    samples = list(v.header.samples)
    logging.info(f'Samples: {samples}')

    # initialize output database, removing if it already exists
    humanatee.silent_remove(f'{prefix}.db')
    conn = sqlite3.connect(f'{prefix}.db')
    cursor = conn.cursor()

    # read schema from file
    vcf_schema = resources.files('humanatee').joinpath('schemas', 'trgtpop.sql').read_text()
    cursor.executescript(vcf_schema)

    # add samples to Sample table
    for sample in samples:
        cursor.execute(
            """
        INSERT INTO Sample (sampleid, cohort) VALUES (?, ?)
        """,
            (sample, cohort),
        )

    # add records to Locus table
    for record in v:
        trid = record.info['TRID']
        if isinstance(trid, tuple):
            trid = ','.join(record.info['TRID'])
        logging.debug(f'TRID: {trid}')
        al_genotypes = [record.samples[sample]['AL'] for sample in samples]
        sd_genotypes = [record.samples[sample]['SD'] for sample in samples]
        locus_stats = LocusStats(al_genotypes, sd_genotypes)
        insert_locus_as_row(
            cursor,
            trid,
            record.chrom,
            record.start,
            record.stop,
            ','.join(record.info['MOTIFS']),
            record.info['STRUC'],
            locus_stats,
        )

    # write db to tsv
    # results = dbutils.smart_execute(cursor, "SELECT * FROM Locus").fetchall()
    # with open(f'{prefix}.tsv', 'w', newline='') as f:
    #     writer = csv.writer(f, delimiter='\t')
    #     writer.writerow([_[0] for _ in cursor.description])  # header
    #     writer.writerows(results)

    v.close()

    # Commit and close the connection
    conn.commit()
    conn.close()


def merge_dbs(db_list, prefix):
    """Merge list of popDBs into a single DB. Expects sample/cohort combinations to be unique in each DB."""
    logging.info(f'Merging databases: {db_list}')

    def merge_dbs_helper(conn, cursor, db2):
        """Merge second database into cursor."""
        db2_conn = sqlite3.connect(db2)
        db2_cursor = db2_conn.cursor()
        humanatee.validate_db_schema(conn, cursor, db2_conn, db2)

        # update Samples table
        db2_cursor.execute('SELECT * FROM Sample')
        samples = db2_cursor.fetchall()
        for row in samples:
            try:
                cursor.execute('INSERT INTO Sample (sampleid, cohort) VALUES (?, ?)', (row))
            except sqlite3.IntegrityError:
                cursor.execute('SELECT * FROM Sample WHERE sampleid = ?', (row[0],))
                original_cohort = cursor.fetchone()[1]
                logging.warning(
                    f"Sample '{row[0]}' already exists in database associated with cohort '{original_cohort}' and will not be added."
                )
                pass

        # update Locus table
        db2_cursor.execute('SELECT * FROM Locus')
        column_index = list(map(lambda x: x[0], db2_cursor.description)).index('allele_length_genotypes')
        for row in db2_cursor.fetchall():
            try:
                cursor.execute(f"INSERT INTO Locus VALUES ({', '.join(['?' for _ in row])})", row)
            except sqlite3.IntegrityError:
                # handle primary key conflict by replace row
                cursor.execute('SELECT allele_length_genotypes FROM Locus where trid = ?', (row[0],))
                results = cursor.fetchone()
                al_genotypes = humanatee.reverse_counter(row[column_index]) + humanatee.reverse_counter(results[0])
                locus_stats = LocusStats(al_genotypes)
                insert_locus_as_row(cursor, row[0], row[1], row[2], row[3], row[4], row[5], locus_stats, replace=True)
        db2_conn.close()

    # initialize output database, removing if it already exists
    humanatee.silent_remove(f'{prefix}.db')
    conn = sqlite3.connect(f'{prefix}.db')
    cursor = conn.cursor()
    vcf_schema = resources.files('humanatee').joinpath('schemas', 'trgtpop.sql').read_text()
    cursor.executescript(vcf_schema)

    # merge databases one by one
    for db in db_list:
        merge_dbs_helper(conn, cursor, db)
        conn.commit()

    # close the connection
    conn.close()


def trgtpop_main(cmdargs):
    """Run from command line."""
    parser = argparse.ArgumentParser(description=__doc__, prog='humanatee trgtpop')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('--prefix', metavar='PREFIX', required=True, type=str, help='Prefix for output')

    mainInput = parser.add_mutually_exclusive_group(required=True)
    mainInput.add_argument('--vcf', metavar='VCF', type=str, help='Multi-sample TRGT VCF of population')
    mainInput.add_argument('--dbs', metavar='DB,DB,...', type=str, nargs='+', help='Databases to merge')
    parser.add_argument('--cohort-id', metavar='TEXT', type=str, help='Cohort name')
    parser.add_argument('--verbose', action='store_true', default=False, help='Verbose logging')
    parser.add_argument('--logfile', metavar='FILENAME', type=str, help='Write log to FILENAME as well as stdout')
    parser._action_groups.reverse()
    args = parser.parse_args(cmdargs)

    if args.vcf and (args.cohort_id is None):
        parser.error('The --vcf parameter requires --cohort-id.')

    humanatee.setup_logging(
        args.verbose,
        sys.stderr if args.logfile is None else humanatee.LogFileStderr(args.logfile),
        show_version=True,
    )

    if args.dbs:
        merge_dbs(args.dbs, args.prefix)
    else:
        vcf_to_db(args.vcf, args.cohort_id, args.prefix)
