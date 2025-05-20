"""Get reference information for BED/VCF.
Much of this is adapted shamelessly from https://github.com/PacificBiosciences/svpack.
"""

import argparse
import collections
import heapq
import logging
import os
import shutil
import sqlite3
import sys

import pysam

from humanatee import dbutils, utils


class Variant:
    """A genomic location."""

    def __init__(self, locus_id, source, chrom, chrom_idx, start, stop, misc_dict={}):
        self.id = locus_id
        self.source = source
        self.chrom = chrom
        self.chrom_idx = chrom_idx
        self.start = int(start)  # TRGT VCFs var.pos equals TRGT BED start
        self.stop = int(stop)
        self.misc = misc_dict

    def update_row(self, misc_dict):
        """Update the consequence of the repeat."""
        self.misc |= misc_dict

    def insert_csq(self, cursor, csq_list):
        """Update the consequence of the repeat."""
        if csq_list is None:
            return
        for csq in csq_list:
            csq.update({'variant_id': self.id, 'source': self.source, 'allele_index': None})
            dbutils.add_table_row(
                cursor, 'Consequence', csq, update_conflict='variant_id, source, ensembl_id', replace=False
            )

    def insert_variant(self, cursor):
        """Insert the repeat locus into the database."""
        row_dict = {
            'variant_id': self.id,
            'source': self.source,
            'chrom': self.chrom,
            'start': self.start,
            'end': self.stop,
        }
        row_dict |= self.misc
        dbutils.add_table_row(cursor, 'Variant', row_dict, update_conflict='variant_id, source', replace=False)


class ParsedVariants:
    """Stream of sorted variants."""

    def __init__(self, vcf, func, sort):
        logging.info(f'Loading VCF file: {vcf}')
        self._variant_file = pysam.VariantFile(vcf)
        self.func = func
        self.sort = sort
        self.chrom_sort_order = None
        if not sort:
            try:
                self.chrom_sort_order = self._variant_file.header.contigs.keys()
            except AttributeError:
                logging.warning('VCF file does not have contig information for sort order. Sorting will be used.')
                self.sort = True
        if sort:
            variants = [func(record, self.chrom_sort_order) for record in self._variant_file]
            sorted_variants = list(sorted(variants, key=lambda x: (x.chrom_idx, x.start, x.stop)))
            self._variant_file.close()
            self._variants = iter(sorted_variants)
        else:
            self._variants = self._variant_file

    def __iter__(self):
        """Return the iterator object (self)."""
        return self

    def __next__(self):
        """Return the next element of the sequence or raises StopIteration if exhausted."""
        if self.sort:
            return next(self._variants)
        else:
            return self.func(next(self._variants), self.chrom_sort_order)

    def close(self):
        """Close the VCF file."""
        if hasattr(self._variants, 'close'):
            self._variants.close()


class AnnoSource:
    """Abstract class for annotation sources."""

    def __init__(self, label, features, distance, func):
        """Initialize annotation source."""
        self.label = label
        self.features = collections.deque(features)  # iterable of objects with (chrom,start,stop), should be sorted
        self.distance = distance
        self.func = func
        self.deque = None
        self.activefeatures = []
        self.serialnumber = 0  # required as a tie-breaker for heapq if chrom and end are the same

    def update_active_features(self, target_variant):
        """Update active features for current repeat."""
        # add genes that start before the current repeat ends to the active heap
        while len(self.features) and (self.features[0].chrom_idx, self.features[0].start) <= (
            target_variant.chrom_idx,
            target_variant.stop + self.distance,
        ):
            feature = self.features.popleft()
            heapq.heappush(self.activefeatures, ((feature.chrom_idx, feature.stop, self.serialnumber), feature))
            self.serialnumber += 1
        # pop genes that end before the current repeat starts, can't overlap any later repeats
        # WARNING: requires that both genes and repeats be sorted!
        while len(self.activefeatures) and self.activefeatures[0][0] < (
            target_variant.chrom_idx,
            target_variant.start - self.distance,
        ):
            heapq.heappop(self.activefeatures)


class Bed:
    """BED file to be used as an annotation source."""

    class BedLine:
        """BED line."""

        def __init__(self, chrom_idx, start, stop):
            # logging.debug(f'Processing BED line: {line}')
            self.chrom_idx = chrom_idx
            self.start = int(start) + 1  # TODO: confirm
            self.stop = int(stop)

    def __init__(self, bedfile, chrom_sort_order=None):
        logging.info(f'Loading BED file: {bedfile}')
        bedlines = list()
        f = utils.open_file(bedfile)
        for line in f:
            line = line.rstrip('\n')
            if line == '' or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            try:
                chrom_idx = chrom_sort_order.index(fields[0]) if chrom_sort_order else fields[0]
                bedlines.append(self.BedLine(chrom_idx, fields[1], fields[2]))
            except ValueError:
                # chrom not in VCF header, so ignore
                continue

        f.close()

        sorted_lines = list(sorted(bedlines, key=lambda x: (x.chrom_idx, x.start, x.stop)))

        # connect overlapping intervals
        # TODO: compare to bedtools
        # TODO: maybe don't merge overlapping intervals?, just get overlap and list labels overlapping from 4th column
        # TODO: or alternatively, merge overlapping intervals but group by 4th column label first
        merged = [sorted_lines[0]]
        for line in sorted_lines[1:]:
            last_merged = merged[-1]
            if line.chrom_idx == last_merged.chrom_idx and line.start <= last_merged.stop:
                # overlap, merge
                merged[-1] = self.BedLine(last_merged.chrom_idx, last_merged.start, max(last_merged.stop, line.stop))
            else:
                # no overlap, add to merged
                merged.append(line)
        self.lines = merged


class Gff:
    """GFF file to be used as an annotation source."""

    class GffLine:
        """GFF line."""

        def __init__(self, line, chrom_idx):
            logging.debug(f'Processing GFF line: {line}')
            chrom, source, feature, start, end, score, strand, phase, attributes = line.rstrip('\n').split('\t')[:9]
            self.chrom_idx = chrom_idx
            self.source = source
            self.feature = feature
            self.start = int(start)  # TODO: make sure you don't have to subtract 1 from start
            self.stop = int(end)
            self.score = None if score == '.' else float(score)
            self.strand = strand
            self.phase = phase
            self.attributes = dict()
            for attribute in attributes.split(';'):
                k, v = attribute.split('=')
                self.attributes[k] = v

    class GffGene:
        """GFF lines that are genes; child=transcript."""

        def __init__(self, gff_line):
            self.chrom_idx = gff_line.chrom_idx
            self.source = gff_line.source
            self.feature = gff_line.feature
            self.start = gff_line.start
            self.stop = gff_line.stop
            self.score = gff_line.score
            self.strand = gff_line.strand
            self.phase = gff_line.phase
            self.attributes = gff_line.attributes
            self.id = self.attributes['ID'].split(':')[1] if ':' in self.attributes['ID'] else self.attributes['ID']
            self.gene_id = self.attributes['gene_id']
            try:
                self.symbol = self.attributes['gene_name']
            except KeyError:
                self.symbol = self.attributes['gene_id']

            self.transcripts = []

        def add_transcript(self, gff_transcript):
            """Add child transcript to gene."""
            self.transcripts.append(gff_transcript)

    class GffTranscript:
        """GFF lines that are transcripts; parent=gene, child=exon."""

        def __init__(self, gff_line):
            self.chrom_idx = gff_line.chrom_idx
            self.source = gff_line.source
            self.feature = gff_line.feature
            self.start = gff_line.start
            self.stop = gff_line.stop
            self.score = gff_line.score
            self.strand = gff_line.strand
            self.phase = gff_line.phase
            self.attributes = gff_line.attributes
            self.id = self.attributes['ID'].split(':')[1] if ':' in self.attributes['ID'] else self.attributes['ID']
            self.gene = (
                self.attributes['Parent'].split(':')[1]
                if ':' in self.attributes['Parent']
                else self.attributes['Parent']
            )
            self.exons = []

        def add_exon(self, gff_exon):
            """Add child exon to transcript."""
            self.exons.append(gff_exon)

    class GffExon:
        """GFF lines that are exons, parent=transcript."""

        def __init__(self, gff_line):
            self.chrom_idx = gff_line.chrom_idx
            self.source = gff_line.source
            self.feature = gff_line.feature
            self.start = gff_line.start
            self.stop = gff_line.stop
            self.score = gff_line.score
            self.strand = gff_line.strand
            self.phase = gff_line.phase
            self.attributes = gff_line.attributes
            self.transcript = (
                self.attributes['Parent'].split(':')[1]
                if ':' in self.attributes['Parent']
                else self.attributes['Parent']
            )

    def __init__(self, gff_file, chrom_sort_order=None):
        logging.info(f'Loading GFF file: {gff_file}')
        self.chrom_sort_order = chrom_sort_order
        self.genes = dict()
        self.transcripts = dict()
        self.exons = list()

        f = utils.open_file(gff_file)
        for line in f:
            line = line.rstrip('\n')
            if line == '' or line.startswith('#'):
                continue
            self.__process_gff_line__(line)
        f.close()

        # link exons to transcripts
        for exon in self.exons:
            if exon.transcript in self.transcripts:
                self.transcripts[exon.transcript].add_exon(exon)
        # link transcripts to genes
        for transcript in self.transcripts.values():
            self.genes[transcript.gene].add_transcript(transcript)

        # provide genes sorted by (chrom,start,end)
        self.genes = list(sorted(self.genes.values(), key=lambda x: (x.chrom_idx, x.start, x.stop)))

        # group genes by chromosome, sorting by (start,end) within a chromosome
        self.genes_by_chrom = collections.defaultdict(list)
        for gene in self.genes:
            self.genes_by_chrom[gene.chrom_idx].append(gene)  # could be used for parallelizing later

    def __process_gff_line__(self, line):
        chrom_idx = line.rstrip('\n').split('\t')[0]
        if self.chrom_sort_order:
            try:
                chrom_idx = self.chrom_sort_order.index(chrom_idx)
            except ValueError:
                # chrom not in VCF header, so ignore
                return
        gff_line = Gff.GffLine(line, chrom_idx)
        # features = gene, transcript, CDS, five_prime_UTR, three_prime_UTR,
        # exon, start_codon, stop_codon, stop_codon_redefined_as_selenocysteine
        if gff_line.feature == 'gene':
            gff_gene = Gff.GffGene(gff_line)
            self.genes[gff_gene.id] = gff_gene
        elif gff_line.feature == 'transcript':
            gff_transcript = Gff.GffTranscript(gff_line)
            # only retain protein-coding transcripts
            if gff_transcript.attributes['transcript_type'] == 'protein_coding':
                self.transcripts[gff_transcript.id] = gff_transcript
        # https://www.cell.com/cms/10.1016/j.xgen.2023.100296/asset/7a873e39-c366-4186-8c04-cd41e1381b4d/main.assets/gr1_lrg.jpg
        elif gff_line.feature in ('five_prime_UTR', 'three_prime_UTR', 'CDS'):
            gff_exon = Gff.GffExon(gff_line)
            self.exons.append(gff_exon)


def get_bp_distance(gene_start, gene_end, bed_start, bed_end, strand=None):
    """Get the distance between a gene and another genomic feature."""
    if gene_start > gene_end or bed_start > bed_end:
        raise ValueError('Gene start/end or bed start/end are invalid.')
    if gene_start > bed_end:
        distance = bed_end - gene_start
    elif gene_end < bed_start:
        distance = bed_start - gene_end
    else:
        distance = 0
    if strand == '-':
        distance = -distance
    if strand is None:
        distance = abs(distance)
    return distance


def annotate_variants(cursor, variants, anno_sources):
    """Get gene consequences for all loci in TRGT BED or VCF."""
    logging.info('Annotating variants')
    prev = None
    for variant in variants:
        if prev and (variant.chrom_idx, variant.start, variant.stop) < (prev.chrom_idx, prev.start, prev.stop):
            logging.debug(
                f'prev: {prev.chrom_idx}:{prev.start}-{prev.stop} current: {variant.chrom_idx}:{variant.start}-{variant.stop}'
            )
            raise ValueError('Variants are not sorted. Use --sort option to sort them.')
        logging.debug(f'Processing variant: {variant.id}')
        for anno_source in anno_sources:
            anno_source.update_active_features(variant)
            update_value = anno_source.func(anno_source.activefeatures, variant)
            if anno_source.func == get_csq:
                # CSQ is a special case, it is a list of consequences
                variant.insert_csq(cursor, update_value)
            else:
                variant.update_row({anno_source.label: update_value})
        variant.insert_variant(cursor)
        prev = variant


def parse_trgt_variant(record, chrom_sort_order=None):
    """Parse a TRGT VCF record."""
    chrom_idx = chrom_sort_order.index(record.chrom) if chrom_sort_order else record.chrom
    variant_id = record.info['TRID']
    if isinstance(variant_id, tuple):
        variant_id = ','.join(variant_id)
    return Variant(
        variant_id,
        'trgt',
        record.chrom,
        chrom_idx,
        record.pos,
        record.stop,
        {
            'motifs': ','.join(record.info['MOTIFS']),
            'struc': record.info['STRUC'],
        },
    )


class Consequence:
    """Consequence of a repeat on a gene."""

    def __init__(self, gene_name, gene_id, gene_type, strand, features):
        self.gene = gene_name
        self.gene_id = gene_id
        self.gene_type = gene_type
        self.strand = strand
        self.features = features

    def __str__(self):
        return '|'.join([self.gene, self.gene_id, self.gene_type, self.strand, ','.join(self.features)])


def intersect_intervals(A, B):
    """Get intersection of two sets of sorted, disjoint, closed intervals."""
    # A and B are lists of closed intervals [start, end]
    # A and B contain disjoint intervals (non-overlapping)
    # A and B are sorted by start position
    i = 0
    j = 0

    result = []
    while i < len(A) and j < len(B):
        a_start, a_end = A[i]
        b_start, b_end = B[j]
        if a_start <= b_end and b_start <= a_end:  # overlap exists
            result.append([max(a_start, b_start), min(a_end, b_end)])

        if a_end <= b_end:  # move to next interval in A
            i += 1
        else:
            j += 1  # move to next interval in B
    return result


def get_overlap_or_distance(activefeatures, variant):
    """Identify bp overlap or distance between variant and features."""
    if not activefeatures:
        return None
    intersection = intersect_intervals(
        [[variant.start, variant.stop]], [[feature.start, feature.stop] for k, feature in activefeatures]
    )
    overlap = sum([stop - start + 1 for start, stop in intersection]) if intersection else 0
    if overlap > 0:
        overlap_percent = int(overlap / (variant.stop - variant.start + 1) * 100)
        return f'overlap={overlap}bp_{overlap_percent}%'
    else:
        distance = min(
            get_bp_distance(feature.start, feature.stop, variant.start, variant.stop) for k, feature in activefeatures
        )
        return f'distance={distance}bp'


def get_csq(activegenes, variant):
    """Identify gene consequences for a repeat."""
    if not activegenes:
        return None
    consequence_order = ['CDS', 'five_prime_UTR', 'three_prime_UTR', 'intron', 'upstream', 'downstream']
    consequences = []
    for k, gene in activegenes:
        # check if repeat overlaps gene body
        if gene.chrom_idx != variant.chrom_idx:
            continue
        features = set()
        distance = 0
        if min(variant.stop, gene.stop) >= max(variant.start, gene.start):  # overlap
            # check if repeat overlaps any exon (defined here as cds, 5_prime_utr, 3_prime_utr) of the gene
            for transcript in gene.transcripts:
                for exon in transcript.exons:
                    if min(variant.stop, exon.stop) > max(variant.start, exon.start):
                        features.add(exon.feature)
            # if it overlaps a protein coding gene (i.e. one with transcripts recorded here),
            # but doesn't overlap an exon, then it overlaps an intron
            if not features and gene.transcripts:
                features.add('intron')
        else:
            # check if repeat is upstream or downstream of gene
            distance = get_bp_distance(gene.start, gene.stop, variant.start, variant.stop, gene.strand)
            features.add(f'upstream={abs(distance)}bp' if distance < 0 else f'downstream={abs(distance)}bp')
        sorted_features = sorted([f.split('=')[0] for f in features], key=consequence_order.index)
        ranked_feature = (
            f'{consequence_order.index(sorted_features[0])+1}_{sorted_features[0]}' if sorted_features else None
        )
        consequences.append(
            {
                'gene': gene.symbol,
                'ensembl_id': gene.gene_id,
                'ranked_consequence': ranked_feature,
                'ranked_impact': None,
                'csq': str(
                    Consequence(gene.symbol, gene.gene_id, gene.attributes.get('gene_type', ''), gene.strand, features)
                ),
            }
        )
    # csq = ';'.join([str(c) for c in consequences])
    return consequences


class ValidateAnnotation(argparse.Action):
    """Argparse class to validate annotation input."""

    def __call__(self, parser, namespace, values, option_string=None):
        """Validate annotation input."""
        items = getattr(namespace, self.dest, None) or []
        file, label, distance = values
        if not os.path.exists(file):
            raise argparse.ArgumentTypeError(f'Annotation file must exist: {file}.')
        if not any([file.endswith(ext) for ext in ['.gff', '.gff.gz', '.gff3.gz', '.gff3', '.bed', '.bed.gz']]):
            raise argparse.ArgumentTypeError(f'Unsupported annotation file format. Must be GFF or BED: {file}.')
        int(distance)
        Anno = collections.namedtuple('Anno', 'bed label distance')
        items.append(Anno(file, label, int(distance)))
        setattr(namespace, self.dest, items)


def refanno_main(cmdargs):
    """Run from command line."""
    parser = argparse.ArgumentParser(description=__doc__, prog='humanatee trgtref')
    requiredNamed = parser.add_argument_group('required arguments')

    parser.add_argument('--db', metavar='DB', type=str, help='Existing SQLite database to update')
    requiredNamed.add_argument('--prefix', metavar='PREFIX', type=str, required=True, help='Prefix for output')
    requiredNamed.add_argument('--vcf', metavar='BED/VCF', required=True, type=str, help='TRGT VCF')
    requiredNamed.add_argument(
        '--source', metavar='SOURCE', required=True, type=str, default='trgt', help='VCF source, e.g. trgt'
    )
    requiredNamed.add_argument(
        '--annotation',
        nargs=3,
        required=True,
        action=ValidateAnnotation,
        metavar=('FILE', 'LABEL', 'DISTANCE'),
        help='Additional BED files to load',
    )
    parser.add_argument(
        '--sort',
        action='store_true',
        default=False,
        help='Sort variants by chrom/start/stop, not recommended for large VCFs',
    )
    parser.add_argument('--verbose', action='store_true', default=False, help='Verbose logging')
    parser.add_argument('--logfile', metavar='FILENAME', type=str, help='Write log to FILENAME as well as stdout')
    parser._action_groups.reverse()
    args = parser.parse_args(cmdargs)
    utils.setup_logging(
        args.verbose,
        sys.stderr if args.logfile is None else utils.LogFileStderr(args.logfile),
        show_version=True,
    )

    if args.db:
        logging.info(f'Copying existing database: {args.db}')
        temp_conn = sqlite3.connect(args.db)
        dbutils.validate_db_tables(temp_conn, {'Variant': ['variant_id', 'source', 'chrom', 'start', 'end']})
        temp_conn.close()
        shutil.copyfile(args.db, f'{args.prefix}.db')
        conn, cursor = dbutils.initialize_db(args.prefix, ['variant', 'consequence'], delete=False)
    else:
        conn, cursor = dbutils.initialize_db(args.prefix, ['variant', 'consequence'], delete=True)

    if args.source.lower() == 'trgt':
        dbutils.add_columns(cursor, 'Variant', {'motifs': 'TEXT', 'struc': 'TEXT'})
        variants = ParsedVariants(args.vcf, parse_trgt_variant, args.sort)
    else:
        raise ValueError('Unsupported VCF source. Only "trgt" is currently supported.')

    anno_sources = []
    for file, label, distance in args.annotation:
        if any([file.endswith(ext) for ext in ['.gff', '.gff.gz', '.gff3.gz', '.gff3']]):
            if label != 'csq':
                raise ValueError('Annotation source label must be "csq" for GFF files.')
            genes = Gff(file, chrom_sort_order=variants.chrom_sort_order)
            anno_sources.append(AnnoSource('csq', genes.genes, distance, get_csq))
        elif any([file.endswith(ext) for ext in ['.bed', '.bed.gz']]):
            bed = Bed(file, chrom_sort_order=variants.chrom_sort_order)
            dbutils.add_columns(cursor, 'Variant', {label: 'TEXT'})
            anno_sources.append(AnnoSource(label, bed.lines, distance, get_overlap_or_distance))

    annotate_variants(cursor, variants, anno_sources)
    conn.commit()
    conn.close()
