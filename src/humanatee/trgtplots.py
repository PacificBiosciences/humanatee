"""TRGT plots."""

import itertools
import logging
import math

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.ticker import MaxNLocator

from humanatee import utils

plt.set_loglevel('warning')

# TODO: allow this to be run standalone


class PopHist:
    """Plot allele length with histogram on top and heatmap of genotypes on bottom.

    Inspired by examples:
    https://gnomad.broadinstitute.org/short-tandem-repeat/ATXN10?dataset=gnomad_r4
    https://gnomad.broadinstitute.org/short-tandem-repeat/RFC1?dataset=gnomad_r4
    """

    def __init__(
        self,
        title,
        pop_al,
        filename,
        short_allele_cutoff,
        long_allele_cutoff,
        bin_width=None,
        figsize=None,
        sample_genotypes=[],
        sample_colors=[],
        sample_labels=[],
        pop_color='tab:gray',
        rugplot_color='dimgray',
        cutoff_linestyle=(0, (5, 5)),
        cutoff_linewidth=0.5,
        cutoff_linealpha=0.7,
        sample_marker='x',
        heatmap=True,
        top_xlabel='Length of allele (bp)',
        bottom_xlabel='Length of long allele (bp)',
        bottom_ylabel='Length of short allele (bp)',
    ):
        self.samples = []
        self.axis_labels = {'top_x': top_xlabel, 'bottom_x': bottom_xlabel, 'bottom_y': bottom_ylabel}
        self.cutoff_line = {'style': cutoff_linestyle, 'width': cutoff_linewidth, 'alpha': cutoff_linealpha}
        self.sample_marker = sample_marker
        self.rugplot_color = rugplot_color
        self.pop_color = pop_color
        self.long_allele_cutoff = long_allele_cutoff
        self.short_allele_cutoff = short_allele_cutoff

        if not pop_al:
            logging.warning(f'Cannot plot allele length histogram for {title} because population distribution is None.')
            return

        logging.debug(f'Plotting allele length histogram for {title}')
        if len(sample_colors) != len(sample_genotypes) != len(sample_labels):
            raise ValueError('Length of sample_colors must match length of sample_genotypes and sample_labels.')

        for gt, color, label in zip(sample_genotypes, sample_colors, sample_labels):
            self.samples.append(self.Sample(gt, color, label))

        self.population = Population(
            pop_al, utils.flatten([sample.genotype for sample in self.samples])
        )  # before axis params
        self.axis_params = AxisParams(self.population.min_allele, self.population.max_allele, bin_width)

        # if haploid set up figure to only have simple histogram on top
        self.figsize = [10, 2.5] if not figsize else figsize
        self.fig, self.ax1 = plt.subplots(figsize=self.figsize)

        # if diploid, make heatmap on bottom
        if not self.population.haploid and heatmap:  # if diploid locus
            self.figsize = [8.5, 10] if not figsize else figsize
            self.fig, (self.ax1, self.ax2) = plt.subplots(
                nrows=2, figsize=self.figsize, gridspec_kw={'height_ratios': [1, 3]}, sharex=True
            )

            # ax2: heatmap on bottom with rugplot and colorbar
            self.repeats_heatmap(self.ax2)

        # ax1: histogram on top
        self.repeats_histogram(self.ax1)

        self.fig.suptitle(title[0:45], y=0.9, verticalalignment='bottom')
        logging.debug(f'Making plot: {filename}')
        plt.savefig(filename, bbox_inches='tight')
        plt.close()  # close figure to avoid memory leak

    class Sample:
        """Keep sample info together."""

        def __init__(self, genotype, color, label):
            self.genotype = parse_genotype(genotype)
            self.color = color
            self.label = label

    def repeats_heatmap(self, ax):
        """Plot heatmap of diploid genotypes."""
        half_bin = self.axis_params.bin_width / 2  # used to center items in a bin
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list('', ['#FCFCFC', 'black'], N=1000)
        cmap.set_under(color='white')  # set nan to white

        # vmin=0.0000001 to to make everything that is 0 white
        short_alleles = self.population.diploid_genotypes[0]
        long_alleles = self.population.diploid_genotypes[1]
        h, xedges, yedges, image = ax.hist2d(
            long_alleles,
            short_alleles,
            bins=[
                np.arange(self.axis_params.min, self.axis_params.max, self.axis_params.bin_width),
                np.arange(self.axis_params.min, self.axis_params.max, self.axis_params.bin_width),
            ],
            cmap=cmap,
            vmin=0.0000001,
        )

        ycenters, y_labels = get_ticks(yedges, self.axis_params.bin_width)
        # make sure axis ticks land on center of bin
        # only display some portion of ticks and labels to avoid overcrowding
        ax.set_yticks(ycenters[0 :: self.axis_params.label_gap])
        ax.set_yticklabels(y_labels[:: self.axis_params.label_gap])
        ax.set_xlabel(self.axis_labels['bottom_x'])
        ax.set_ylabel(self.axis_labels['bottom_y'])
        # plot extras
        ax.set_ylim(
            [
                self.population.min_allele - self.axis_params.buffer * 2,
                self.population.max_allele + max(self.axis_params.bin_width, self.axis_params.buffer),
            ]
        )
        ax.text(
            self.population.min_allele - self.axis_params.buffer,
            self.population.max_allele
            + max(self.axis_params.bin_width, self.axis_params.buffer)
            - self.axis_params.buffer,
            f'# Genotypes: {len(short_alleles)}',
            verticalalignment='top',
            horizontalalignment='left',
        )
        if self.short_allele_cutoff:
            ax.axhline(
                y=self.short_allele_cutoff + half_bin,
                color=self.rugplot_color,
                linewidth=self.cutoff_line['width'],
                linestyle=self.cutoff_line['style'],
                alpha=self.cutoff_line['alpha'],
            )
        if self.long_allele_cutoff:
            ax.axvline(
                x=self.long_allele_cutoff + half_bin,
                color=self.rugplot_color,
                linewidth=self.cutoff_line['width'],
                linestyle=self.cutoff_line['style'],
                alpha=self.cutoff_line['alpha'],
            )
        sns.rugplot(
            x=[allele + half_bin for allele in long_alleles],
            y=[allele + half_bin for allele in short_alleles],
            alpha=0.25,
            color=self.rugplot_color,
        )
        for sample in self.samples:
            if len(sample.genotype) == 2:
                plt.scatter(
                    x=max(sample.genotype) + half_bin,
                    y=min(sample.genotype) + half_bin,
                    color=sample.color,
                    marker=self.sample_marker,
                    s=100,
                )
        cbar = self.fig.colorbar(image, ax=ax, orientation='horizontal', aspect=100, fraction=0.02)
        cbar.set_label('Number of individuals')

    def repeats_histogram(self, ax):
        """Plot histogram of allele lengths."""
        h, xedges, patches = ax.hist(
            self.population.alleles,
            bins=np.arange(self.axis_params.min, self.axis_params.max, self.axis_params.bin_width),
            color=self.pop_color,
            alpha=0.25,
        )
        ax.yaxis.set_major_locator(MaxNLocator(integer=True, nbins=4))  # make sure y axis is integers
        handles = []
        ybinsize = (h.max() - h.min()) / max(18, len(self.samples))
        allele_idx = 0
        for sample in self.samples:
            for allele in sample.genotype:
                allele_idx += 1
                ax.axvline(x=allele + self.axis_params.bin_width / 2, color=sample.color, linewidth=1, zorder=1)
                ax.scatter(
                    x=allele + self.axis_params.bin_width / 2,
                    y=ybinsize / 2 + allele_idx * ybinsize,
                    color=sample.color,
                    linewidth=1,
                    zorder=2,
                )
            patch = mpatches.Patch(color=sample.color, label=sample.label)
            handles.append(patch)
            # if (sample_idx + 1) == len(sample_genotypes):
        plt.legend(handles=handles, loc='upper center', ncol=max([1, len(self.samples)]), framealpha=0)
        if self.long_allele_cutoff:
            ax.axvline(
                x=self.long_allele_cutoff + self.axis_params.bin_width / 2,
                color=self.rugplot_color,
                linewidth=self.cutoff_line['width'],
                linestyle=self.cutoff_line['style'],
                alpha=self.cutoff_line['alpha'],
            )
        ax.set_xlabel(self.axis_labels['top_x'])
        ax.set_ylabel('Alleles')
        ax.margins(x=0)  # remove weird histogram padding on extremes of x axis

        xcenters, x_labels = get_ticks(xedges, self.axis_params.bin_width)
        ax.set_xlim(
            [
                self.population.min_allele - self.axis_params.buffer * 2,
                self.population.max_allele + max(self.axis_params.bin_width, self.axis_params.buffer),
            ]
        )
        try:
            self.ax2
            ypadding = 0
        except AttributeError:
            ypadding = ybinsize * 3
        ax.set_ylim([0, max(h) + ybinsize + ypadding])
        ax.text(
            self.population.min_allele - self.axis_params.buffer,
            h.max() + ypadding,
            f'# Alleles: {len(self.population.alleles)}',
            verticalalignment='top',
        )
        ax.text(
            self.population.min_allele - self.axis_params.buffer,
            h.max() - ybinsize * 2 + ypadding,
            f'# Samples: {len(self.population.genotypes)}',
            verticalalignment='top',
        )
        ax.set_xticks(xcenters[:: self.axis_params.label_gap])
        ax.set_xticklabels(x_labels[:: self.axis_params.label_gap])
        ax.xaxis.set_tick_params(which='both', labelbottom=True)  # share ticks but show labels on both


class AxisParams:
    """Set axis limits and ticks based on data."""

    def __init__(self, min_allele, max_allele, bin_width_override):
        data_range = max([1, max_allele - min_allele])
        self.buffer = data_range / 40
        self.bin_width = int(np.ceil(self.buffer)) if not bin_width_override else bin_width_override

        num_bins = int(math.ceil(data_range / self.bin_width))
        max_n_labels = 10 if self.bin_width != 1 else 20  # TODO: determine this number based on fig size?
        self.label_gap = int(math.ceil(num_bins / max_n_labels))
        if max_allele > 1000:
            self.label_gap = self.label_gap + 2  # give more space for larger numbers

        # set same axis for x & y dependent on min and max of data plus a buffer based on bin_width
        self.min, self.max = max(min_allele, 0), max_allele + self.bin_width * 2


class Population:
    """Parse population allele lengths and genotypes."""

    def __init__(self, population_alleles, sample_alleles):
        self.haploid = False

        # pop_al is a string representing a counter dict
        # convert pop_al string to array of arrays where each subarray represents a sample genotype
        self.genotypes = [genotype for genotype in utils.reverse_counter(population_alleles) if set(genotype) != {None}]
        self.alleles = [al for al in utils.flatten(self.genotypes) if al]  # haploid and diploid included
        self.min_allele = np.min(self.alleles + sample_alleles)
        self.max_allele = np.max(self.alleles + sample_alleles)

        self.diploid_genotypes = [g for g in self.genotypes if len(g) == 2]  # for heatmap, can only use diploid
        # transpose so each subarray is long alleles or short alleles
        self.diploid_genotypes = list(map(list, itertools.zip_longest(*self.diploid_genotypes, fillvalue=np.nan)))

        if not self.diploid_genotypes and self.alleles:
            self.haploid = True


def parse_genotype(genotype):
    """Parse genotype string or list into sorted list of integers."""
    # logging.debug(f"genotype: {genotype}")
    if isinstance(genotype, list):
        parsed_genotype = sorted(genotype)
    elif isinstance(genotype, str):
        try:
            parsed_genotype = sorted([int(allele) for allele in genotype.split(',')])
        except ValueError:
            parsed_genotype = []
    else:
        parsed_genotype = []
    return parsed_genotype


def get_ticks(edges, bin_width):
    """Get x or y ticks and labels for histogram."""
    if bin_width > 1:
        # set axis labels to bin ranges
        centers = (edges[:-1] + edges[1:]) / 2
        labels = [f'{int(edges[i])}-{int(edges[i+1])-1}' for i in range(len(edges) - 1)]
    else:
        centers = (edges[:-1] + edges[1:]) / 2
        labels = [f'{int(edges[i])}' for i in range(len(edges) - 1)]
    return centers, labels
