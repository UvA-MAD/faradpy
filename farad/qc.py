from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import logging

# setup module logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def plot_qualities_along_read(fqin, plotout):
    """Make a boxplot with quality distribution along the read.

    :param fqin: input reads in fastq file
    :type fqin: file handle or path
    :param plotout: file with generated plot
    :type plotout: file handle or path
    """

    # get read qualities
    try:
        read_quality_strs = [r.letter_annotations['phred_quality']
                             for r in SeqIO.parse(fqin, 'fastq')]
    except FileNotFoundError:
        log.error('Input %s fastq file not found' % fqin)
        raise

    # preallocate the list of length of longest read
    max_rlen = max(len(r) for r in read_quality_strs)
    quals_per_pos = [[] for _ in range(max_rlen)]

    for read_qualities in read_quality_strs:
        for pos, base_quality in enumerate(read_qualities):
            quals_per_pos[pos].append(base_quality)

    # make a plot
    ## make ticks every 10 nucleotides
    xticks = np.arange(0, len(quals_per_pos) + 1, 20)
    plt.boxplot(quals_per_pos)
    plt.xticks(xticks)
    plt.title('Phread quality distribution along the read')
    plt.xlabel('Position in the read')
    plt.ylabel('Phred quality score')
    try:
        plt.savefig(plotout)
    except FileNotFoundError:
        log.error('Output filepath "%s" is not valid.' % plotout)
        raise


def plot_mean_quality_distribution(fqin, plotout):
    """Make a histogram mean quality of the reads

    :param fqin: input reads in fastq file
    :type fqin: file handle or path
    :param plotout: file with generated plot
    :type plotout: file handle or path
    """

    # get read qualities
    try:
        mean_read_qualities = [np.mean(r.letter_annotations['phred_quality'])
                               for r in SeqIO.parse(fqin, 'fastq')]
    except FileNotFoundError:
        log.error('Input %s fastq file not found' % fqin)
        raise

    nbins = np.ceil(max(mean_read_qualities) - min(mean_read_qualities))
    plt.hist(mean_read_qualities, nbins)
    plt.title('Mean read phred quality distribution.')
    plt.xlabel('Mean read quality')
    plt.ylabel('Number of reads')
    try:
        plt.savefig(plotout)
    except FileNotFoundError:
        log.error('Output filepath "%s" is not valid.' % plotout)
        raise
