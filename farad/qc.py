from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np


def plot_read_quality(fqin, plotout):
    """Make a boxplot with quality distribution along the read.

    :param fqin: input reads in fastq file
    :type fqin: file handle or path
    :param plotout: file with generated plot
    :type plotout: file handle or path
    """

    # read in read qualities
    read_quality_strs = [r.letter_annotations['phred_quality']
                         for r in SeqIO.parse(fqin, 'fastq')]
    quals_per_pos = []
    for read_qualities in read_quality_strs:
        for pos, base_quality in enumerate(read_qualities):
            try:
                quals_per_pos[pos].append(base_quality)
            except IndexError:
                quals_per_pos.append([base_quality])

    # make a plot
    ## make ticks every 10 nucleotides
    xticks = np.arange(0, len(quals_per_pos) + 1, 20)
    plt.boxplot(quals_per_pos)
    plt.xticks(xticks)
    plt.title('Phread quality distribution along the read')
    plt.xlabel('Position in the read')
    plt.ylabel('Phred quality score')
    plt.savefig(plotout)
