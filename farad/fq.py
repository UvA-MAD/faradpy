from Bio import SeqIO
from .utils import window
import logging
import random


# setup module logger
logger = logging.getLogger(__name__)
FORMAT = '%(asctime)s - %(levelname)s - %(message)s'
# logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
formatter = logging.Formatter(FORMAT)
handler.setFormatter(formatter)
logger.addHandler(handler)


def downsample(fqin, N, fqout):
    """Downsample fastq file to N reads.

    Sample rendomly N reads from fqin file
    and write results in fqout.

    :param fqin: input fastq file
    :type fqin: file handle or path
    :param N: sample size
    :type N: int.
    :param fqout: output fastq file
    :type fqout: file handle or path
    :raises: FileNotFoundError, ValueError, TypeError
    """

    # read in quality reads into a list
    try:
        input_reads = [r for r in SeqIO.parse(fqin, 'fastq')]
    except FileNotFoundError:
        logger.error('File with input reads not found')
        raise

    # fetch a sample
    try:
        sample_reads = random.sample(input_reads, N)
    except ValueError:
        logger.error('Sample size should be between zero'
                     ' and the number of reads in fq file.')
        raise
    except TypeError:
        logger.error('Sample size should be an integer.')
        raise

    # write sample to file fqout
    try:
        SeqIO.write(sample_reads, fqout, 'fastq')
    except FileNotFoundError:
        logger.error('Provided path for sample reads is not valid.')
        raise


def quality_trim_read(sr, window_len, min_qual):
    """Trim read with the moving window from the end.

    Trim one base at a time untill all bases in the window
    have quality equal or above ``min_qual``

    :param sr: Sequencing read to be trimmed.
    :type sr: Bio.SeqRecord.SeqRecord.
    :param window_len: lenght of the sliding window
    :type window_len: int.
    :param min_qual: Minimal quality used for trimming.
    :type min_qual: int.
    :returns: Bio.SeqRecord.SeqRecord -- trimmed sequencing read.
    """
    # get quality values of the read
    quals = sr.letter_annotations['phred_quality']

    # invert it to use moving window going from the end
    inv_quals = quals[::-1]
    n_trimmed = 0
    for w in window(inv_quals, window_len):
        if all([q >= min_qual for q in w]):
            break
        else:
            n_trimmed += 1
    if n_trimmed > 0:
        trimmed = sr[:-n_trimmed]
    else:
        trimmed = sr
    return(trimmed)


def quality_trim_fastq(fastq_in, fastq_out, window_len, min_qual):
    """Trim all reads in provided fastq file.

    Iterate over all reads in the file.
    Trim trailing bases which using sliding window.

    :param fastq_in: Input fastq file.
    :type fastq_in: File handle or path.
    :param fastq_out: Output fastq file.
    :type fastq_out: File handle or path.
    :param window_len: Size of sliding window.
    :type window_len: int.
    """
    # parse input fastq
    seqs = SeqIO.parse(fastq_in, 'fastq')
    trimmed_seqs = [quality_trim_read(s, window_len, min_qual)
                    for s in seqs]
    SeqIO.write(trimmed_seqs, fastq_out, 'fastq')


def filter_by_length_fastq(fastq_in, fastq_out, min_len=1, max_len=None):
    """Fetch only reads of length between min_len and max_len.

    Can be used to either fetching the reads of lentght between min_len
    and max_len or by specifying only one of the borders to fetch reads longer
    or shoter than specified value.

    :param fastq_in: Input fastq file.
    :type fastq_in: File handle or path.
    :param fastq_out: Output fastq file.
    :type fastq_out: File handle or path.
    :param min_len: Minimal length of the read. Defaults to 1.
    :type min_len: int.
    :param max_len: Maximal length of the read. Defaults to length of the read.
    """
    # get an iterator from the reads file
    seqs = SeqIO.parse(fastq_in, 'fastq')

    # if min_len and max_len have not been specified
    # write reads as they are
    if (min_len == 1 and max_len is None):
        filt_seqs = [s for s in seqs]
    elif (max_len is None):
        # filter longer or equal to min_len
        filt_seqs = [s for s in seqs if len(s) >= min_len]
    else:
        # filter reads within a range between min and max _len
        filt_seqs = [s for s in seqs
                     if (len(s) >= min_len and len(s) <= max_len)]
    SeqIO.write(filt_seqs, fastq_out, 'fastq')
