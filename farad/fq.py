# from Bio import SeqIO
from utils import window


def quality_trim_read(sr, window_len, min_qual):
    """Trim read with the moving window from the end.

    Trim one base at a time untill all bases in the window
    have quality equal or above ``min_qual``

    :param sr: Sequencing read to be trimmed.
    :type name: Bio.SeqRecord.SeqRecord.
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