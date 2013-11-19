# from Bio import SeqIO
from utils import window


def quality_trim_read(sr, window_len, min_qual):
    """Trim read with the moving window from the end.

    Trim one base at a time untill all bases in the window
    have quality equal or above ``min_qual``

    :param sr: BioPython's ``SeqRecord`` object
    :param window_len: lenght of the sliding window
    :param min_qual: minimal quality used for trimming
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
