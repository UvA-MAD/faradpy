from itertools import islice


def window(seq, n):
    """Return a sliding window of length n from seq"""

    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result
