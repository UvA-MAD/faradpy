"""Tests for farad."""
import os
import re
import sys
import Bio
from Bio import SeqIO
import nose.tools as nt

# import tested package
sys.path.append('../')
# import farad.utils as futils
import farad.fq as ffq


# test data
TEST_DIR = os.path.split(os.path.abspath(__file__))[0]
SAMPLE_FQ1 = os.path.join(TEST_DIR, 'sample1.fq')


class TestDownsample:
    """Test downsample function in fq module"""

    def setup(self):
        """Create objects needed for testing."""
        pass

    def teardown(self):
        """Teardown.

        Remove files created during the test.
        """

        # all files created during test should start with 'test_'
        file_pattern = r'^test_.*\.fq$'
        for_removal = [f for f in os.listdir(TEST_DIR)
                       if re.search(file_pattern, f)]
        [os.remove(os.path.join(TEST_DIR, f)) for f in for_removal]

    def test_raised_error_when_invalid_input_fastq(self):
        """Should raise an exeption when input fastq not valid"""
        fqin = 'nonexistent.fq'
        N = 10
        fqout = os.path.join(TEST_DIR, 'test_invalid.fq')
        nt.assert_raises(FileNotFoundError, ffq.downsample, fqin, N, fqout)

    def test_raised_error_when_sample_larger_than_input(self):
        """Should raise an exception if sample larger than input."""

        N = 15
        fqout = os.path.join(TEST_DIR, 'test_downsample.fq')
        nt.assert_raises(ValueError,
                         ffq.downsample, SAMPLE_FQ1, N, fqout)

    def test_error_if_N_not_integer(self):
        """Should raise an exception if N is not an integer."""
        N = 'not_an_int'
        fqout = os.path.join(TEST_DIR, 'test_sample_not_int.fq')
        nt.assert_raises(TypeError,
                         ffq.downsample, SAMPLE_FQ1, N, fqout)

    def test_error_when_output_path_invalid(self):
        """ Should raise error when invalid output path."""
        N = 1
        fqout = os.path.join(TEST_DIR, 'invalid', 'test_invalid_path.fq')
        nt.assert_raises(FileNotFoundError,
                         ffq.downsample, SAMPLE_FQ1, N, fqout)


class TestQualityTrimRead:
    """ Test quality_trim_read function in fq module"""

    # read in the reads from fastq as list
    sreads = [r for r in SeqIO.parse(SAMPLE_FQ1, 'fastq')]
    # trimming settings
    window_len = 3
    min_qual = 21

    def test_length_corectness(self):
        """Test if reads are trimmed correctly"""

        # first read of len 82 expected last 2 bases to be trimmed
        # with min_qual=21
        trimmed_sread1 = ffq.quality_trim_read(self.sreads[0],
                                               self.window_len,
                                               self.min_qual)
        nt.assert_equal(len(trimmed_sread1), 80)

        # second read of len 29 expected to get one base trimmed
        trimmed_sread2 = ffq.quality_trim_read(self.sreads[1],
                                               self.window_len,
                                               self.min_qual)
        nt.assert_equal(len(trimmed_sread2), 28)

    def test_return_value(self):
        """test if an object returned by quality_trim_read is
        of class Bio.SeqRecord.SeqRecord"""

        trimmed_sread = ffq.quality_trim_read(self.sreads[0],
                                              self.window_len,
                                              self.min_qual)
        nt.assert_is_instance(trimmed_sread, Bio.SeqRecord.SeqRecord)


class TestQualityTrimFastq:
    """Test quality_trim_fastq function in fq module."""

    def setup(self):
        """Create objects needed for testing."""
        pass

    def teardown(self):
        """Teardown.

        Remove files created during the test.
        """

        # all files created during test should start with 'test_'
        file_pattern = r'^test_.*\.fq$'
        for_removal = [f for f in os.listdir(TEST_DIR)
                       if re.search(file_pattern, f)]
        [os.remove(os.path.join(TEST_DIR, f)) for f in for_removal]

    def test_output(self):
        """Test if quality_trim_fastq creates file with results"""
        # run trimming on example file
        fastq_output = os.path.join(TEST_DIR, 'test_trimmed.fq')
        window_len = 3
        min_qual = 20
        ffq.quality_trim_fastq(SAMPLE_FQ1, fastq_output,
                               window_len, min_qual)
        # check if file has been created
        file_exists = os.path.exists(fastq_output)
        nt.assert_true(file_exists)


class TestFilterByLengthFastq:
    """ Test test_filter_by_length_fastq function in fq module."""

    fastq_output = os.path.join(TEST_DIR, 'test_filt.fq')

    def setup(self):
        """Create objects needed for testing."""
        pass

    def teardown(self):
        """Teardown.

        Remove files created during the test.
        """

        # all files created during test should start with 'test_'
        file_pattern = r'^test_.*\.fq$'
        for_removal = [f for f in os.listdir(TEST_DIR)
                       if re.search(file_pattern, f)]
        [os.remove(os.path.join(TEST_DIR, f)) for f in for_removal]

    def test_run_without_min_len_and_max_len(self):
        """Test filter_by_length_fastq without min_len and max_len
        provided.

        Should write to new file with all same reads.
        """

        # if max_len and min_len are not specified all the reads without
        # filtering should be written to output file
        ffq.filter_by_length_fastq(SAMPLE_FQ1, self.fastq_output)

        # check if the file was created
        file_exists = os.path.exists(self.fastq_output)
        nt.assert_true(file_exists)

        # check if it has the same number of reads
        with open(SAMPLE_FQ1, 'r') as fh:
            n_input_reads = len([line for line in fh
                                 if line.startswith('@')])

        with open(self.fastq_output, 'r') as fh:
            n_output_reads = len([line for line in fh
                                 if line.startswith('@')])

        nt.assert_equal(n_input_reads, n_output_reads)

    def test_filtering_with_upper_bound(self):
        """Test filter_by_length_fastq with upper bound
        (shorter or equal to)."""

        max_len = 75
        ffq.filter_by_length_fastq(SAMPLE_FQ1,
                                   self.fastq_output,
                                   max_len=max_len)
        filt_seqs = [seq for seq in SeqIO.parse(self.fastq_output, 'fastq')]
        filt_lengths = [len(seq) for seq in filt_seqs]
        # nb of reads above threshold should be zero
        n_long = len([l for l in filt_lengths if l > max_len])
        nt.assert_equal(n_long, 0)
        # also knowing what's in the fq file there should be 2 seqs
        # shorter than 75
        nt.assert_equal(len(filt_seqs), 2)

    def test_filtering_with_lower_bound(self):
        """Test filter_by_length_fastq with lower bound
        (longer or equal to)."""

        min_len = 70
        ffq.filter_by_length_fastq(SAMPLE_FQ1,
                                   self.fastq_output,
                                   min_len=min_len)
        filt_seqs = [seq for seq in SeqIO.parse(self.fastq_output, 'fastq')]
        filt_lengths = [len(seq) for seq in filt_seqs]
        n_short = len([l for l in filt_lengths if l < min_len])
        nt.assert_equal(n_short, 0)
        # also knowing what's in the fq file there should be 2 seqs
        # longer than 70
        nt.assert_equal(len(filt_seqs), 2)

    def test_filtering_with_range(self):
        """ Test filter_by_length_fastq with a range
        (longer or equal to and shorter or equal to)
        """
        min_len = 20
        max_len = 50
        ffq.filter_by_length_fastq(SAMPLE_FQ1,
                                   self.fastq_output,
                                   min_len=min_len,
                                   max_len=max_len)
        filt_seqs = [seq for seq in SeqIO.parse(self.fastq_output, 'fastq')]
        filt_lengths = [len(seq) for seq in filt_seqs]
        n_out_range = len([l for l in filt_lengths
                          if (l < min_len and l > max_len)])
        nt.assert_equal(n_out_range, 0)
        # also knowing what's in the fq file there should be 1 seq
        # in range
        nt.assert_equal(len(filt_seqs), 1)
