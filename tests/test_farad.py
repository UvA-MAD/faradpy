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


class TestFaradFq:

    # load test fastq
    test_dir = os.path.split(os.path.abspath(__file__))[0]
    sample_fq = os.path.join(test_dir, 'sample1.fq')

    def setup(self):
        """Create objects needed for testing."""
        pass

    def teardown(self):
        """Teardown."""

        # remove files created during the test
        # all files created during test should start with 'test_'
        file_pattern = r'^test_.*\.fq$'
        for_removal = [f for f in os.listdir(self.test_dir)
                       if re.search(file_pattern, f)]

        [os.remove(os.path.join(self.test_dir, f)) for f in for_removal]

    def test_qualit_trim_read(self):
        """Test quality trimming with sliding window"""
        sreads = [r for r in SeqIO.parse(self.sample_fq, 'fastq')]

        # trimming options
        window_len = 3
        min_qual = 21
        # first read of len 82 expected last 2 bases to be trimmed
        # with min_qual=21
        trimmed_sread1 = ffq.quality_trim_read(sreads[0], window_len, min_qual)
        nt.assert_equal(len(trimmed_sread1), 80)
        # second read of len 29 expected to get one base trimmed
        trimmed_sread2 = ffq.quality_trim_read(sreads[1], window_len, min_qual)
        nt.assert_equal(len(trimmed_sread2), 28)

        # make sure that the return object is Bio.SeqRecord.SeqRecord
        nt.assert_is_instance(trimmed_sread2, Bio.SeqRecord.SeqRecord)

    def test_quality_trim_fastq(self):
        """Test quality trimming of all reads in fastq file."""
        # run trimming on example file
        fastq_output = os.path.join(self.test_dir, 'test_trimmed.fq')
        window_len = 3
        min_qual = 20
        ffq.quality_trim_fastq(self.sample_fq, fastq_output,
                               window_len, min_qual)
        # check if file has been created
        file_exists = os.path.exists(fastq_output)
        nt.assert_true(file_exists)

    def test_filter_by_length_fastq(self):
        """Test filtering reads by their length"""

        # if max_len and min_len are not specified all the reads without
        # filtering should be written to output file
        fastq_output = os.path.join(self.test_dir, 'test_filt.fq')
        ffq.filter_by_length_fastq(self.sample_fq, fastq_output)

        # check if the file was created
        file_exists = os.path.exists(fastq_output)
        nt.assert_true(file_exists)

        # check if it has the same number of reads
        n_input_reads = len([line for line
                             in open(self.sample_fq, 'r').readlines()
                             if line.startswith('@')])

        n_output_reads = len([line for line
                             in open(fastq_output, 'r').readlines()
                             if line.startswith('@')])
        nt.assert_equal(n_input_reads, n_output_reads)
