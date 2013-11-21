"""Tests for farad."""
import os
import sys
import Bio
from Bio import SeqIO
from nose.tools import assert_equal, assert_is_instance, assert_true

# import tested package
sys.path.append('../')
# import farad.utils as futils
import farad.fq as ffq


class TestFaradFq:

    def setup(self):
        """Create objects needed for testing."""
        pass

    def teardown(self):
        """Teardown."""
        pass

    # load test fastq
    test_script_path = os.path.split(os.path.abspath(__file__))[0]
    fq_path = os.path.join(test_script_path, 'sample1.fq')

    def test_qualit_trim_read(self):
        """Test quality trimming with sliding window"""
        sreads = [r for r in SeqIO.parse(self.fq_path, 'fastq')]

        # test trimming
        window_len = 3
        min_qual = 21
        # first read of len 82 expected last 2 bases to be trimmed
        # with min_qual=21
        trimmed_sread1 = ffq.quality_trim_read(sreads[0], window_len, min_qual)
        assert_equal(len(trimmed_sread1), 80)
        # second read of len 29 expected to get one base trimmed
        trimmed_sread2 = ffq.quality_trim_read(sreads[1], window_len, min_qual)
        assert_equal(len(trimmed_sread2), 28)

        # make sure that the return object is Bio.SeqRecord.SeqRecord
        assert_is_instance(trimmed_sread2, Bio.SeqRecord.SeqRecord)

    def test_quality_trim_fastq(self):
        """Test quality trimming of all reads in fastq file."""
        # run trimming on example file
        fastq_output = os.path.join(self.test_script_path, 'trimmed_reads.fq')
        window_len = 3
        min_qual = 20
        ffq.quality_trim_fastq(self.fq_path, fastq_output,
                               window_len, min_qual)
        # check if file has been created
        file_exists = os.path.exists(fastq_output)
        assert_true(file_exists)

        # clean up
        os.remove(fastq_output)
