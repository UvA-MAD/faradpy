"""Tests for farad."""
import os
import sys
import Bio
from Bio import SeqIO
from nose.tools import assert_equal, assert_is_instance

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

    def test_qualit_trim_read(self):
        """Test quality trimming with sliding window"""
        # load test fastq
        test_script_path = os.path.split(os.path.abspath(__file__))[0]
        fq_path = os.path.join(test_script_path, 'sample1.fq')
        sreads = [r for r in SeqIO.parse(fq_path, 'fastq')]

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
