"""Tests for farad.qc module"""
import os
import re
import sys
import nose.tools as nt

# import tested package
sys.path.append('../')
import farad.qc as fqc

# test data

class TestFaradQCDefault():
    """Basic farad.qc testing class.

    Setups and treardown for farad.qc testing.
    """
    def setup(self):
        """Create objects needed for testing"""
        self.test_dir = os.path.split(os.path.abspath(__file__))[0]
        self.sample_fq1 = os.path.join(self.test_dir, 'sample1.fq')

    def teardown(self):
        """Teardown everything created during tests

        Remove created files.
        """
        file_pattern = r'^test_.*\.fq$'
        for_removal = [f for f in os.listdir(self.test_dir)
                       if re.search(file_pattern, f)]
        [os.remove(os.path.join(self.test_dir, f)) for f in for_removal]

class TestPlotQUalitiesAlongRead(TestFaradQCDefault):
    """Test plot_qualities_along_read function in farad.qc module"""
    def test_exception_on_wrong_input(self):
        """Should throw an exception on wrong input file"""





# test for exception if input is not valid fastq file
#
# test for exception if output path is not valid
#
# test if output file is created in expected path
#

