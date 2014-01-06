"""Tests for farad.qc module"""
import os
import sys
import nose.tools as nt

# import tested package
sys.path.append('../')
import farad.qc as fqc


class FaradQCPlottingTest():
    """Basic farad.qc testing class.

    Setups and treardown for farad.qc testing.
    """
    def setup(self):
        """Create objects needed for testing"""
        self.test_dir = os.path.split(os.path.abspath(__file__))[0]
        self.sample_fq1 = os.path.join(self.test_dir, 'sample1.fq')
        self.invalid_input_fq = 'nonexistent.fq'
        self.invalid_output_png = os.path.join(self.test_dir,
                                               'not_there',
                                               'output.png')
        self.output_plot = os.path.join(self.test_dir,
                                        'test_plot.png')

    def teardown(self):
        """Teardown everything created during tests

        Remove created files.
        """
        output_prefix = 'test_'
        output_suffixes = ('.fq', 'png')
        for_removal = [f for f in os.listdir(self.test_dir)
                       if f.startswith(output_prefix)
                       and f.endswith(output_suffixes)]
        [os.remove(os.path.join(self.test_dir, f)) for f in for_removal]

    def test_exception_on_wrong_input(self):
        """Exception should be trown in case of invalid (nonexisting) input."""
        nt.assert_raises(FileNotFoundError,
                         self.tested_function,
                         self.invalid_input_fq,
                         self.output_plot)
        print('in inherited method')

    def test_exception_on_wrong_output_path(self):
        """Should trow an exception if path for output file does not exist."""

        nt.assert_raises(FileNotFoundError,
                         self.tested_function,
                         self.sample_fq1,
                         self.invalid_output_png)

    def test_plot_saved(self):
        """If all goes well a plot file should be saved."""

        fqc.plot_qualities_along_read(self.sample_fq1, self.output_plot)
        file_exists = os.path.exists(self.output_plot)
        nt.assert_true(file_exists)


class TestPlotQUalitiesAlongRead(FaradQCPlottingTest):
    """Test plot_qualities_along_read function in farad.qc module"""

    def __init__(self):
        # each class tests different functions but usually the same way
        self.tested_function = fqc.plot_qualities_along_read


class TestPlotMeanQualityDistribution(FaradQCPlottingTest):
    """Test plot_qualities_along_read function in farad.qc module"""

    def __init__(self):
        # each class tests different functions but usually the same way
        self.tested_function = fqc.plot_mean_quality_distribution
