import unittest

class InitializationTests(unittest.TestCase):

    def test_initialization(self):
        """
        Check that the test suite runs.
        """
        self.assertEqual(2+2, 4)

    def test_import(self):
        """
        Check that the modules can be imported
        """
        try:
            from chemical_curation import curate
        except ImportError:
            self.fail("Unable to import module curate from chemical_curation")
