import unittest
import p3
from errors import ConfigError


class TestConfigParser(unittest.TestCase):
    """Tests for cfgparse method in p3.py."""

    def test_no_stages_no_save(self):
        """Raise ConfigError if running no stages & not saving."""
        self.assertRaises(ConfigError, p3.cfgparse,
                          'test/config_files/branch1.yaml')


    def test_bad_monthdir(self):
        """Raise ConfigError if path to monthly directory does not exist."""
        self.assertRaises(ConfigError, p3.cfgparse,
                          'test/config_files/baddirs.yaml')


    def test_bad_days(self):
        """Raise ConfigError if setup: day: input is wrong."""
        self.assertRaises(ConfigError, p3.cfgparse,
                          'test/config_files/baddays.yaml')


    def test_skycatalogDB(self):
        """Raise ConfigError if path to SkyCatalog DB does not exist."""
        self.assertRaises(ConfigError, p3.cfgparse,
                          'test/config_files/badskyDB.yaml')


    def test_skycatalogs(self):
        """Raise ConfigError if sky catalogs are incorrect."""
        self.assertRaises(ConfigError, p3.cfgparse,
                          'test/config_files/badskycats.yaml')


#class TestImageInit(unittest.TestCase):
    """Tests for iminit method in p3.py."""


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestConfigParser)
    unittest.TextTestRunner(verbosity=2).run(suite)
