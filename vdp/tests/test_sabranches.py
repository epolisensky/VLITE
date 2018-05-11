"""Integration tests for all the possible paths
through the VDP logic framework that involve stages
up through source association (branches 11, 14, 18, 19). 

A single test can be run by calling unittest.main()
and by typing the following at the command prompt:
$ python test_sabranches.py TestSABranches.test_sfsa_new_image

"""
import unittest
import os
import sys
sys.path.append('../')
import vdp


class TestSABranches(unittest.TestCase):
    """Tests for correct logic flow based on stages & options."""

    def setUp(self):
        """Define variables to be used by all tests and
        connect the database, overwriting existing tables
        every time."""
        self.dirs = ['/home/erichards/work/data/test/2540-06/03/Images/']
        self.files = [[]]
        self.catalogs = ['NVSS']
        self.sfparams = {'mode' : 'default', 'thresh' : 'hard', 'scale' : 1.0}
        self.qaparams = {'min nvis' : 1000.,
                         'max sensitivity metric' : 3000.,
                         'max beam axis ratio' : 4.,
                         'max source count metric' : 10.}
        self.conn = vdp.dbinit('branchtest', 'erichards', True,
                               self.qaparams, True)


    def tearDown(self):
        """Disconnect the database."""
        self.conn.close()


    def test_sfsa_new_image(self):
        """Branch 11 - SF + SA on new image"""
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : False}
        opts = {'save to database' : False, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        # Pipeline should stop after source association
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 11)


    def test_sfsa_reprocess(self):
        """Branch 14 - SF + SA reprocess"""
        # Add sources to database
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : True,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Reprocess with SA
        opts['save to database'] = False
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 14)
        

    def test_saonly_already_done(self):
        """Branch 18 - SA already done, so do nothing"""
        # Process through SA first
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : True,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Now try SA only
        stages = {'source finding' : False, 'source association' : True,
                  'catalog matching' : False}
        opts['save to database'] = False
        # Code should exit
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 18)


    def test_saonly(self):
        """Branch 19 - SA only"""
        # Add sources to database
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Now do SA only
        stages = {'source finding' : False, 'source association' : True,
                  'catalog matching' : False}
        opts['save to database'] = False
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 19)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSABranches)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
