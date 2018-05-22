"""Integration tests for all the possible paths
through the VDP logic framework that involve the
sky catalog cross-matching stage
(branches 10, 12, 13, 15, 16, 17, 20, 21)

A single test can be run by calling unittest.main()
and by typing the following at the command prompt:
$ python test_cmbranches.py TestCMBranches.test_sfcm_new_image

"""
import unittest
import os
import sys
sys.path.append('../')
import vdp


class TestCMBranches(unittest.TestCase):
    """Tests for correct logic flow based on stages & options."""

    def setUp(self):
        """Define variables to be used by all tests and
        connect the database, overwriting existing tables
        every time."""
        self.dirs = ['/data3/vpipe/processed/2015-03/31/Images/',
                     '/nfsshare/vpipe/processed/2017-07/25/Images/']
        self.files = [['13B-266.deepfield.IPln1.fits'],
                      ['1.5GHz.J0205-206.IPln1.fits']]
        self.catalogs = ['NVSS']
        self.sfparams = {'mode' : 'default', 'thresh' : 'hard', 'scale' : 1.0}
        self.qaparams = {'min nvis' : 1000.,
                         'max sensitivity metric' : 3000.,
                         'max beam axis ratio' : 4.,
                         'max source count metric' : 10.}
        self.conn = vdp.dbinit('branchtest', 'vpipe', True, self.qaparams, True)


    def tearDown(self):
        """Disconnect the database."""
        self.conn.close()


    def test_sfcm_new_image(self):
        """Branch 10 - SF + CM on new image"""
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : True}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 10)


    def test_all_new_image(self):
        """Branch 12 - SF + SA + CM on new image"""
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : True}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 12)
        

    def test_sfcm_reprocess(self):
        """Branch 13 - SF + CM with reprocessing"""
        # Add images & sources to DB first
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : True,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Re-do SF & add in CM
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : True}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 13)


    def test_all_reprocess(self):
        """Branch 15 - SF + SA + CM with reprocessing"""
        # Add image & sources to DB first
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : True,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Re-do SF & add SA + CM
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : True}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 15)


    def test_cmonly_fail(self):
        """Branch 16 - CM only, but stage < 2"""
        # Add images & sources to database
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Now try to run CM only
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : True}
        opts['save to database'] = False
        # Should fail due to stage < 2
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 16)


    def test_cmonly(self):
        """Branch 17.1 - new CM only"""        
        # Run SF + SA first
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}        
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Now run CM only
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : True}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 17.1)


    def test_cmredo(self):
        """Branch 17.2 - redo existing CM only"""
        # Run SF + SA + CM first
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : True}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Now redo CM
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : True}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : True, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 17.2)


    def test_cmupdate(self):
        """Branch 17.3 - add new catalog to existing CM results"""
        # Run SF + SA + CM with NVSS catalog only
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : True}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Now run CM again adding the WENSS catalog
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : True}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : True}
        self.catalogs = ['WENSS']
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 17.3)


    def test_sacm(self):
        """Branch 20 - SA + CM"""
        # Run SF first
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Now run SA + CM
        stages = {'source finding' : False, 'source association' : True,
                  'catalog matching' : True}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 20)


    def test_sapass_cmonly(self):
        """Branch 21 - SA already done, CM only"""
        # Run SF + SA first
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Now try to run SA + CM
        stages = {'source finding' : False, 'source association' : True,
                  'catalog matching' : True}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 21)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCMBranches)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
