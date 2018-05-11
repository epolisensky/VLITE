"""Integration tests for all the possible paths
through the VDP logic framework that involve either
no stages or source finding only (branches 2-9). 

A single test can be run by calling unittest.main()
and by typing the following at the command prompt:
>> python test_sfbranches.py TestSFBranches.test_no_stages_add_new_image

"""
import unittest
import os
import sys
sys.path.append('../')
import vdp


class TestSFBranches(unittest.TestCase):
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


    def test_no_stages_add_new_image(self):
        """Branch 2 - add new image to DB"""
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        # Pipeline should stop after adding image to DB
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 2)


    def test_no_stages_no_reprocess(self):
        """Branch 3 - image already in DB, no reprocess"""
        # Need to add images to DB first
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Now try to add again; nothing should happen
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 3)


    def test_no_stages_reprocess(self):
        """Branch 4 - image already in DB, reprocess"""
        # Need to add images to DB first
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : True,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Now re-add
        opts['save to database'] = False
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 4)


    def test_no_sf_new_image(self):
        """Branch 5 - trying to sa/cm before sf, new image"""
        stages = {'source finding' : False, 'source association' : True,
                  'catalog matching' : True}
        opts = {'save to database' : False, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 5)


    def test_sfonly_new_image(self):
        """Branch 6 - source findng on new image"""
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : False, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        # Pipeline should stop after source finding
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 6)
        

    def test_no_sf_stage1(self):
        """Branch 7 - trying to sa/cm before sf, existing image"""
        # Add images to DB first
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Now try to sa/cm before sf
        stages = {'source finding' : False, 'source association' : True,
                  'catalog matching' : True}
        opts['save to database'] = False
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 7)


    def test_sfonly_reprocess(self):
        """Branch 8 - redo source finding on existing image"""
        # Add sources to database
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : True,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Now redo source finding
        opts['save to database'] = False
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 8)


    def test_sf_no_reprocess(self):
        """Branch 9 - SF, but image already in DB & no reprocess"""
        # Need to add images & sources to DB first
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        # Now try to reprocess; nothing should happen
        opts['save to database'] = False
        vdp.process(self.conn, stages, opts, self.dirs, self.files,
                    self.catalogs, self.sfparams, self.qaparams)
        self.assertEqual(vdp.branch, 9)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSFBranches)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
