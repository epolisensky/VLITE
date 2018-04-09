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
from vdp import dbinit, process


class TestSFBranches(unittest.TestCase):
    """Tests for correct logic flow based on stages & options."""

    def setUp(self):
        """Define variables to be used by all tests and
        connect the database, overwriting existing tables
        every time."""
        self.dirs = ['/home/erichards/work/data/test/2540-06/04/Images/']
        self.files = [[]]
        self.catalogs = ['NVSS']
        self.sfparams = {'mode' : 'default', 'thresh' : 'hard', 'scale' : 0.5}
        self.qaparams = {'min time on source (s)' : 60.,
                         'max noise (mJy/beam)' : 1000.,
                         'max beam axis ratio' : 4.,
                         'min problem source separation (deg)' : 20.,
                         'max source metric' : 10.}
        self.conn = dbinit('branchtest', 'erichards', True, self.qaparams, True)


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
        # Uses fresh DB, no need to add images first
        # Pipeline should stop after adding image to DB
        process(self.conn, stages, opts, self.dirs, self.files,
                self.catalogs, self.sfparams, self.qaparams)
        # Check DB - image table should have 2 rows
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        rows = self.cur.fetchall()
        sorted_rows = sorted(rows, key=lambda tup: tup[0])
        self.assertEqual(sorted_rows, [(1, 1), (2, 1)])
        self.cur.close()


    def test_no_stages_no_reprocess(self):
        """Branch 3 - image already in DB, no reprocess"""
        # Need to add images to DB first
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        # Add images to DB
        process(self.conn, stages, opts, self.dirs, self.files,
                self.catalogs, self.sfparams, self.qaparams)
        # Now try to reprocess; nothing should happen
        self.assertIsNone(process(self.conn, stages, opts, self.dirs,
                                  self.files, self.catalogs, self.sfparams,
                                  self.qaparams))


    def test_no_stages_reprocess(self):
        """Branch 4 - image already in DB, reprocess"""
        # Need to add images to DB first
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        # Add images to DB
        process(self.conn, stages, opts, self.dirs, self.files,
                self.catalogs, self.sfparams, self.qaparams)
        # Now reprocess
        process(self.conn, stages, opts, self.dirs, self.files,
                self.catalogs, self.sfparams, self.qaparams)
        # Check DB - image table should have 2 rows
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        rows = self.cur.fetchall()
        sorted_rows = sorted(rows, key=lambda tup: tup[0])
        self.assertEqual(sorted_rows, [(1, 1), (2, 1)])
        self.cur.close()


    def test_no_sf_new_image(self):
        """Branch 5 - trying to sa/cm before sf, new image"""
        stages = {'source finding' : False, 'source association' : True,
                  'catalog matching' : True}
        opts = {'save to database' : False, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        self.assertIsNone(process(self.conn, stages, opts, self.dirs,
                                  self.files, self.catalogs, self.sfparams,
                                  self.qaparams))


    def test_sfonly_new_image(self):
        """Branch 6 - source findng on new image"""
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        # Pipeline should stop after source finding
        process(self.conn, stages, opts, self.dirs, self.files,
                self.catalogs, self.sfparams, self.qaparams)
        # Check DB
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        img_rows = self.cur.fetchall()
        sorted_img_rows = sorted(img_rows, key=lambda tup: tup[0])
        self.cur.execute('''SELECT COUNT(1) FROM detected_island
            WHERE image_id = 1''')
        nisls1 = self.cur.fetchone()[0]
        self.cur.execute('''SELECT COUNT(1) FROM detected_island
            WHERE image_id = 2''')
        nisls2 = self.cur.fetchone()[0]
        self.cur.execute('''SELECT COUNT(1) FROM detected_source
            WHERE image_id = 1''')
        nsrcs1 = self.cur.fetchone()[0]
        self.cur.execute('''SELECT COUNT(1) FROM detected_source
            WHERE image_id = 2''')
        nsrcs2 = self.cur.fetchone()[0]
        result = [sorted_img_rows, nisls1, nisls2, nsrcs1, nsrcs2]
        self.assertEqual(result, [[(1, 2), (2, 2)], 28, 15, 29, 15])
        self.cur.close()


    def test_no_sf_stage1(self):
        """Branch 7 - trying to sa/cm before sf, existing image"""
        # Add images to DB first - sf, sa, cm
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        process(self.conn, stages, opts, self.dirs, self.files,
                self.catalogs, self.sfparams, self.qaparams)
        # Now try to sa/cm before sf
        stages = {'source finding' : False, 'source association' : True,
                  'catalog matching' : True}
        self.assertIsNone(process(self.conn, stages, opts, self.dirs,
                                  self.files, self.catalogs, self.sfparams,
                                  self.qaparams))


    def test_sfonly_reprocess(self):
        """Branch 8 - redo source finding on existing image"""
        # Process through sf so I can make sure sources get deleted
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : True,
                'redo match' : False, 'update match' : False}
        process(self.conn, stages, opts, self.dirs, self.files,
                self.catalogs, self.sfparams, self.qaparams)
        # Use different scale so I know the originals were removed
        self.sfparams = {'mode' : 'default', 'thresh' : 'hard', 'scale' : 0.3}
        # Now redo source finding
        process(self.conn, stages, opts, self.dirs, self.files,
                self.catalogs, self.sfparams, self.qaparams)
        # Check DB
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        img_rows = self.cur.fetchall()
        sorted_img_rows = sorted(img_rows, key=lambda tup: tup[0])
        self.cur.execute('''SELECT COUNT(1) FROM detected_island
            WHERE image_id = 1''')
        nisls1 = self.cur.fetchone()[0]
        self.cur.execute('''SELECT COUNT(1) FROM detected_island
            WHERE image_id = 2''')
        nisls2 = self.cur.fetchone()[0]
        self.cur.execute('''SELECT COUNT(1) FROM detected_source
            WHERE image_id = 1''')
        nsrcs1 = self.cur.fetchone()[0]
        self.cur.execute('''SELECT COUNT(1) FROM detected_source
            WHERE image_id = 2''')
        nsrcs2 = self.cur.fetchone()[0]
        result = [sorted_img_rows, nisls1, nisls2, nsrcs1, nsrcs2]
        self.assertEqual(result, [[(1, 2), (2, 2)], 15, 8, 16, 8])
        self.cur.close()


    def test_sf_no_reprocess(self):
        """Branch 9 - sf, but image already in DB & no reprocess"""
        # Need to add images to DB first
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        # Add images to DB
        process(self.conn, stages, opts, self.dirs, self.files,
                self.catalogs, self.sfparams, self.qaparams)
        # Now try to reprocess; nothing should happen
        self.assertIsNone(process(self.conn, stages, opts, self.dirs,
                                  self.files, self.catalogs, self.sfparams,
                                  self.qaparams))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSFBranches)
    #unittest.TextTestRunner(verbosity=2).run(suite)
    unittest.main()
