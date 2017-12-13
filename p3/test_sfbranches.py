"""Unit tests (sort of) for all the possible paths
through the P3 logic framework that involve either
no stages or source finding only (branches 2-9). 

A single test can be run by calling unittest.main()
and by typing the following at the command prompt:
$ python test_sfbranches.py TestSFBranches.test_no_stages_add_new_image

"""
import unittest
import os
from p3 import dbinit, process


class TestSFBranches(unittest.TestCase):
    """Tests for correct logic flow based on stages & options."""

    def setUp(self):
        """Define variables to be used by all tests and
        connect the database, overwriting existing tables
        every time."""
        self.dirs = ['/home/erichards/work/p3/test/2540-06/B/Images/']
        self.catalogs = ['NVSS']
        self.params = {'mode' : 'default', 'thresh' : 'hard', 'scale' : 0.5}
        self.res_tol = 5.0
        self.conn = dbinit('branchtest', 'erichards', True, True)


    def tearDown(self):
        """Disconnect the database."""
        self.conn.close()


    def test_no_stages_add_new_image(self):
        """Branch 2 - add new image to DB"""
        # sf, sa, cm
        stages = (False, False, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        # Uses fresh DB, no need to add images first
        # Pipeline should stop after adding image to DB
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
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
        # sf, sa, cm
        stages = (False, False, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        # Add images to DB
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Now try to reprocess; nothing should happen
        self.assertIsNone(process(self.conn, stages, opts, self.dirs,
                                  self.catalogs, self.params, self.res_tol))


    def test_no_stages_reprocess(self):
        """Branch 4 - image already in DB, reprocess"""
        # Need to add images to DB first
        # sf, sa, cm
        stages = (False, False, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, True, False, False)
        # Add images to DB
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Now reprocess
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Check DB - image table should have 2 rows
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        rows = self.cur.fetchall()
        sorted_rows = sorted(rows, key=lambda tup: tup[0])
        self.assertEqual(sorted_rows, [(1, 1), (2, 1)])
        self.cur.close()


    def test_no_sf_new_image(self):
        """Branch 5 - trying to sa/cm before sf, new image"""
        # sf, sa, cm
        stages = (False, True, True)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (False, False, False, False, False, False)
        self.assertIsNone(process(self.conn, stages, opts, self.dirs,
                                  self.catalogs, self.params, self.res_tol))


    def test_sfonly_new_image(self):
        """Branch 6 - source findng on new image"""
        # sf, sa, cm
        stages = (True, False, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        # Pipeline should stop after source finding
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
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
        stages = (False, False, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Now try to sa/cm before sf
        stages = (False, True, True)
        self.assertIsNone(process(self.conn, stages, opts, self.dirs,
                                  self.catalogs, self.params, self.res_tol))


    def test_sfonly_reprocess(self):
        """Branch 8 - redo source finding on existing image"""
        # Process through sf so I can make sure sources get deleted
        stages = (True, False, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, True, False, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Use different scale so I know the originals were removed
        self.params = {'mode' : 'default', 'thresh' : 'hard', 'scale' : 0.3}
        # Now redo source finding
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
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
        # sf, sa, cm
        stages = (True, False, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        # Add images to DB
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Now try to reprocess; nothing should happen
        self.assertIsNone(process(self.conn, stages, opts, self.dirs,
                                  self.catalogs, self.params, self.res_tol))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSFBranches)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
