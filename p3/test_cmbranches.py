"""Unit tests (sort of) for all the possible paths
through the P3 logic framework that involve the
sky catalog cross-matching stage
(branches 10, 12, 13, 15, 16, 17, 20, 21)

A single test can be run by calling unittest.main()
and by typing the following at the command prompt:
$ python test_cmbranches.py TestCMBranches.test_sfcm_reprocess

"""
import unittest
import os
from p3 import dbinit, process


class TestCMBranches(unittest.TestCase):
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


    def test_sfcm_new_image(self):
        """Branch 10 - SF + CM on new image"""
        # sf, sa, cm
        stages = (True, False, True)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # No results past SF should be saved to the DB
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        img_rows = self.cur.fetchall()
        sorted_img_rows = sorted(img_rows, key=lambda tup: tup[0])
        self.cur.execute('SELECT COUNT(1) FROM catalog_match')
        matches = self.cur.fetchone()[0]
        result = [sorted_img_rows, matches]
        self.assertEqual(result, [[(1, 2), (2, 2)], 0])
        self.cur.close()


    def test_all_new_image(self):
        """Branch 12 - SF + SA + CM on new image"""
        # sf, sa, cm
        stages = (True, True, True)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Check DB
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        img_rows = self.cur.fetchall()
        sorted_img_rows = sorted(img_rows, key=lambda tup: tup[0])
        self.cur.execute('SELECT COUNT(1) FROM assoc_source WHERE nmatches = 1')
        assoc_matches = self.cur.fetchone()[0] # 21
        self.cur.execute('SELECT COUNT(1) FROM catalog_match')
        catalog_matches = self.cur.fetchone()[0] # 21
        self.cur.execute('SELECT COUNT(1) FROM vlite_unique WHERE detected')
        unique_detections = self.cur.fetchone()[0] # 13
        result = [sorted_img_rows, assoc_matches, catalog_matches,
                  unique_detections]
        self.assertEqual(result, [[(1, 4), (2, 4)], 20, 20, 14])
        self.cur.close()


    def test_sfcm_reprocess(self):
        """Branch 13 - SF + CM with reprocessing"""
        # sf, sa, cm
        # Add image to DB first
        stages = (False, False, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, True, False, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Now do SF + CM
        stages = (True, False, True)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # No results past SF should be saved to the DB
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        img_rows = self.cur.fetchall()
        sorted_img_rows = sorted(img_rows, key=lambda tup: tup[0])
        self.cur.execute('SELECT COUNT(1) FROM catalog_match')
        matches = self.cur.fetchone()[0]
        result = [sorted_img_rows, matches]
        self.assertEqual(result, [[(1, 2), (2, 2)], 0])
        self.cur.close()


    def test_all_reprocess(self):
        """Branch 15 - SF + SA + CM with reprocessing"""
        # sf, sa, cm
        # Add image to DB first
        stages = (False, False, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, True, False, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Now run all stages
        stages = (True, True, True)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Check DB
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        img_rows = self.cur.fetchall()
        sorted_img_rows = sorted(img_rows, key=lambda tup: tup[0])
        self.cur.execute('SELECT COUNT(1) FROM assoc_source WHERE nmatches = 1')
        assoc_matches = self.cur.fetchone()[0] # 20
        self.cur.execute('SELECT COUNT(1) FROM catalog_match')
        catalog_matches = self.cur.fetchone()[0] # 20
        self.cur.execute('SELECT COUNT(1) FROM vlite_unique WHERE detected')
        unique_detections = self.cur.fetchone()[0] # 14
        result = [sorted_img_rows, assoc_matches, catalog_matches,
                  unique_detections]
        self.assertEqual(result, [[(1, 4), (2, 4)], 20, 20, 14])
        self.cur.close()


    def test_cmonly_fail(self):
        """Branch 16 - CM only, but stage < 2"""
        # sf, sa, cm
        # Populate DB with sources first
        stages = (True, False, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Now try to run CM only
        stages = (False, False, True)
        # Should fail due to stage < 2
        self.assertIsNone(process(self.conn, stages, opts, self.dirs,
                                  self.catalogs, self.params, self.res_tol))


    def test_cmonly(self):
        """Branch 17a - new CM only"""
        # sf, sa, cm
        # Run SF + SA first
        stages = (True, True, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Now run CM only
        stages = (False, False, True)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Check DB - should be same as branch 12/15 (SF + SA + CM)
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        img_rows = self.cur.fetchall()
        sorted_img_rows = sorted(img_rows, key=lambda tup: tup[0])
        self.cur.execute('SELECT COUNT(1) FROM assoc_source WHERE nmatches = 1')
        assoc_matches = self.cur.fetchone()[0] # 20
        self.cur.execute('SELECT COUNT(1) FROM catalog_match')
        catalog_matches = self.cur.fetchone()[0] # 20
        self.cur.execute('SELECT COUNT(1) FROM vlite_unique WHERE detected')
        unique_detections = self.cur.fetchone()[0] # 14
        result = [sorted_img_rows, assoc_matches, catalog_matches,
                  unique_detections]
        self.assertEqual(result, [[(1, 4), (2, 4)], 20, 20, 14])
        self.cur.close()


    def test_cmredo(self):
        """Branch 17b - redo existing CM only"""
        # sf, sa, cm
        # Run SF + SA + CM first
        stages = (True, True, True)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Now redo CM
        stages = (False, False, True)
        opts = (True, False, False, False, True, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Check DB - should be same as branch 12/15 (SF + SA + CM)
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        img_rows = self.cur.fetchall()
        sorted_img_rows = sorted(img_rows, key=lambda tup: tup[0])
        self.cur.execute('SELECT COUNT(1) FROM assoc_source WHERE nmatches = 1')
        assoc_matches = self.cur.fetchone()[0] # 20
        self.cur.execute('SELECT COUNT(1) FROM catalog_match')
        catalog_matches = self.cur.fetchone()[0] # 20
        self.cur.execute('SELECT COUNT(1) FROM vlite_unique WHERE detected')
        unique_detections = self.cur.fetchone()[0] # 14
        result = [sorted_img_rows, assoc_matches, catalog_matches,
                  unique_detections]
        self.assertEqual(result, [[(1, 4), (2, 4)], 20, 20, 14])
        self.cur.close()


    def test_cmupdate(self):
        """Branch 17c - add new catalog to existing CM results"""
        # sf, sa, cm
        # Run SF + SA + CM with NVSS catalog first
        stages = (True, True, True)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Now run CM again adding the FIRST catalog
        stages = (False, False, True)
        opts = (True, False, False, False, False, True)
        self.catalogs = ['FIRST']
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Check DB - should have results for both NVSS & FIRST
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        img_rows = self.cur.fetchall()
        sorted_img_rows = sorted(img_rows, key=lambda tup: tup[0])
        self.cur.execute('SELECT COUNT(1) FROM assoc_source WHERE nmatches > 0')
        assoc_matches = self.cur.fetchone()[0] # 23
        self.cur.execute('''SELECT COUNT(1) FROM catalog_match
            WHERE catalog_id = 7''')
        nvss_matches = self.cur.fetchone()[0] # 20
        self.cur.execute('''SELECT COUNT(1) FROM catalog_match
            WHERE catalog_id = 2''')
        first_matches = self.cur.fetchone()[0] # 21
        self.cur.execute('SELECT COUNT(1) FROM vlite_unique WHERE detected')
        unique_detections = self.cur.fetchone()[0] # 10
        result = [sorted_img_rows, assoc_matches, nvss_matches,
                  first_matches, unique_detections]
        self.assertEqual(result, [[(1, 4), (2, 4)], 23, 20, 21, 10])
        self.cur.close()


    def test_sacm(self):
        """Branch 20 - SA + CM"""
        # sf, sa, cm
        # Run SF first
        stages = (True, False, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Now run SA + CM
        stages = (False, True, True)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Check DB - should be same as branch 12/15 (SF + SA + CM)
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        img_rows = self.cur.fetchall()
        sorted_img_rows = sorted(img_rows, key=lambda tup: tup[0])
        self.cur.execute('SELECT COUNT(1) FROM assoc_source WHERE nmatches = 1')
        assoc_matches = self.cur.fetchone()[0] # 20
        self.cur.execute('SELECT COUNT(1) FROM catalog_match')
        catalog_matches = self.cur.fetchone()[0] # 20
        self.cur.execute('SELECT COUNT(1) FROM vlite_unique WHERE detected')
        unique_detections = self.cur.fetchone()[0] # 14
        result = [sorted_img_rows, assoc_matches, catalog_matches,
                  unique_detections]
        self.assertEqual(result, [[(1, 4), (2, 4)], 20, 20, 14])
        self.cur.close()


    def test_sapass_cmonly(self):
        """Branch 21 - SA already done, CM only"""
        # sf, sa, cm
        # Run SF + SA first
        stages = (True, True, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Now try to run SA + CM
        stages = (False, True, True)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Check DB - should be same results as branch 17
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        img_rows = self.cur.fetchall()
        sorted_img_rows = sorted(img_rows, key=lambda tup: tup[0])
        self.cur.execute('SELECT COUNT(1) FROM assoc_source WHERE nmatches = 1')
        assoc_matches = self.cur.fetchone()[0] # 20
        self.cur.execute('SELECT COUNT(1) FROM catalog_match')
        catalog_matches = self.cur.fetchone()[0] # 20
        self.cur.execute('SELECT COUNT(1) FROM vlite_unique WHERE detected')
        unique_detections = self.cur.fetchone()[0] # 14
        result = [sorted_img_rows, assoc_matches, catalog_matches,
                  unique_detections]
        self.assertEqual(result, [[(1, 4), (2, 4)], 20, 20, 14])
        self.cur.close()


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCMBranches)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
