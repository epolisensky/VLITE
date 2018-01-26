"""Integration tests for all the possible paths
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
        self.conn = dbinit('branchtest', 'erichards', True, True)


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
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
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
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : True}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
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


    def test_sfcm_reprocess(self):
        """Branch 13 - SF + CM with reprocessing"""
        # Add image to DB first
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : True,
                'redo match' : False, 'update match' : False}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
        # Now do SF + CM
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : True}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
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
        # Add image to DB first
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : True,
                'redo match' : False, 'update match' : False}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
        # Now run all stages
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : True}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
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
        # Populate DB with sources first
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
        # Now try to run CM only
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : True}
        # Should fail due to stage < 2
        self.assertIsNone(process(self.conn, stages, opts, self.dirs,
                                  self.catalogs, self.params))


    def test_cmonly(self):
        """Branch 17a - new CM only"""        
        # Run SF + SA first
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}        
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
        # Now run CM only
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : True}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
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
        # Run SF + SA + CM first
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : True}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
        # Now redo CM
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : True}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : True, 'update match' : False}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
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
        # Run SF + SA + CM with NVSS catalog first
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : True}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
        # Now run CM again adding the FIRST catalog
        stages = {'source finding' : False, 'source association' : False,
                  'catalog matching' : True}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : True}
        self.catalogs = ['NVSS', 'FIRST']
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
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
        first_matches = self.cur.fetchone()[0] # 22
        self.cur.execute('SELECT COUNT(1) FROM vlite_unique WHERE detected')
        unique_detections = self.cur.fetchone()[0] # 10
        result = [sorted_img_rows, assoc_matches, nvss_matches,
                  first_matches, unique_detections]
        self.assertEqual(result, [[(1, 4), (2, 4)], 23, 20, 22, 10])
        self.cur.close()


    def test_sacm(self):
        """Branch 20 - SA + CM"""
        # Run SF first
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
        # Now run SA + CM
        stages = {'source finding' : False, 'source association' : True,
                  'catalog matching' : True}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
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
        # Run SF + SA first
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
        # Now try to run SA + CM
        stages = {'source finding' : False, 'source association' : True,
                  'catalog matching' : True}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
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
