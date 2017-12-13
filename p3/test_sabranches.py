"""Unit tests (sort of) for all the possible paths
through the P3 logic framework that involve stages
up through source association (branches 11, 14, 18, 19). 

A single test can be run by calling unittest.main()
and by typing the following at the command prompt:
$ python test_sabranches.py TestSABranches.test_sfsa_new_image

"""
import unittest
import os
from p3 import dbinit, process


class TestSABranches(unittest.TestCase):
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


    def test_sfsa_new_image(self):
        """Branch 11 - SF + SA on new image"""
        # sf, sa, cm
        stages = (True, True, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        # Pipeline should stop after source association
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Check DB
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        img_rows = self.cur.fetchall()
        sorted_img_rows = sorted(img_rows, key=lambda tup: tup[0])
        self.cur.execute('''SELECT COUNT(1) FROM assoc_source
            WHERE ndetect = 1''')
        single_detections = self.cur.fetchone()[0]
        self.cur.execute('''SELECT COUNT(1) FROM assoc_source
            WHERE ndetect = 2''')
        associations = self.cur.fetchone()[0]
        result = [sorted_img_rows, single_detections, associations]
        self.assertEqual(result, [[(1, 3), (2, 3)], 24, 10])
        self.cur.close()


    def test_sfsa_reprocess(self):
        """Branch 14 - SF + SA reprocess"""
        # sf, sa, cm
        stages = (True, True, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, True, False, False)
        # Process once
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Use different scale so results are different
        self.params = {'mode' : 'default', 'thresh' : 'hard', 'scale' : 0.3}
        # Process twice
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Check DB
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        img_rows = self.cur.fetchall()
        sorted_img_rows = sorted(img_rows, key=lambda tup: tup[0])
        self.cur.execute('''SELECT COUNT(1) FROM assoc_source
            WHERE ndetect = 1''')
        single_detections = self.cur.fetchone()[0]
        self.cur.execute('''SELECT COUNT(1) FROM assoc_source
            WHERE ndetect = 2''')
        associations = self.cur.fetchone()[0]
        result = [sorted_img_rows, single_detections, associations]
        self.assertEqual(result, [[(1, 3), (2, 3)], 14, 5])
        self.cur.close()


    def test_saonly_already_done(self):
        """Branch 18 - SA already done, so do nothing"""
        # Process through SA first
        # sf, sa, cm
        stages = (True, True, False)
        # save, qa, overwrite, reprocess, redo match, update match
        # Testing with reprocess True to make sure it doesn't matter
        opts = (True, False, False, True, False, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Now try SA only
        stages = (False, True, False)
        # Code should exit
        self.assertIsNone(process(self.conn, stages, opts, self.dirs,
                                  self.catalogs, self.params, self.res_tol))


    def test_saonly(self):
        """Branch 19 - SA only"""
        # Process through SF first
        # sf, sa, cm
        stages = (True, False, False)
        # save, qa, overwrite, reprocess, redo match, update match
        opts = (True, False, False, False, False, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Now do SA only
        stages = (False, True, False)
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params, self.res_tol)
        # Check DB - should be same as SF + SA (branch 11)
        self.cur = self.conn.cursor()
        self.cur.execute('SELECT id, stage FROM image')
        img_rows = self.cur.fetchall()
        sorted_img_rows = sorted(img_rows, key=lambda tup: tup[0])
        self.cur.execute('''SELECT COUNT(1) FROM assoc_source
            WHERE ndetect = 1''')
        single_detections = self.cur.fetchone()[0]
        self.cur.execute('''SELECT COUNT(1) FROM assoc_source
            WHERE ndetect = 2''')
        associations = self.cur.fetchone()[0]
        result = [sorted_img_rows, single_detections, associations]
        self.assertEqual(result, [[(1, 3), (2, 3)], 24, 10])
        self.cur.close()


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSABranches)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
