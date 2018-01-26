"""Integration tests for all the possible paths
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
        self.conn = dbinit('branchtest', 'erichards', True, True)


    def tearDown(self):
        """Disconnect the database."""
        self.conn.close()


    def test_sfsa_new_image(self):
        """Branch 11 - SF + SA on new image"""
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        # Pipeline should stop after source association
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
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
        self.assertEqual(result, [[(1, 3), (2, 3)], 22, 11])
        self.cur.close()


    def test_sfsa_reprocess(self):
        """Branch 14 - SF + SA reprocess"""
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : True,
                'redo match' : False, 'update match' : False}
        # Process once
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
        # Use different scale so results are different
        self.params = {'mode' : 'default', 'thresh' : 'hard', 'scale' : 0.3}
        # Process twice
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
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
        self.assertEqual(result, [[(1, 3), (2, 3)], 12, 6])
        self.cur.close()


    def test_saonly_already_done(self):
        """Branch 18 - SA already done, so do nothing"""
        # Process through SA first
        stages = {'source finding' : True, 'source association' : True,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : True,
                'redo match' : False, 'update match' : False}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
        # Now try SA only
        stages = {'source finding' : False, 'source association' : True,
                  'catalog matching' : False}
        # Code should exit
        self.assertIsNone(process(self.conn, stages, opts, self.dirs,
                                  self.catalogs, self.params))


    def test_saonly(self):
        """Branch 19 - SA only"""
        # Process through SF first
        stages = {'source finding' : True, 'source association' : False,
                  'catalog matching' : False}
        opts = {'save to database' : True, 'quality checks' : True,
                'overwrite' : False, 'reprocess' : False,
                'redo match' : False, 'update match' : False}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
        # Now do SA only
        stages = {'source finding' : False, 'source association' : True,
                  'catalog matching' : False}
        process(self.conn, stages, opts, self.dirs,
                self.catalogs, self.params)
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
        self.assertEqual(result, [[(1, 3), (2, 3)], 22, 11])
        self.cur.close()


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSABranches)
    unittest.TextTestRunner(verbosity=2).run(suite)
    #unittest.main()
