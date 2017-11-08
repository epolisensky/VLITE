import psycopg2
import numpy as np
from database import createdb


conn = psycopg2.connect(host='localhost', database='test', user='erichards')
cur = conn.cursor()

createdb.create(cur, safe=True)

imagevalues = [(1, 'A'), (2, 'B'), (3, 'C')]
sql = 'INSERT INTO image (id, filename) VALUES (%s, %s);'
cur.executemany(sql, imagevalues)

assocvalues = [(1, np.mean([10., 10.1, 10.]), 3, 3),
               (2, np.mean([20., 20.1, 20.]), 3, 3),
               (3, np.mean([15., 15.1]), 2, 3), (4, 30.0, 1, 3)]
cur.executemany('''INSERT INTO assoc_source (id, ra, ndetect, nopp)
    VALUES (%s, %s, %s, %s);''', assocvalues)

rawisland = [(0, 1), (1, 1), (2, 1), (3, 1),
             (0, 2), (1, 2), (2, 2), (0, 3)]
cur.executemany('''INSERT INTO raw_island (isl_id, image_id)
    VALUES (%s, %s);''', rawisland)

rawvalues = [(0, 0, 1, 10.0, 1), (1, 1, 1, 20.0, 2), (2, 2, 1, 15.0, 3),
             (3, 3, 1, 30.0, 4), (0, 0, 2, 10.1, 1), (1, 1, 2, 20.1, 2),
             (2, 2, 2, 15.1, 3), (0, 0, 3, 10.0, 1), (1, 0, 3, 20.0, 2)]
cur.executemany('''INSERT INTO raw_source (
    src_id, isl_id, image_id, ra, assoc_id)
    VALUES (%s, %s, %s, %s, %s)''', rawvalues)

nullvalues = [(1, 4, 2), (2, 3, 3), (3, 4, 3)]
cur.executemany('''INSERT INTO null_detections (
    id, assoc_id, image_id)
    VALUES (%s, %s, %s)''', nullvalues)

conn.commit()
cur.close()

cur = conn.cursor()

# test remove image id 1
print('\nRemoving image id = 1...')
cur.execute('DELETE FROM image WHERE id = 1')

cur.execute('SELECT id, ra, ndetect, nopp FROM assoc_source')
rows = cur.fetchall()
sorted_rows = sorted(rows, key=lambda tup: tup[0])
print('\nassoc_source rows: {}'.format(sorted_rows))
exp = [(1, 10.05, 2, 2), (2, 20.05, 2, 2), (3, 15.1, 1, 2)]
if sorted_rows != exp:
    print('\nFailed on assoc_source.')

cur.execute('SELECT id, assoc_id, image_id FROM null_detections')
rows = cur.fetchall()
print('\nnull_detections rows: {}'.format(rows))
exp = [(2, 3, 3), ]
if rows != exp:
    print('\nFailed on null_detections.')

conn.rollback()

# test remove image id 2
print('\nRemoving image id = 2...')
cur.execute('DELETE FROM image WHERE id = 2')

cur.execute('SELECT id, ra, ndetect, nopp FROM assoc_source')
rows = cur.fetchall()
sorted_rows = sorted(rows, key=lambda tup: tup[0])
print('\nassoc_source rows: {}'.format(sorted_rows))
exp = [(1, 10., 2, 2), (2, 20., 2, 2), (3, 15., 1, 2), (4, 30., 1, 2)]
if sorted_rows != exp:
    print('\nFailed on assoc_source.')

cur.execute('SELECT id, assoc_id, image_id FROM null_detections')
rows = cur.fetchall()
print('\nnull_detections rows: {}'.format(rows))
exp = [(2, 3, 3), (3, 4, 3)]
if rows != exp:
    print('\nFailed on assoc_source.')

conn.rollback()

# test remove image id 3
print('\nRemoving image id = 3...')
cur.execute('DELETE FROM image WHERE id = 3')

cur.execute('SELECT id, ra, ndetect, nopp FROM assoc_source')
rows = cur.fetchall()
sorted_rows = sorted(rows, key=lambda tup: tup[0])
print('\nassoc_source rows: {}'.format(sorted_rows))
exp = [(1, 10.05, 2, 2), (2, 20.05, 2, 2), (3, 15.05, 2, 2), (4, 30., 1, 2)]
if sorted_rows != exp:
    print('\nFailed on assoc_source.')

cur.execute('SELECT id, assoc_id, image_id FROM null_detections')
rows = cur.fetchall()
print('\nnull_detections rows: {}'.format(rows))
exp = [(1, 4, 2)]
if rows != exp:
    print('\nFailed on assoc_source.')

conn.rollback()

cur.close()
conn.close()


