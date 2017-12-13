import psycopg2
import sys
from errors import ConfigError
from yaml import load
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader


def cluster():
    try:
        cf = sys.argv[1]
    except IndexError:
        raise ConfigError('Please provide a configuration file.')

    with open(cf, 'r') as stream:
        data = load(stream, Loader=Loader)
        dbname = (data['setup'])['dbname']
        dbusr = (data['setup'])['dbuser']
    
    try:
        conn = psycopg2.connect(host='localhost', database=dbname, user=dbusr)
    except:
        raise ConfigError('Could not connect to database.')
        
    try:
        cur = conn.cursor()
        cur.execute(
            'CLUSTER detected_source_q3c_ang2ipix_idx ON detected_source;')
        cur.execute('ANALYZE detected_source;')
        cur.execute('CLUSTER assoc_source_q3c_ang2ipix_idx ON assoc_source;')
        cur.execute('ANALYZE assoc_source;')
        cur.close()
    except:
        raise ConfigError('Tables could not be clustered.')


if __name__ == '__main__':
    cluster()
