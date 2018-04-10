from setuptools import setup


setup(name='vdp',
      version='1.8',
      description='VLITE Database Pipeline',
      url='https://github.com/erichards/VLITE',
      author='Emily Richards',
      packages=['vdp', 'vdp/database', 'vdp/matching',
                'vdp/skycatalog', 'vdp/sourcefinding'])
