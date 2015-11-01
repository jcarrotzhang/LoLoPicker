from setuptools import setup

setup(name='lolopicker',
      version='0.1',
      description='somatic mutation caller',
      url='http://github.com/jcarrotzhang/LoLoPicker',
      author='Jian Carrot-Zhang',
      author_email='jian.carrot.zhang@gmail.com',
      license='McGill',
      packages=['lolopicker'],
      install_requires=[
          'pysam',
	  'pysamstats',
      ],
      zip_safe=False)
