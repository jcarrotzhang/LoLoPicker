from setuptools import setup

try:
    import pysam
except ImportError:
    raise Exception('pysam not found; please install pysam first')
try:
    import pysamstats
except ImportError:
    raise Exception('pysamstats not found; please install pysamstats first')
try:
	import scipy
except ImportError:
	raise Exception('scipy not found; please install scipy first')

from distutils.version import StrictVersion
required_pysam_version = '0.8.4'
if StrictVersion(pysam.__version__) < StrictVersion(required_pysam_version):
    raise Exception('pysam version >= %s is required; found %s' %
                    (required_pysam_version, pysam.__version__))

setup(name='lolopicker',
      version='0.4',
      description='somatic mutation caller',
      url='http://github.com/jcarrotzhang/LoLoPicker',
      author='Jian Carrot-Zhang',
      author_email='jian.carrot.zhang@gmail.com',
      license='McGill',
      scripts=[ 'scripts/LoLoPicker_control.py',
		'scripts/LoLoPicker_stats.py',
        	'scripts/LoLoPicker_somatic.py'
    		],
      packages=['lolopicker'],
      install_requires=[
          'pysam',
	  'pysamstats'
      ],
      zip_safe=False)


