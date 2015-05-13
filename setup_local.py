from distutils.core import setup
from Cython.Build import cythonize

setup(name='parabam',
      description='Parallel BAM File Analysis',
      version='0.1.7dev',
      author="JHR Farmery",
      license='BSD',
      author_email='jhrf2@cam.ac.uk',
      packages=['parabam', 'parabam.command'],
      package_dir = {'parabam':'parabam','parabam.command':'parabam/command'},
      requires=['cython', 'numpy', 'argparse', 'pysam'],
      scripts=['parabam/bin/parabam'],
      ext_modules=cythonize(('parabam/chaser.pyx','parabam/core.pyx',
                             'parabam/merger.pyx','parabam/command/core.pyx',
                             'parabam/command/stat.pyx',
                             'parabam/command/subset.pyx')))
