import os
from distutils.core import setup
from Cython.Build import cythonize

setup(name='parabam',
      description='Parallel BAM File Analysis',
      version='0.1.4dev',
      author="JHR Farmery",
      license='BSD',
      author_email = 'jhrf2@cam.ac.uk',
      packages = ['parabam','parabam.command'],
      package_dir = {'parabam': '','parabam.command': 'command'},
      package_data = {'parabam':['chaser.pyx',"support.pyx","core.pyx",],
                  'parabam.command':["subset.pyx","stat.pyx",]},
      requires = ['cython','numpy','argparse','pysam'],
      scripts = ['bin/parabam'],
      ext_modules=cythonize(("chaser.pyx","support.pyx","core.pyx","command/subset.pyx","command/stat.pyx"))#""
      )
