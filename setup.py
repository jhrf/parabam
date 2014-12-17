import os
from distutils.core import setup

setup(name='parabam',
	  description='Parallel BAM File Analysis',
      version='0.1.1dev',
      author="JHR Farmery",
      license='BSD',
      author_email = 'jhrf2@cam.ac.uk',
      package_dir = {'parabam': '','parabam.interface': 'interface'},
      packages = ['parabam','parabam.interface'],
      requires = ['numpy','argparse','pysam'],
      package_data = {'parabam.interface':['master',]},
      scripts = ['bin/parabam']
      )
