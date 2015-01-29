import os
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from distutils.command.sdist import sdist as _sdist

cmdclass = {}
ext_modules = [
      Extension("parabam.core", [ "core.c" ]),
      Extension("parabam.chaser", ["chaser.c"]),
      Extension("parabam.support", ["support.c"]),
      Extension("parabam.command.subset", ["command/subset.c"]),
      Extension("parabam.command.stat", ["command/stat.c"])
]

class sdist(_sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        cythonize(['chaser.pyx','core.pyx','support.pyx','command/stat.pyx','command/subset.pyx'])
        _sdist.run(self)
cmdclass['sdist'] = sdist

setup(name='parabam',
      description='Parallel BAM File Analysis',
      version='0.1.11dev',
      author="JHR Farmery",
      license='BSD',
      author_email = 'jhrf2@cam.ac.uk',
      packages = ['parabam','parabam.command'],
      package_dir = {'parabam':'','parabam.command':'command'},
      #package_data = {'parabam':['chaser.pyx',"support.pyx","core.pyx",],
      #            'parabam.command':["subset.pyx","stat.pyx",]},
      requires = ['cython','numpy','argparse','pysam'],
      scripts = ['bin/parabam'],
      cmdclass = cmdclass,
      ext_modules=ext_modules#""
      )
