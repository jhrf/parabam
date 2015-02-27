from setuptools import setup
from setuptools.extension import Extension
from setuptools.command.sdist import sdist as _sdist

cmdclass = {}
ext_modules = [
      Extension("parabam.core", [ "core.c" ]),
      Extension("parabam.chaser", ["chaser.c"]),
      Extension("parabam.merger", ["merger.c"]),
      Extension("parabam.command.subset", ["command/subset.c"]),
      Extension("parabam.command.stat", ["command/stat.c"]),
      Extension("parabam.command.core", ["command/core.c"])] 

class sdist(_sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        from Cython.Build import cythonize
        cythonize(['chaser.pyx','core.pyx','merger.pyx','command/core.pyx','command/stat.pyx','command/subset.pyx'])
        _sdist.run(self)
cmdclass['sdist'] = sdist

setup(name='parabam',
      description='Parallel BAM File Analysis',
      version='0.2.dev5',
      author="JHR Farmery",
      license='GPL',
      author_email = 'jhrf2@cam.ac.uk',
      packages = ['parabam','parabam.command'],
      package_dir = {'parabam':'','parabam.command':'command'},
      requires = ['numpy','argparse','pysam'],
      scripts = ['bin/parabam'],
      cmdclass = cmdclass,
      ext_modules=ext_modules
      )
