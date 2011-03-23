from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy


sourcefiles = ['pyMinc.pyx']

ext_modules = [
	Extension("mincFile", sourcefiles,
						libraries=['minc','netcdf','curl'],
						include_dirs = ['/usr/local/include','/usr/local/mni/include'],
						library_dirs = ['/usr/local/lib','/opt/local32/lib'],
						)
	]

setup(
	name = 'pyMinc',
	version = "0.1",
	description="A Python module to interface with MINC",
	author = 'Robert Brown',
	author_email="robb@robbtech.com",
	url='http://www.robbtech.com',
  cmdclass = {'build_ext': build_ext},
	include_dirs = [numpy.get_include()],
  ext_modules = ext_modules
)
