from distutils.core import setup
from distutils.extension import Extension
import numpy


sourcefiles = ['pyMinc.c']

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
	include_dirs = [numpy.get_include()],
	ext_package = 'pyMinc',
  	ext_modules = ext_modules,
	py_modules = ['pyMinc.examples','pyMinc.mincConstants']

)
