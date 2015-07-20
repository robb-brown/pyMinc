from distutils.core import setup
from distutils.extension import Extension
import numpy


sourcefiles = ['pyMinc.c','pyMinctracc.c']
libraries = ['minc2','volume_io2','netcdf','hdf5','curl','minctracc']
include_dirs = ['.','pyMinctracc','/usr/local/include','/usr/local/minc-toolkit/include']
library_dirs = ['.','pyMinctracc','/usr/lib','/usr/local/lib','/usr/local/minc-toolkit/lib']

ext_modules = [
	Extension("pyMinc", ['pyMinc.c'],
		libraries=libraries,
		include_dirs = include_dirs,
		library_dirs = library_dirs,
	),
	Extension("pyMinctracc", ['pyMinctracc/pyMinctracc.c'],
		libraries=libraries,
		include_dirs = include_dirs,
		library_dirs = library_dirs,
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
	py_modules = ['pyMinc.examples','pyMinc.utilities']
)