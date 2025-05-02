from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy,os


#libraries = ['minc2','netcdf','hdf5','curl','minctracc']# 'volume_io2' - volume_io2 is in minc2 in newer versions
libraries = ['minc2','netcdf','hdf5','curl']# 'volume_io2' - volume_io2 is in minc2 in newer versions
include_dirs = ['.','pyMinctracc',os.path.expandvars('$MINC_TOOLKIT/include'),'/usr/local/include','/usr/local/minc-toolkit/include']
library_dirs = ['.','pyMinctracc',os.path.expandvars('$MINC_TOOLKIT/lib'),'/usr/lib','/usr/local/lib','/usr/local/minc-toolkit/lib']

extensions = [
	Extension("pyMinc", ['pyMinc.pyx'],
						libraries=libraries,
						include_dirs = include_dirs,
						library_dirs = library_dirs,
						),
						
	# Extension("pyMinctracc", ['pyMinctracc/pyMinctracc.pyx'],
	# 					libraries=libraries,
	# 					include_dirs = include_dirs,
	# 					library_dirs = library_dirs,
	# 					),
]

setup(
	name = 'pyMinc',
	version = "0.2",
	description="A Python module to interface with MINC",
	author = 'Robert Brown',
	author_email="robb@robbtech.com",
	url='http://www.robbtech.com',
	#packages = ['pyMinctracc'],
	include_dirs = [numpy.get_include()],
	ext_modules = extensions,
  	cmdclass = {'build_ext': build_ext},
)
