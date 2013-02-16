import numpy as np
cimport numpy as np

from pyMinctracc cimport *
from builtins cimport *
from cython import PyString_AsString

class Minctracc(object):
	
	def exectute(command):
		argv = command.split(' ')
		argc = len(argv)
		cdef char **strings = <char**>malloc(argc*sizeof(char*))
		for i,s in enumerate(argv):
			strings[i] = argv[i]
		
		ret = minctracc(argc,strings)
		return ret
		

