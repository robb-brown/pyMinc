import numpy as np
cimport numpy as np

from pyMinctracc cimport *
from builtins cimport *

class Minctracc(object):
	
	def minctracc(self,command):
		argv = command.split(' ')
		argc = len(argv)
		cdef char **strings = <char**>malloc(argc*sizeof(char*))
		for i,s in enumerate(argv):
			strings[i] = argv[i]
		
		ret = minctracc(argc,strings)
		return ret
		

