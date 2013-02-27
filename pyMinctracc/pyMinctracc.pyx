import numpy as np
cimport numpy as np

from pyMinctracc cimport *
from builtins cimport *
from volume_io cimport VIO_Volume


class Minctracc(object):
	
	def minctracc(self,source,target,sourceMask=None,targetMask=None,initialXFM=None,command=''):
		argv = command.split(' ')
		argc = len(argv)
		cdef char **strings = <char**>malloc(argc*sizeof(char*))
		for i,s in enumerate(argv):
			strings[i] = argv[i]
			
		
		cdef VIO_General_transform *initial = NULL
		cdef VIO_General_transform *final = NULL
		cdef VIO_Volume sourceC, targetC, sourceMaskC, targetMaskC = NULL
		
		sourceC = source.getVolumePtr(); targetC = target.getVolumePtr()
		
		sourceC = sourceMask.getVolumePtr() if sourceMask else NULL
		targetC = targetMask.getVolumePtr() if targetMask else NULL
		
		finalXFM = minctracc(sourceC,targetC,sourceMaskC,targetMaskC,initial,final,argc,strings)
		
		return finalXFM
		

