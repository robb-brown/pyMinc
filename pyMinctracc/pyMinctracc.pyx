import numpy as np
cimport numpy as np

from pyMinctracc cimport *
from builtins cimport *
from volume_io cimport ALLOC, FREE

from pyMinc import VIOGeneralTransform

from cpython cimport PyCapsule, PyCapsule_New, PyCapsule_GetPointer


cdef class Minctracc(object):
	
	cdef copyArgs(self,Arg_Data *cArgs, args):
		for k,v in args.items():
			if k == 'debug':
				cArgs.flags.debug = v
			if k == 'verbose':
				cArgs.flags.verbose = v
			if k == 'transformType':
				if v == 'lsq3': cArgs.trans_info.transform_type = TRANS_LSQ3
				if v == 'lsq6': cArgs.trans_info.transform_type = TRANS_LSQ6
				if v == 'lsq7': cArgs.trans_info.transform_type = TRANS_LSQ7
				if v == 'lsq9': cArgs.trans_info.transform_type = TRANS_LSQ9
				if v == 'lsq12': cArgs.trans_info.transform_type = TRANS_LSQ12
			if k == 'optimizer':
				if v == 'simplex': cArgs.optimize_type = OPT_SIMPLEX
				if v == 'BFGS': cArgs.optimize_type = OPT_BFGS
	
	
	def minctracc(self,source,target,sourceMask=None,targetMask=None,initialXFM=None,**args):
		
		cdef Arg_Data *cArgs
		ALLOC(cArgs,1)
		initializeArgs(cArgs)
		
		self.copyArgs(cArgs,args)

		cdef VIO_General_transform *initial = NULL
		cdef VIO_General_transform *final = NULL
			
		if initialXFM:
			initial = <VIO_General_transform *>PyCapsule_GetPointer(initialXFM.transformPtr,NULL)
			
		cdef VIO_Volume sourceC, targetC, sourceMaskC, targetMaskC = NULL
		
		sourceC = <VIO_Volume>PyCapsule_GetPointer(source.volumePtr,NULL); targetC = <VIO_Volume>PyCapsule_GetPointer(target.volumePtr,NULL)
		
		sourceMaskC = <VIO_Volume>PyCapsule_GetPointer(sourceMask.volumePtr,NULL) if sourceMask else NULL
		targetMaskC = <VIO_Volume>PyCapsule_GetPointer(targetMask.volumePtr,NULL) if targetMask else NULL
		
		final = minctracc(sourceC,targetC,sourceMaskC,targetMaskC,initial,cArgs)

		finalXFM = VIOGeneralTransform(PyCapsule_New(<void*>final,NULL,NULL))

		return finalXFM
		

