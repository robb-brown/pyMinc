import numpy as np
cimport numpy as np

from pyMinctracc cimport *
from builtins cimport *
from volume_io cimport ALLOC, FREE

from pyMinc import VIOGeneralTransform, VIOVolume

from cpython cimport PyCapsule, PyCapsule_New, PyCapsule_GetPointer


cdef class Minctracc(object):

	translator = {
		# transformType
		'lsq3' : TRANS_LSQ3,
		'lsq6' : TRANS_LSQ6,
		'lsq7' : TRANS_LSQ7,
		'lsq9' : TRANS_LSQ9,
		'lsq10' : TRANS_LSQ10,
		'lsq12' : TRANS_LSQ12,
		'nonlinear' : TRANS_NONLIN,
		'PAT' : TRANS_PAT,
		'pat' : TRANS_PAT,
		
		# interpolating type
		'linear' : TRILINEAR,
		'cubic' : TRICUBIC,
		'nearest' : N_NEIGHBOUR,
		
		# metric
		'xcorr' : XCORR,
		'mi' : MUTUAL_INFORMATION,
		'nmi' : NORMALIZED_MUTUAL_INFORMATION,
		
		# optimizer
		'simplex' : OPT_SIMPLEX,
		'BFGS' : OPT_BFGS,
	}
			
	cdef copyArgs(self,Arg_Data *cArgs, args):
		for k,v in args.items():
			if k == 'debug':
				cArgs.flags.debug = v
			if k == 'verbose':
				cArgs.flags.verbose = v
			if k == 'transformType':
				cArgs.trans_info.transform_type = self.translator[v]
			if k == 'interpolation' :
				cArgs.interpolant_type = self.translator[v]
			if k == 'metric' :
				cArgs.obj_function_type = self.translator[v]
			if k == 'step' :
				cArgs.step[0] = float(v[0]); cArgs.step[1] = float(v[1]); cArgs.step[2] = float(v[2])
			if k == 'optimizer':
				cArgs.optimize_type = self.translator[v]
			
	
	
	def minctracc(self,source,target,sourceMask=None,targetMask=None,initialXFM=None,iterations=4,**args):
		
		cdef Arg_Data *cArgs = NULL
		ALLOC(cArgs,1)
		initializeArgs(cArgs)
		
		self.copyArgs(cArgs,args)
		
		if cArgs.trans_info.transform_type == TRANS_NONLIN:
			smd = source.metadata
			if not smd['dtype'] == np.float64:
				source = VIOVolume(source.data.astype(np.float64),**smd)
			tmd = target.metadata
			if not tmd['dtype'] == np.float64:
				target = VIOVolume(target.data.astype(np.float64),**tmd)

		cdef VIO_General_transform *initial = NULL
		cdef VIO_General_transform *final = NULL
			
		if initialXFM:
			initial = <VIO_General_transform *>PyCapsule_GetPointer(initialXFM.transformPtr,NULL)
			
		cdef VIO_Volume sourceC, targetC, sourceMaskC, targetMaskC = NULL
		
		sourceC = <VIO_Volume>PyCapsule_GetPointer(source.volumePtr,NULL); targetC = <VIO_Volume>PyCapsule_GetPointer(target.volumePtr,NULL)
		
		sourceMaskC = <VIO_Volume>PyCapsule_GetPointer(sourceMask.volumePtr,NULL) if sourceMask else NULL
		targetMaskC = <VIO_Volume>PyCapsule_GetPointer(targetMask.volumePtr,NULL) if targetMask else NULL
		
		final = minctracc(sourceC,targetC,sourceMaskC,targetMaskC,initial,iterations,cArgs)

		finalXFM = VIOGeneralTransform(PyCapsule_New(<void*>final,NULL,NULL))
		
		if not cArgs.trans_info.transform_type == TRANS_NONLIN:
			finalXFM.calculateInverseLinearTransform()

		return finalXFM
		

