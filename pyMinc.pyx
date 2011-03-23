import numpy as np
cimport numpy as np

from pyMinc.pyMinc cimport *

from mincConstants import *


ncTypeToNumpy = {	NC_BYTE : {True:np.int8, False:np.uint8},
					NC_CHAR : {True:np.int8, False:np.uint8},
					NC_SHORT : {True:np.int16, False:np.uint16},
					NC_INT : {True:np.int32, False:np.uint32},
					NC_FLOAT : {True:np.float32, False:np.float32},
					NC_DOUBLE : {True:np.float64, False:np.float64},
				}

class mincFile(object):

	def __init__(self,fname=None):
		self.fname = fname
		self.mincFile = -1
		self.icv = -1
		self.validRange = None
		self.actualRange = None
		self.dims = None
		self.data = None
		if (fname):
			self.loadFile()

	def setupICV(self):
		cdef int mdatatype
		cdef int issigned
		cdef np.ndarray validRange = np.zeros(2,np.float64)
		cdef np.ndarray actualRange = np.zeros(2,np.float64)

#		print "Getting datatype"
#		miget_datatype(self.mincFile,self.imgid,<nc_type *>mdatatype,<int *>issigned)
#		issigned = not issigned
#		print "Getting valid range"
#		miget_valid_range(self.mincFile,self.imgid,<double *>validRange.data)
#		print "getting image range"
#		miget_image_range(self.mincFile,<double *>actualRange.data)
		
		# Specify the type of image to get.  Just use float and fix later
		miicv_setint(self.icv,MI_ICV_TYPE,NC_FLOAT)
		
		# Specify image dimensions.  We want positive values (do we?)
		miicv_setint(self.icv, MI_ICV_DO_DIM_CONV, True);
		miicv_setint(self.icv, MI_ICV_XDIM_DIR, MI_ICV_POSITIVE);
		miicv_setint(self.icv, MI_ICV_YDIM_DIR, MI_ICV_POSITIVE);
		miicv_setint(self.icv, MI_ICV_ZDIM_DIR, MI_ICV_POSITIVE);
			
		self.mdatatype = mdatatype
		self.issigned = issigned
		self.validRange = validRange
		self.actualRange = actualRange
		
#		print mdatatype
#		print issigned
#		print validRange
#		print actualRange
		

	def loadFile(self,fname=None):
		cdef int imgid = -1
		cdef int ndims = -1
		cdef np.ndarray dims = np.zeros(10,np.int32)
		cdef nc_type nctype = -1
		cdef char name[50]
		cdef long dimLength = -1
		dimLengths = []
		dimNames = []
		cdef double min
		cdef double max
		cdef int ret
		
		self.imgid = imgid
		if (fname):
			self.fname = fname
		self.mincFile = miopen(self.fname,NC_NOWRITE)
		self.icv = miicv_create()
		self.setupICV()
		miicv_attach(self.icv,self.mincFile,ncvarid(self.mincFile,MIimage))
		ret = miicv_inqint(self.icv, MI_ICV_VARID, &imgid)
		# get num of dimensions
		ncvarinq(self.mincFile,imgid,name,&nctype,&ndims,<int*>dims.data,NULL)

		cdef np.ndarray origin = np.zeros(MAX_VAR_DIMS,np.int32)
		cdef np.ndarray count = np.zeros(MAX_VAR_DIMS,np.int32)
		for i in range(0,ndims):
			ncdiminq(self.mincFile,dims[i],name,&dimLength)
			dimLengths.append(dimLength)
			dimNames.append(name)
			count[i] = dimLength
#		print dims
#		print dimLengths
#		print dimNames
		
		# get the min and max values
		miicv_inqdbl(self.icv,MI_ICV_NORM_MIN,&min)
		miicv_inqdbl(self.icv,MI_ICV_NORM_MIN,&max)
#		print min
#		print max
		
		cdef int mdatatype
		cdef char sigd[10]
		
		miicv_inqint(self.icv, MI_ICV_TYPE, &mdatatype);
		miicv_inqstr(self.icv, MI_ICV_SIGN, sigd);
		
#		print mdatatype
#		print sigd
		
		datatype = None
		typesize = -1

		datatype = ncTypeToNumpy[mdatatype][sigd == MI_SIGNED]
		
		miset_coords(ndims,0,<long*>origin.data)
		cdef np.ndarray data = np.zeros(dimLengths,datatype)
		miicv_get(self.icv,<long*>origin.data,<long*>count.data,<void*>data.data)
		
		self.data = data
		
		miclose(self.mincFile)
		miicv_free(self.icv)
		

	def getNumpyArray(self):
		return self.data