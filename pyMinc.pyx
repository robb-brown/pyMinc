import numpy as np
cimport numpy as np

from cpython.ref cimport PyObject
from cpython cimport PyCapsule

from builtins cimport *
from netCDF cimport *
from libminc cimport *
from volume_io cimport *

np.import_array()


ncTypeToNumpy = {	NC_BYTE : {True:np.int8, False:np.uint8},
					NC_CHAR : {True:np.int8, False:np.uint8},
					NC_SHORT : {True:np.int16, False:np.uint16},
					NC_INT : {True:np.int32, False:np.uint32},
					NC_FLOAT : {True:np.float32, False:np.float32},
					NC_DOUBLE : {True:np.float64, False:np.float64},
				}		
				


# libminc MINC file reading

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
		outputType = NC_FLOAT
		outputSigned = 1
		cdef nc_type mdatatype = -1
		cdef int issigned = -1
		cdef np.ndarray validRange = np.zeros(2,np.float64)
		cdef np.ndarray defaultRange = np.zeros(2,np.float64)
		cdef np.ndarray actualRange = np.zeros(2,np.float64)
		
#		print "In setupICV"
		self.imgid = ncvarid(self.mincFile,MIimage);
#		print self.imgid
#		print "Getting datatype"
		miget_datatype(self.mincFile,self.imgid,&mdatatype,&issigned)
#		issigned = not issigned
#		print "Getting valid range"
		miget_valid_range(self.mincFile,self.imgid,<double *>validRange.data)
		miget_default_range(outputType,outputSigned,<double *>defaultRange.data)
#		print "getting image range"
		miget_image_range(self.mincFile,<double *>actualRange.data)
		
		# Specify the type of image to get.  Just use float and fix later
		miicv_setint(self.icv,MI_ICV_TYPE,outputType)
#		miicv_setstr(self.icv,MI_ICV_SIGN,outputSigned)
		miicv_setdbl(self.icv, MI_ICV_VALID_MIN,defaultRange[0])
		miicv_setdbl(self.icv, MI_ICV_VALID_MAX,defaultRange[1])
		miicv_setint(self.icv, MI_ICV_DO_NORM,True)
		miicv_setint(self.icv, MI_ICV_USER_NORM,True)
		
		# Specify image dimensions.  We want positive values (do we?)
		miicv_setint(self.icv, MI_ICV_DO_DIM_CONV, True);
		miicv_setint(self.icv, MI_ICV_XDIM_DIR, MI_ICV_POSITIVE);
		miicv_setint(self.icv, MI_ICV_YDIM_DIR, MI_ICV_POSITIVE);
		miicv_setint(self.icv, MI_ICV_ZDIM_DIR, MI_ICV_POSITIVE);
			
		self.mdatatype = mdatatype
		self.issigned = issigned
		self.validRange = validRange
		self.actualRange = actualRange
		
#		print "MDatatype: ", mdatatype
#		print "Is signed: ", issigned
#		print "valid range: ",validRange
#		print "actual range: ",actualRange
#		print "default range: ",defaultRange
		

	def loadFile(self,fname=None):
		if sizeof(long) == 4:
			longType = np.int32
		else:
			longType = np.int64
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
		
#		print "Sizeof long: %d" % sizeof(long)
		self.imgid = imgid
		if (fname):
			self.fname = fname
		self.mincFile = miopen(self.fname,NC_NOWRITE)
#		print "Creating ICV"
		self.icv = miicv_create()
		self.setupICV()

#		print "Done setting up ICV"
		miicv_attach(self.icv,self.mincFile,ncvarid(self.mincFile,MIimage))
#		print "Attaching ICV"
		ret = miicv_inqint(self.icv, MI_ICV_VARID, &imgid)
		# get num of dimensions
#		print "Getting num of dimensions"
		ncvarinq(self.mincFile,imgid,name,&nctype,&ndims,<int*>dims.data,NULL)

		cdef np.ndarray origin = np.zeros(MAX_VAR_DIMS,longType)
		cdef np.ndarray count = np.zeros(MAX_VAR_DIMS,longType)
#		print "Getting dimlengths and names"
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
		self.min = min
		self.max = max
		
		cdef int mdatatype
		cdef char sigd[10]
		
		miicv_inqint(self.icv, MI_ICV_TYPE, &mdatatype);
		miicv_inqstr(self.icv, MI_ICV_SIGN, sigd);
		
#		print mdatatype
#		print sigd
		
		datatype = None
		typesize = -1

		datatype = ncTypeToNumpy[mdatatype][sigd == MI_SIGNED]
		self.datatype = datatype
		
		miset_coords(ndims,0,<long*>origin.data)
		cdef np.ndarray data = np.zeros(dimLengths,datatype)

#		print "Getting data"		
		miicv_get(self.icv,<long*>origin.data,<long*>count.data,<void*>data.data)
#		print "Done getting data"
		self.data = data
		
#		print "Closing file"
		miclose(self.mincFile)
#		print "Freeing icv"
		miicv_free(self.icv)
		

	def getNumpyArray(self):
		return self.data
		
		
		
		
# VIO_Volume wrapper		
cdef class VIOVolume:
	
	cdef VIO_Volume volume
	cdef minc_input_options *options
		
	def __init__(self):
		pass


	cdef alloc_volume_data(self,VIO_Volume volume):
		cdef unsigned int data_size = get_volume_total_n_voxels( volume )
		cdef unsigned int type_size = get_type_size( get_volume_data_type( volume ) )
		volume.is_cached_volume = False
#		self.alloc_multidim_array(volume.array)
	    		
		
	def invent(self):
		cdef VIO_Volume vol = NULL
		cdef char **dim_names
		cdef np.ndarray sizes = np.array([10,10,10],np.int32)
		
		dim_names = get_default_dim_names(3)
		ALLOC(vol,1)
		vol = create_volume(3,dim_names,NC_BYTE,False,0,255)
		set_volume_sizes(vol,<int*>(sizes.data))
		
		voxels = get_volume_total_n_voxels(vol)
		size = get_type_size(get_volume_data_type(vol))
		dtype = ncTypeToNumpy[vol.nc_data_type][vol.signed_flag]
		
		alloc_volume_data(vol)

		cdef np.ndarray arr = np.zeros(sizes,dtype)
		arr[0,0,0] = 1; arr[5,5,5] = 2; arr[9,9,9] = 3
		cdef void *dataPtr = NULL
		GET_VOXEL_PTR(dataPtr,vol,0,0,0,0,0)
		
		memcpy(<char*>dataPtr,<char*>arr.data,1000)
		
#		set_volume_separations
#		set_volume_starts
#		set_volume_type
		
		status = output_volume("/temp/scrap.mnc",NC_BYTE,False,0,0,vol,NULL,NULL)
		return status
		
		
		
	def read(self,fname,dtype=NC_UNSPECIFIED,dsigned=False,dmin=0,dmax=255,create=True):
		cdef VIO_Volume vol = NULL
		cdef minc_input_options *options = NULL

		ALLOC(vol,1)
		dtype = MI_ORIGINAL_TYPE; dsigned = False
		status = input_volume(fname,0,NULL,dtype,dsigned,dmin,dmax,create,&vol,options)
		self.volume = vol
		self.options = options
		
		status = output_volume("/temp/scrap.mnc",dtype,dsigned,dmin,dmax,vol,NULL,NULL)
		
		return status
		
	property volumePtr:
	def __get__(self):
		return PyCapsule_New(<void *>self.volume, NULL, NULL)

	property data:
		def __get__(self):
			cdef int i
			volume = {}
			volume['dtype'] = ncTypeToNumpy[self.volume.nc_data_type][self.volume.signed_flag]
			volume['min'] = self.volume.voxel_min
			volume['max'] = self.volume.voxel_max
			volume['dimensions'] = self.volume.array.n_dimensions
			volume['names'] = []
			cdef np.ndarray tempArr = np.zeros(3,np.int64)
			volume['shape'] = tempArr
			volume['spacing'] = []
			volume['starts'] = []
			volume['cosines'] = []
			for i in range(0,volume['dimensions']):
				volume['shape'][i] = self.volume.array.sizes[i]
				volume['names'].append(self.volume.dimension_names[i])
				volume['spacing'].append(self.volume.separations[i])
				volume['starts'].append(self.volume.starts[i])
				volume['cosines'].append([])
				for m in range(0,VIO_N_DIMENSIONS):
					volume['cosines'][i].append(self.volume.direction_cosines[i][m])

			cdef void *dataPtr = NULL
			GET_VOXEL_PTR(dataPtr,self.volume,0,0,0,0,0)

			typenum = np.dtype(volume['dtype']).num

			cdef np.ndarray data = np.PyArray_SimpleNewFromData(volume['dimensions'],<int *>tempArr.data,typenum,dataPtr)

			print np.shape(data)
			print volume['shape']

			# USE THIS TO COPY THE DATA TO PYTHON INSTEAD OF JUST WRAPPING A NUMPY ARRAY AROUND IT
	#		np.zeros(np.product(volume['shape']),volume['dtype'])
	#		count=data.dtype.itemsize*np.product(volume['shape'])
	#		memcpy(<char*>data.data,<char*>dataPtr,count)
			volume['data'] = data

			return volume
		
		
