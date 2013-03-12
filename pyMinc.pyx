import numpy as np
cimport numpy as np

from cpython.ref cimport PyObject
from cpython cimport PyCapsule, PyCapsule_New, PyCapsule_GetPointer, PyCapsule_IsValid

from builtins cimport *
from netCDF cimport *
from libminc cimport *
from volume_io cimport *
import traceback

np.import_array()


ncTypeToNumpy = {	NC_BYTE : {True:np.int8, False:np.uint8},
					NC_CHAR : {True:np.int8, False:np.uint8},
					NC_SHORT : {True:np.int16, False:np.uint16},
					NC_INT : {True:np.int32, False:np.uint32},
					NC_FLOAT : {True:np.float32, False:np.float32},
					NC_DOUBLE : {True:np.float64, False:np.float64},
				}
				
numpyToNCType = {	np.int8	: (NC_BYTE,True), 	np.uint8	: (NC_BYTE,False), 
					np.int16 : (NC_SHORT,True), 	np.uint16	: (NC_SHORT,False), 
					np.int32 : (NC_INT,True), 	np.uint32	: (NC_INT,False), 
					np.float32 : (NC_FLOAT,True),
					np.float64 : (NC_DOUBLE,True),
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
	cdef int owner
		
	def __init__(self,data=None,**args):
		self.owner = True
		if not data == None:
			if data.__class__ == str:
				self.read(data,**args)
			else:
				self.createWithData(data,args.get('spacing',None),args.get('starts',None),args.get('names',None))
				
	def __del__(self):
		if not self.volume == NULL and self.owner:
			FREE(self.volume)
		if not self.options == NULL:
			FREE(self.options)
			
			
	cdef setVolumePtr(self,VIO_Volume vol,int owner=1):
		self.volume = vol
		self.owner = owner
		

	cdef alloc_volume_data(self,VIO_Volume volume):
		cdef unsigned int data_size = get_volume_total_n_voxels( volume )
		cdef unsigned int type_size = get_type_size( get_volume_data_type( volume ) )
		volume.is_cached_volume = False
#		self.alloc_multidim_array(volume.array)


	cdef createWithData(self,np.ndarray data,spacing=None,starts=None,names=None):
		cdef VIO_Volume vol = NULL
		cdef char **dim_names = NULL
		cdef np.ndarray vec
		
		
		dimN = len(np.shape(data))

		if spacing == None:
			spacing = np.ones(dimN,np.float64)
		else:
			spacing = np.array(spacing,np.float64)

		if starts == None:
			starts = np.zeros(dimN,np.float64)
		else:
			starts = np.array(starts,np.float64)
		
		if names == None:
			dim_names = get_default_dim_names(dimN)
		else:
			ALLOC(dim_names,dimN)
			for i,s in enumerate(names):
				dim_names[i] = s
		type,signed = numpyToNCType[data.dtype.type]
		
		# Setup volume
		ALLOC(vol,1)
		vol = create_volume(dimN,dim_names,type,signed,data.min(),data.max())
		vec = np.array(np.shape(data),np.int32); set_volume_sizes(vol,<int*>(vec.data))
		vec = spacing; set_volume_separations(vol,<VIO_Real*>(vec.data))
		vec = starts; set_volume_starts(vol,<VIO_Real*>(vec.data))
		vol.voxel_to_world_transform_uptodate = False
		alloc_volume_data(vol)
		
		# Do the copy
		
		voxels = get_volume_total_n_voxels(vol)
		size = get_type_size(get_volume_data_type(vol))
		cdef void *dataPtr = NULL
		GET_VOXEL_PTR(dataPtr,vol,0,0,0,0,0)
		
		data = np.require(data,requirements='C')
		memcpy(<char*>dataPtr,<char*>data.data,voxels*size)
		
		self.volume = vol
		self.owner = True

		return True


		
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
		
#		alloc_volume_data(vol)

		cdef np.ndarray arr = np.zeros(sizes,dtype)
		arr[0,0,0] = 1; arr[5,5,5] = 2; arr[9,9,9] = 3
		cdef void *dataPtr = NULL
		GET_VOXEL_PTR(dataPtr,vol,0,0,0,0,0)

		arr = np.require(arr,requirements='F')
		memcpy(<char*>dataPtr,<char*>arr.data,1000)
		
#		set_volume_separations
#		set_volume_starts
#		set_volume_type
		
		status = output_volume("/temp/scrap.mnc",NC_BYTE,False,0,0,vol,NULL,NULL)
		return status
		
		
	def write(self,fname):
		type = self.volume.nc_data_type
		signed = self.volume.signed_flag
		status = output_volume(fname,type,signed,0,0,self.volume,NULL,NULL)
		
		
	def read(self,fname,type=MI_ORIGINAL_TYPE,dsigned=False,min=0,max=0,create=True):
		cdef VIO_Volume vol = NULL
		cdef minc_input_options *options = NULL
		type,signed = numpyToNCType.get(type,(type,dsigned))
		ALLOC(vol,1)
		status = input_volume(fname,0,NULL,type,signed,min,max,create,&vol,options)
		self.volume = vol
		self.options = options		
		return status
		
		
	def getVoxelToWorldTransform(self,cosines=False):
		cdef VIO_General_transform *xfm = NULL
		ALLOC(xfm,1)
		info = self.metadata
		compute_world_transform(self.volume.spatial_axes,self.volume.separations,self.volume.direction_cosines,self.volume.starts,xfm)
		transform = VIOGeneralTransform(); transform.setTransformPtr(xfm)
		return transform

	
	cdef voxelToWorldFast(self,np.ndarray point,np.ndarray transformed):
		convert_voxel_to_world(self.volume,<VIO_Real*>point.data,<VIO_Real*>transformed.data,<VIO_Real*>transformed.data+1,<VIO_Real*>transformed.data+2)


	def voxelToWorld(self,point):
		cdef np.ndarray point2 = np.array(point,np.float64)
		cdef np.ndarray transformed = np.zeros(3,np.float64)
		convert_voxel_to_world(self.volume,<VIO_Real*>point2.data,<VIO_Real*>transformed.data,<VIO_Real*>transformed.data+1,<VIO_Real*>transformed.data+2)
		return transformed
		

	def worldToVoxel(self,point):
		cdef np.ndarray transformed = np.zeros(self.volume.array.n_dimensions,np.float64)
		convert_world_to_voxel(self.volume,point[0],point[1],point[2],<VIO_Real*>transformed.data)
		return transformed
		
		
	property voxelToWorldTransform:
		def __get__(self):
			cdef VIO_General_transform *xfm = NULL
			xfm = get_voxel_to_world_transform(self.volume)
			transform = VIOGeneralTransform()
			transform.setTransformPtr(xfm); transform.owner = False
			return transform
		
		
	property volumePtr:
		def __get__(self):
			return PyCapsule_New(<void *>self.volume, NULL, NULL)
			

	property metadata:
		def __get__(self):
			cdef int i
			volume = {}
			volume['dtype'] = ncTypeToNumpy[self.volume.nc_data_type][self.volume.signed_flag]
			volume['min'] = self.volume.voxel_min
			volume['max'] = self.volume.voxel_max
			volume['dimensions'] = self.volume.array.n_dimensions
			volume['names'] = []
			cdef np.ndarray tempArr = np.zeros(volume['dimensions'],np.int64)
			volume['shape'] = tempArr
			volume['spacing'] = []
			volume['starts'] = []
			volume['cosines'] = []
			volume['spatialAxes'] = []
			for i in range(0,volume['dimensions']):
				volume['shape'][i] = self.volume.array.sizes[i]
				volume['names'].append(self.volume.dimension_names[i])
				volume['spacing'].append(self.volume.separations[i])
				volume['starts'].append(self.volume.starts[i])
				volume['cosines'].append([])
				volume['spatialAxes'].append(self.volume.spatial_axes[i])
				for m in range(0,VIO_N_DIMENSIONS):
					volume['cosines'][i].append(self.volume.direction_cosines[i][m])
			return volume			
			
			
	property data:
		def __get__(self):
			cdef int i
			volume = self.metadata
			cdef void *dataPtr = NULL
			cdef np.ndarray tempArr = np.zeros(volume['dimensions'],np.int64)
			GET_VOXEL_PTR(dataPtr,self.volume,0,0,0,0,0)
			typenum = np.dtype(volume['dtype']).num
			for i in range(0,volume['dimensions']):
				tempArr[i] = self.volume.array.sizes[i]
			cdef np.ndarray data = np.PyArray_SimpleNewFromData(volume['dimensions'],<np.npy_intp*>tempArr.data,typenum,dataPtr)

			# USE THIS TO COPY THE DATA TO PYTHON INSTEAD OF JUST WRAPPING A NUMPY ARRAY AROUND IT
	#		np.zeros(np.product(volume['shape']),volume['dtype'])
	#		count=data.dtype.itemsize*np.product(volume['shape'])
	#		memcpy(<char*>data.data,<char*>dataPtr,count)

			return data
	

		
# VIO_General_transform wrapper		
cdef class VIOGeneralTransform:

	cdef VIO_General_transform *transform
	cdef int owner


	def __init__(self,ptr=None,owner=True):
		if not ptr == None:
			if PyCapsule_IsValid(ptr,NULL):
				self.transform = <VIO_General_transform *>PyCapsule_GetPointer(ptr,NULL)
			elif ptr.__class__ == str:
				self.read(ptr)
			else:
				ALLOC(self.transform,1)
				self.setupTransform(self.transform,ptr)
		self.owner = owner


	cdef setupTransform(self,VIO_General_transform *xfmPtr,data,inverseFlag=False):		
		cdef VIO_Transform *ltransform = NULL
		cdef np.ndarray xfm
		cdef int i,inverse_flag
		if data.__class__ == VIOGeneralTransform:
			self.setupTransform(xfmPtr,data.getData(noInverse=True),data.inverseFlag)
		if data.__class__ == VIOVolume:
			xfmPtr.type = GRID_TRANSFORM
			xfmPtr.inverse_flag = inverseFlag
			xfmPtr.displacement_volume = <VIO_Volume>PyCapsule_GetPointer(data.volumePtr,NULL)
			xfmPtr.n_transforms = 1
		elif np.shape(data) == (4,4):
			xfmPtr.type = LINEAR
			xfmPtr.inverse_flag = inverseFlag
			xfm = np.array(data,np.float64)
			ALLOC(ltransform,1)
			xfm = np.require(xfm,requirements='F')
			memcpy(<char*>ltransform.m[0],<char*>xfm.data,16*8)	
			xfmPtr.linear_transform = ltransform
			ALLOC(ltransform,1)
			xfm = np.linalg.inv(xfm)
			memcpy(<char*>ltransform.m[0],<char*>xfm.data,16*8)	
			xfmPtr.inverse_linear_transform = ltransform
			xfmPtr.n_transforms = 1
		elif len(np.shape(data)) == 1:
			xfmPtr.type = CONCATENATED_TRANSFORM
			xfmPtr.inverse_flag = inverseFlag
			ALLOC(xfmPtr.transforms,len(data))
			for i,x in enumerate(data):
				inverse_flag = x.inverseFlag
				self.setupTransform(xfmPtr.transforms+i,x.getData(noInverse=True),inverse_flag)
			xfmPtr.n_transforms = len(data)


	def __del__(self):
		if not self.transform == NULL and self.owner:
			FREE(self.transform)
			
	
	cdef setTransformPtr(self,VIO_General_transform *ptr,int owner=1):
		self.transform = ptr
		self.owner = owner
			

	def write(self,fname,comments=''):
		status = output_transform_file(fname,comments,self.transform)
		return status


	def read(self,fname):
		cdef VIO_General_transform *xfm = NULL
		ALLOC(xfm,1)
		status = input_transform_file(fname,xfm)
		self.transform = xfm
		return status
		
		
	cdef transformPointFast(self,np.ndarray point,np.ndarray transformed):
		general_transform_point(self.transform,point.data[0],point.data[1],point.data[2],<VIO_Real*>transformed.data,<VIO_Real*>transformed.data+1,<VIO_Real*>transformed.data+2)
		
		
	def transformPoint(self,point):
		cdef np.ndarray transformed = np.zeros(3,np.float64)
		cdef np.ndarray point2 = np.array(point,np.float64)
		general_transform_point(self.transform,point2[0],point2[1],point2[2],<VIO_Real*>transformed.data,<VIO_Real*>transformed.data+1,<VIO_Real*>transformed.data+2)
		return transformed
		
	
	def getDeformation(self,invert=False,coordinateMap=False,source=None,target=None,like=None):
		cdef VIO_Volume vol
		cdef VIO_General_transform *xfm
		cdef VIO_Real vcoord[4], wcoord[3], wcoord_t[3], value
		cdef np.ndarray data,output
		cdef int x,y,z,v,invert_flag, coord_map_flag
		
#		if self.transformType == 'concatenated':
			# We need to make a (light) copy because we might have to invert some grid transforms
			# and we want to do it on the low res ones before we sample
#			transforms2 = []
#			for x in self.transforms:
#				if x.inverseFlag == invert:
#					transforms2.append(x)
#				else:
#					transforms2.append(x.inverse)
					
					
#				i.evaluateGrid()
#		else:
#			transform = self
		
		transform = self
		
		if target:
			like = target
		
		if not like:
			for i in transform.transforms:
				if i.transformType == 'grid':
					like = i.getData(noInverse=True)
					break
					
		vol = <VIO_Volume>PyCapsule_GetPointer(like.volumePtr,NULL)
		metadata = like.metadata
		shp = metadata['shape']
		if len(shp) == 3:
			shp = [3,shp[0],shp[1],shp[2]]
			spacing = [1.0,metadata['spacing'][0],metadata['spacing'][1],metadata['spacing'][2]]
			starts = [0.0,metadata['starts'][0],metadata['starts'][1],metadata['starts'][2]]
			names = ['vector_dimension',metadata['names'][0],metadata['names'][1],metadata['names'][2]]
		else:
			spacing = metadata['spacing']; starts = metadata['starts']; names = metadata['names']

		output = np.zeros(shp,np.float64)
		xfm = transform.transform
		invert_flag = invert
		coord_map_flag = coordinateMap
		
		if coordinateMap:
			# We're making a coordinate map so we need the transform including the world transforms
			transforms = transform.transforms
			transforms.insert(0,source.getVoxelToWorldTransform())
			transforms.append(target.getVoxelToWorldTransform().inverse)
			tempTransform = VIOGeneralTransform(transforms)
			xfm = tempTransform.transform

		for x in range(0,shp[-1]):
			for y in range(0,shp[-2]):
				for z in range(0,shp[-3]):
					if not coord_map_flag:					# Just make a world-space deformation map
						vcoord[0] = 0; vcoord[1] = z; vcoord[2] = y; vcoord[3] = x
						convert_voxel_to_world(vol,vcoord,wcoord,wcoord+1,wcoord+2)
						transform_or_invert_point(xfm,invert_flag,wcoord[0],wcoord[1],wcoord[2],wcoord_t,wcoord_t+1,wcoord_t+2)
						for v in range(0,3):
							value = wcoord_t[v] - wcoord[v]
							output[v,z,y,x] = value
							
					else:								# Make a voxel to voxel deformation map
						transform_or_invert_point(xfm,invert_flag,x,y,z,wcoord_t,wcoord_t+1,wcoord_t+2)
						output[0,z,y,x] = wcoord_t[2]
						output[1,z,y,x] = wcoord_t[1]
						output[2,z,y,x] = wcoord_t[0]

		newVolume = VIOVolume(output,spacing=spacing,starts=starts,names=names)
		return newVolume
		
		
	def evaluateGrid(self):
		if self.transform.type == GRID_TRANSFORM:
			if self.transform.inverse_flag:
				vol = VIOVolume(self.getData())
				self.transform.displacement_volume = <VIO_Volume>PyCapsule_GetPointer(vol.volumePtr,NULL)
				self.flipInverseFlag()
		
		
	def flipInverseFlag(self):
		self.transform.inverse_flag = not self.transform.inverse_flag
		
		
	def getData(self,noInverse=False):
		cdef np.ndarray xfm = np.zeros((4,4),np.float64)
		if self.transform.type == CONCATENATED_TRANSFORM:
			transforms = self.transforms
			return [i.getData(noInverse) for i in transforms]
		elif self.transform.type == LINEAR:
			memcpy(<char*>xfm.data,<char*>self.transform.linear_transform.m[0],16*8)
			if self.inverseFlag and not noInverse:
				xfm = xfm.I
			return xfm.T
		elif self.transform.type == GRID_TRANSFORM:
			if self.transform.inverse_flag and not noInverse:
				return self.getDeformation(invert=True)
			else:
				grid = VIOVolume();
				grid.setVolumePtr(<VIO_Volume>self.transform.displacement_volume,0)
				return grid
		
		
	property inverse:
		def __get__(self):
			cdef VIO_Transform *temp = NULL
			cdef VIO_General_transform *ntptr = NULL
			if self.transform.type == GRID_TRANSFORM:
				newTransform = VIOGeneralTransform(self.getData(noInverse=True))
				newTransform.inverseFlag = not self.inverseFlag
			elif self.transform.type == LINEAR:
				newTransform = VIOGeneralTransform(self.data)
				ntptr = <VIO_General_transform*>PyCapsule_GetPointer(newTransform.transformPtr,NULL)
				temp = ntptr.linear_transform
				ntptr.linear_transform = ntptr.inverse_linear_transform
				ntptr.inverse_linear_transform = temp
			elif self.transform.type == CONCATENATED_TRANSFORM:
				xfms = [i.inverse for i in self.transforms]
				newTransform = VIOGeneralTransform(xfms)
			else:
				print "I don't know how to handle this kind of transform"
				newTransform = None
			return newTransform
	
	
	property inverseFlag:
		def __get__(self):
			return self.transform.inverse_flag

		def __set__(self,flag):
			self.transform.inverse_flag = flag
		
	property transformType:
		def __get__(self):
			if self.transform.type == LINEAR:
				return 'linear'
			elif self.transform.type == GRID_TRANSFORM:
				return 'grid'
			elif self.transform.type == CONCATENATED_TRANSFORM:
				return 'concatenated'
			else:
				return self.transform.type
		
		
	property transformN:
		def __get__(self):
			return get_n_concated_transforms(self.transform)
	
	
	property transforms:
		def __get__(self):
			if self.transform.type == CONCATENATED_TRANSFORM:	
				transforms = []
				for i in range(0,self.transform.n_transforms):
					xfm = VIOGeneralTransform(); xfm.setTransformPtr(&self.transform.transforms[i])
					transforms.append(xfm)
			else:
				return [self]
			return transforms
			
			
	property data:
		def __get__(self):
			return self.getData()
			

	property transformPtr:
		def __get__(self):
			return PyCapsule_New(<void *>self.transform, NULL, NULL)
		
		def __set__(self,ptr):
			self.transform = <VIO_General_transform *>PyCapsule_GetPointer(ptr,NULL)

	

