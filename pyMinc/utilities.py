from scipy.ndimage.interpolation import *
from scipy.ndimage import *
from numpy import *
import mincUtils as minc
from pyMinc import VIOVolume, VIOGeneralTransform


def blurImage(image,level):
	size = level / 2. / 2.35			# sigma in mm
	spacing = array(image.metadata['spacing'])
	sizes = size / spacing
	return VIOVolume(filters.gaussian_filter(image.data,sigma=sizes,mode='constant'),**image.metadata)
	

def linearResample(source,transform,like,invert=False,order=3,originalSpacing=False):
	if transform.__class__ == VIOGeneralTransform:
		xfm = mat(transform.data)
	else:
		xfm = mat(transform)
		
	if alltrue(xfm == identity(4)):
		return VIOVolume(source)
		
	if invert:
		xfm = xfm.I
	vol = source; target = like;
	metadata = vol.metadata
	im = transpose(vol.data,metadata['spatialAxes'])			# Transforms are stored as x,y,z so our image should be too
	spacing = array(metadata['spacing'])[metadata['spatialAxes']]
	starts = array(metadata['starts'])[metadata['spatialAxes']]
	extents = array(metadata['shape'])[metadata['spatialAxes']]

	# Get the source to world transform
	sourceToWorld = mat(source.voxelToWorldTransform.data)
	
	# Get the target to world transform.
	targetToWorld = mat(target.voxelToWorldTransform.data)
	
	# Idea here is to change the targetToWorld transform so the output spacing is equal to the source spacing.
	# Might have to change the starts too
	# if originalSpacing:
	# 	params = minc.xfmParameters(targetToWorld)
	# 	params['scales'] = spacing
	# 	targetToWorld = minc.xfmFromParameters(**params)
	
	# Total source voxel to target voxel transform, inverted
	# It's important to multiply backwards, and don't forget about order of operations!
	total = (targetToWorld.I*(xfm*sourceToWorld)).I
	
	#total = targetToWorld*xfm.I*sourceToWorld
	
	# These could be useful for doing transform input sampling
	e1 = (total*mat([0,0,0,1]).T).T
	e2 = (total*mat(list(array(source.metadata['shape'])[source.metadata['spatialAxes']])+[1]).T).T
	s1 = (total*mat([1,1,1,1]).T).T - e1
	#outputShape = (array(maximum(e1,e2)+1)[0,0:3]).astype(int16).tolist()   
	
	# Get the shape for our output image
	outputShape = array(target.metadata['shape'])[target.metadata['spatialAxes']]
	
	# Do the actual (linear) transform
	linearResampled = affine_transform(im,total[0:3,0:3],offset=array(total)[0:3,3],output_shape=outputShape,order=order)
	
	# Put the transformed axes back in the MINC order.  It doesn't appear to be necessary to flip if spacing is negative.
	linearResampled = transpose(linearResampled,target.metadata['spatialAxes'])
	
	metadata = dict(like.metadata)
	if originalSpacing:
		metadata['spacing'] = spacing[metadata['spatialAxes']]
	
	return VIOVolume(linearResampled,**metadata)
	
	
def nonlinearResample(source,transform,like,invert=False,order=2):
	if invert:
		xfm = transform.inverse
	else:
		xfm = transform
	cMap = xfm.getDeformation(invert=True,coordinateMap=True,source=source,target=like).data
	nonlinearResampled = map_coordinates(source.data,cMap,order=order)
	return VIOVolume(nonlinearResampled,**like.metadata)
	


def resample(source,transform,like,invert=False,order=3):
	if order.__class__ == str and order.startswith('nearest'):
		order = 0
	allLinear = len([1 for i in transform.transforms if i.transformType == 'grid']) == 0
	if allLinear:
		return linearResample(source,transform,like,invert,order=order)
	else:
		return nonlinearResample(source,transform,like,invert,order=order)
	
	
	