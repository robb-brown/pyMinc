from scipy.ndimage.interpolation import *
from scipy.ndimage import *
from numpy import *
import lingo.mincUtils as minc
from pyMinc import VIOVolume, VIOGeneralTransform


def blurImage(image,level):
	size = level / 2 / 2.35
	return VIOVolume(filters.gaussian_filter(image.data,sigma=size,mode='constant'),**image.metadata)
	

def linearResample(source,transform,like,invert=False,order=2):
	if transform.__class__ == VIOGeneralTransform:
		xfm = mat(transform.data)
	else:
		xfm = mat(transform)
	if invert:
		xfm = xfm.I
	vol = source; target = like;
	metadata = vol.metadata
	im = transpose(vol.data,metadata['spatialAxes'])			# Transforms are stored as x,y,z so our image should be too
	spacing = array(metadata['spacing'])[metadata['spatialAxes']]
	starts = array(metadata['starts'])[metadata['spatialAxes']]
	extents = array(metadata['shape'])[metadata['spatialAxes']]

	# Get the (flipped) source to world transform
	sourceToWorld = mat(minc.xfmFromParameters(translations=starts,scales=spacing))
	
	# Get the target to world transform.  Note that this is assuming the atlas, which doesn't need to be flipped
	targetToWorld = mat(target.voxelToWorldTransform.data)
	
	# Total source voxel to target voxel transform, inverted
	# It's important to multiply backwards, and don't forget about order of operations!
	total = (targetToWorld.I*(xfm*sourceToWorld)).I
	
	# These could be useful for doing transform input sampling
	e1 = (total*mat([0,0,0,1]).T).T
	e2 = (total*mat(list(array(source.metadata['shape'])[source.metadata['spatialAxes']])+[1]).T).T
	s1 = (total*mat([1,1,1,1]).T).T - e1
	#outputShape = (array(maximum(e1,e2)+1)[0,0:3]).astype(int16).tolist()   
	
	# Use this if there's a like (using the target for this one)
	outputShape = array(target.metadata['shape'])[target.metadata['spatialAxes']]
	
	# Do the actual (linear) transform
	linearResampled = affine_transform(im,total[0:3,0:3],offset=array(total)[0:3,3],output_shape=outputShape,order=order)
	
	# Put the transformed axes back in the MINC order.  Note that we should check for negative spacing and flip here too
	linearResampled = transpose(linearResampled,target.metadata['spatialAxes'])
	
	return VIOVolume(linearResampled,**like.metadata)
	
	
def nonlinearResample(source,transform,like,invert=False,order=2):
	if invert:
		xfm = transform.inverse
	else:
		xfm = transform
	cMap = xfm.getDeformation(invert=True,coordinateMap=True,source=source,target=like).data
	nonlinearResampled = map_coordinates(source.data,cMap,order=order)
	return VIOVolume(nonlinearResampled,**like.metadata)
	


def resample(source,transform,like,invert=False,order=2):
	allLinear = len([1 for i in transform.transforms if i.transformType == 'grid']) == 0
	if allLinear:
		return linearResample(source,transform,like,invert,order=order)
	else:
		return nonlinearResample(source,transform,like,invert,order=order)
	
	
	