#!/usr/local/bin/python

from pyMinc import *
from pyMinctracc import *
from pylab import *
from time import time
from scipy.ndimage.interpolation import *
from scipy.ndimage import *
import lingo.mincUtils as minc
from numpy import ma as ma
import os.path as p



def linearResample(source,transform,like):
	vol = source; target = like; xfm = transform
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
	linearResampled = affine_transform(im,total[0:3,0:3],offset=array(total)[0:3,3],output_shape=outputShape,order=1)
	
	# Put the transformed axes back in the MINC order.  Note that we should check for negative spacing and flip here too
	linearResampled = transpose(linearResampled,target.metadata['spatialAxes'])
	
	return linearResampled





sourceF = '/temp/source.mnc.gz'
targetF = '/temp/target.mnc.gz'

source = VIOVolume(sourceF,type=float32)
target = VIOVolume(targetF,type=float32)
sourceMask = VIOVolume('/temp/sourceMask.mnc.gz',type=float32)
targetMask = VIOVolume('/temp/targetMask.mnc.gz',type=float32)

m = Minctracc()

if 0:			# Time comparison
	print "Simplex"
	t1 = time()
	xfm = m.minctracc(source,target,initialXFM=None,transformType='lsq12',debug=False,verbose=0)
	t2 = time()

	print "\n\nBFGS"
	t3 = time()
	xfm = m.minctracc(source,target,initialXFM=None,transformType='lsq12',optimizer='BFGS',debug=False,verbose=0)
	t4 = time()
	print "\nSimplex took %0.2fs.  BFGS took %0.2fs." % (t2-t1,t4-t3)


if 0:
	sourceMeta = source.metadata
	targetMeta = target.metadata
	sourceArgs = {'spacing':sourceMeta['spacing'],'starts':sourceMeta['starts']}
	targetArgs = {'spacing':targetMeta['spacing'],'starts':targetMeta['starts']}

	source16 = VIOVolume(filters.uniform_filter(source.data,size=16,mode='constant'),**sourceArgs)
	source8 = VIOVolume(filters.uniform_filter(source.data,size=8,mode='constant'),**sourceArgs)
	source4 = VIOVolume(filters.uniform_filter(source.data,size=4,mode='constant'),**sourceArgs)
	target16 = VIOVolume(filters.uniform_filter(target.data,size=16,mode='constant'),**targetArgs)
	target8 = VIOVolume(filters.uniform_filter(target.data,size=8,mode='constant'),**targetArgs)
	target4 = VIOVolume(filters.uniform_filter(target.data,size=4,mode='constant'),**targetArgs)

	t1 = time()
	print 'Doing linear registration'
	linear = m.minctracc(source4,target4,sourceMask=sourceMask,targetMask=targetMask,initialXFM=None,transformType='lsq12',debug=False)
	linearI = linear.inverse
	t2 = time()
	print 'Doing non linear 16 registration'
	nonlinear16 = m.minctracc(target4,source4,sourceMask=targetMask,targetMask=sourceMask,initialXFM=linearI,transformType='nonlinear',step=[16.0,16.0,16.0],debug=False)
	t3 = time()
	print 'Doing non linear 8 registration'
	nonlinear8 = m.minctracc(target8,source8,sourceMask=targetMask,targetMask=sourceMask,initialXFM=nonlinear16,transformType='nonlinear',step=[8.0,8.0,8.0],debug=False)
	t4 = time()
	print 'Doing non linear 4 registration'
	nonlinear4 = m.minctracc(target4,source4,sourceMask=targetMask,targetMask=sourceMask,initialXFM=nonlinear8,transformType='nonlinear',step=[4.0,4.0,4.0],debug=False)
	t5 = time()

	nonlinear4 = nonlinear4.transforms[1];
	nonlinear4 = VIOGeneralTransform([linear,nonlinear4])
	nonlinear4.transforms[1].flipInverseFlag()

	print "Registration took: "
	print "	Linear: %0.2f s" % (t2-t1)
	print "	Nonlinear 16: %0.2f s" % (t3-t2)
	print "	Nonlinear 8: %0.2f s" % (t4-t3)
	print "	Nonlinear 4: %0.2f s" % (t5-t4)
	
	transform = nonlinear4
	
	transform.write('/temp/nonlinear.xfm')
	linear.write('/temp/linear.xfm')
	
else:
	transform = VIOGeneralTransform('/temp/nonlinear.xfm')


xfm = mat(transform.transforms[0].data)

param = minc.xfmParameters(xfm); param.pop('shears')
xfm = mat(minc.xfmFromParameters(**param))

# Make a mincresample image
if not p.exists('/temp/sourceResampled.mnc.gz'):
	cmd = 'mincresample %s %s -transformation %s -like %s' % (sourceF,'/temp/sourceResampled.mnc.gz','/temp/nonlinear.xfm',targetF)
	minc.execute(cmd,quiet=False)

mincResampled = VIOVolume('/temp/sourceResampled.mnc.gz',type=float32)

# Resample source
linearResampled = linearResample(source,transform=xfm,like=target)

# Nonlinear resampling
print "Doing nonlinear resampling"
xfm = transform.transforms[1]

if 0:		# Alternate (faster?) way of making the full size deformation field
	t6 = time()
	deformation = xfm.getDeformation(invert=False,coordinateMap=False)
	metadata = deformation.metadata; spacing = metadata['spacing'][1:]; starts = metadata['starts'][1:]; names = metadata['names'][1:]
	planes = [VIOVolume(deformation.data[i],spacing=spacing,starts=starts,names=names) for i in range(0,3)]
	big = array([linearResample(i,transform=identity(4),like=target) for i in planes])
	slc = [slice(0,i) for i in shape(big)[1:]]
	cMap = big + mgrid[slc]
	t7 = time()
	t8 = time()
	nonlinearResampled = map_coordinates(linearResampled,cMap,order=1)
	t9 = time()
	
	figure(6); clf(); imshow(deformation.data[0,:,34]); colorbar()
	figure(7); clf(); imshow(big[0,:,34*4]); colorbar()
	figure(8); clf(); imshow(cMap[0,:,34*4]); colorbar()
	
if 1:			# Same way minc does it?
	t6 = time()
	cMap = transform.getDeformation(invert=False,coordinateMap=True,like=target).data
	t7 = time()
	t8 = time()
	nonlinearResampled = map_coordinates(source,cMap,order=1)
	t9 = time()

	deformation = xfm.getDeformation(invert=False,coordinateMap=False)
	slc = [slice(0,i) for i in shape(cMap)[1:]]
	big = cMap - mgrid[slc]	
	figure(6); clf(); imshow(deformation.data[0,:,34]); colorbar()
	figure(7); clf(); imshow(big[0,:,34*4]); colorbar()
	figure(8); clf(); imshow(cMap[0,:,34*4]); colorbar()


print "Time to do nonlinear resample: %0.2f s to get cMap, %0.2f s to do resample, %0.2f s total" % (t7-t6,t9-t8,t7-t6+t9-t8)

titles = ['Source','Target','Resampled Source (Linear)','Resampled Source (nonlinear)','Minc resampled source (nonlinear)']
sdata = ma.array(source.data); sdata.mask = where(sourceMask.data > 0,0,1)
tdata = ma.array(target.data); tdata.mask = where(targetMask.data > 0,0,1)
for i,data in enumerate([sdata,tdata,linearResampled,nonlinearResampled,mincResampled.data]):
	s = array(shape(data)) / 2
	figure(i+1);
	subplot(221); imshow(data[s[0]]);title(titles[i])
	subplot(222); imshow(data[:,s[1]]);
	subplot(223); imshow(data[:,:,s[2]]);

