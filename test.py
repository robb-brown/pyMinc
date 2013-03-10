#!/usr/local/bin/python

from pyMinc import *
from pyMinctracc import *
from pylab import *
from time import time
from scipy.ndimage.interpolation import *
from scipy.ndimage import *
import lingo.mincUtils as minc



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





sourceF = '/temp/CIS-B_001-HSC-1_pt_00001_v02_t1g.mnc.gz'
targetF = '/temp/icbm152_sym_t1w.mnc.gz'

source = VIOVolume(sourceF,type=float32)
target = VIOVolume(targetF,type=float32)

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
	linear = m.minctracc(source4,target4,initialXFM=None,transformType='lsq12',debug=False)
	t2 = time()
	print 'Doing non linear 16 registration'
	nonlinear16 = m.minctracc(source4,target4,initialXFM=linear,transformType='nonlinear',step=[16.0,16.0,16.0],debug=False)
	t3 = time()
	print 'Doing non linear 8 registration'
	nonlinear8 = m.minctracc(source8,target8,initialXFM=nonlinear16,transformType='nonlinear',step=[8.0,8.0,8.0],debug=False)
	t4 = time()
	print 'Doing non linear 4 registration'
	nonlinear4 = m.minctracc(source4,target4,initialXFM=nonlinear8,transformType='nonlinear',step=[4.0,4.0,4.0],debug=False)
	t5 = time()

	print "Registration took: "
	print "	Linear: %0.2f s" % (t2-t1)
	print "	Nonlinear 16: %0.2f s" % (t3-t2)
	print "	Nonlinear 8: %0.2f s" % (t4-t3)
	print "	Nonlinear 4: %0.2f s" % (t5-t4)
	
	transform = nonlinear4
	
	transform.write('/temp/nonlinear.xfm')
	
else:
	transform = VIOGeneralTransform('/temp/nonlinear.xfm')


xfm = mat(transform.data[0])

param = minc.xfmParameters(xfm); param.pop('shears')
xfm = mat(minc.xfmFromParameters(**param))

# Make a mincresample image
minc.writeXFM(array(xfm),'/temp/temp.xfm')
cmd = 'mincresample %s %s -transformation %s -like %s' % (sourceF,'/temp/sourceResampled.mnc.gz','/temp/temp.xfm',targetF)
minc.execute(cmd)
mincResampled = VIOVolume('/temp/sourceResampled.mnc.gz',type=float32).data

# Resample source
linearResampled = linearResample(source,transform=xfm,like=target)

# Nonlinear resampling
xfm = transform.transforms[1]

deformation = xfm.getDeformation(invert=True,coordinateMap=False)
metadata = deformation.metadata; spacing = metadata['spacing'][1:]; starts = metadata['starts'][1:]; names = metadata['names'][1:]
planes = [VIOVolume(deformation.data[i],spacing=spacing,starts=starts,names=names) for i in range(0,3)]
big = array([linearResample(i,transform=identity(4),like=target) for i in planes])
slc = [slice(0,i) for i in shape(big)[1:]]
cMap = big/4. + mgrid[slc]

figure(5); clf(); imshow(deformation.data[0,:,34]); colorbar()
figure(6); clf(); imshow(big[0,:,34*4]); colorbar()
figure(7); clf(); imshow(cMap[0,:,34*4]); colorbar()


nonlinearResampled = map_coordinates(linearResampled,cMap,order=1)
nonlinearResampledM = VIOVolume('/temp/sourceToTargetNonlinear.mnc')


s = 75
titles = ['Source','Target','Resampled Source (Linear)','Resampled Source (nonlinear)','Minc resampled source (nonlinear)]
for i,data in enumerate([source.data,target.data,linearResampled,nonlinearResampled,nonlinearResampledM]):
	figure(i+1);
	subplot(221); imshow(data[s]);title(titles[i])
	subplot(222); imshow(data[:,s]);
	subplot(223); imshow(data[:,:,s]);

